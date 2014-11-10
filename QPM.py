import numpy
import sharedmem
import hashlib
import os.path
import pickle
from kdcount import correlate

class SubSampleStore(object):
    def __init__(self, template):
        """ filename template
        """
        self.template = template
        rem = os.path.dirname(template)
        u = os.path.basename(template)
        while not '%' in u:
            basename = os.path.basename(rem)
            nrem = os.path.dirname(rem)
            if nrem == rem:
                raise Exception("Filename shall have a %04d pattern somewhere!")
            rem = nrem
            u = '_'.join([basename, u])
        self.uniquename = u
    def query(self, *args):
        try:
            return SubSampleFile(self.template % args)
        except IOError:
            return SubSampleFiles(self.template % args)
    def __reduce__(self):
        return (SubSampleStore, (self.template,))

class Mock(object):
    @staticmethod
    def gethash(file, factory, parameters, N):
        h = hashlib.sha1()
        N = numpy.array([N], dtype='int32')
        h.update(N)
        h.update(os.path.basename(file))
        h.update(factory)
        h.update(numpy.array(parameters))
        return h.hexdigest()

    def __init__(self, sampleid, factory, parameters, N):
        self.sampleid = sampleid
        if hasattr(factory, '__call__'):
            pass
        else:
            factory = globals()[factory]

        self.factory = factory
        self.parameters = parameters
        self.biasfunction = factory(*parameters)
        self.N = N
#        self.hash = self.gethash(sample, factory.__name__, parameters, N)    
        self.hash = os.path.join(factory.__name__, '%.1E' % N, *['P%+0.3f' % p for p in
            parameters] + ['%04d' % sampleid])
    def __reduce__(self):
        return (self.__class__, (self.sampleid, self.factory.__name__, self.parameters, self.N))

class XiFactory(object):
    def __init__(self, samplestore, path):
        self.samplestore = samplestore
        self.path = path
        try:
            os.makedirs(path)
        except OSError:
            pass
        with open(os.path.join(self.path, 'store.dat'), 'wb') as outfile:
            pickle.dump(self, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    def __reduce__(self):
        return (XiFactory, (self.samplestore, self.path))

    def getfile(self, mock):
        hash = mock.hash
        out = os.path.join(self.path, hash + '.dat')
        return out

    def update(self, mock):
        out = self.getfile(mock)
        if os.path.exists(out):
            return
        self.replace(mock)

    def replace(self, mock):
        out = self.getfile(mock)
        f = self.samplestore.query(mock.sampleid)
        mymock = halomock(f, mock.N, mock.biasfunction, continuous=False)
        r, xi_, junk = xi([mymock])
        try:
            os.makedirs(os.path.dirname(out))
        except OSError:
            pass
        with open(out, 'wb') as outfile:
            data = {}
            data['r'] = r
            data['xi'] = xi_
            data['mock'] = mock
            pickle.dump(data, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    def query(self, mock):
        """ fails of mock is not calculated """
        out = self.getfile(mock)
        with open(out, 'rb') as outfile:
            data = pickle.load(outfile)
        return data

    def xi(self, mocks):
        result = [self.query(mock) for mock in mocks]
        xi = numpy.empty((3, len(result[0]['r'])))
        xi[0] = result[0]['r']
        xis = [i['xi'] for i in result]
        xi[1] = numpy.mean(xis, axis=0)
        xi[2] = numpy.std(xis, axis=0) * len(xis) ** -0.5
        return xi

class SubSampleFiles(list):
    def __init__(self, filename):
        i = 0
        self.filename = filename
        while True:
            fn = filename + ('.%02d' % i)
            print fn
            if not os.path.exists(fn): break
            i = i + 1
            self.append(SubSampleFile(fn))
        if len(self) == 0:
            raise Exception("No files found at %s", filename)

class SubSampleFile(object):
    def __init__(self, filename):
        self.ints = numpy.memmap(filename, mode='r', dtype='i4')
        self.eflag = self.ints[0]
        self.hsize = self.ints[1]

        self.npart = self.ints[2]
        self.nsph = self.ints[3]
        self.nstar = self.ints[4]
        self.aa = self.ints[5].view('f4')
        self.softlen = self.ints[6].view('f4')
        npart = self.npart
        self.pos = self.ints[7:7+npart * 3].view('f4').reshape(-1, 3)
        self.vel = self.ints[7+npart*3:7+npart*6].view('f4').reshape(-1, 3)
        self.rho = self.ints[7+npart*6:7+npart*7].view('f4').reshape(-1)
        self.filename = filename
    def __str__(self):
        return '\n'.join([
                'QPMFile: %s' % self.filename,
                'npart: %d' % self.npart,
                'nsph: %d' % self.nsph,
                'nstar: %d' % self.nstar,
                'aa: %g' % self.aa,
                'softlen: %g' % self.softlen,
                'pos: %s' % str(self.pos),
                'vel: %s' % str(self.vel),
                'rho: %s' % str(self.rho)])
"""
"""
def LogNormalBiasFunction(mu, sigma):
    def func(rho):
        lnrho = numpy.log(rho)
        lnrho -= mu
        fac = -1.0 / (2 * sigma ** 2)
        return numpy.exp(lnrho ** 2 * fac)
    return func
def DumbSkewedLogNormalBiasFunction(mu, sigma, epsilon):
    def func(rho):
        lnrho = numpy.log(rho)
        lnrho -= mu
        lnrho /= sigma
        fac = -1.0 / (2)
        t = numpy.exp(lnrho ** 2 * fac) * ( 1 + epsilon * lnrho)
        t[t <= 0] = 0.0
        return t
    return func
def LogWindowBiasFunction(lnrho1, lnrho2):
    def func(rho):
        lnrho = numpy.log(rho)
        mask = lnrho >= lnrho1
        mask &= lnrho < lnrho2
        return mask * 1.0
    return func

def SkewedLogNormalBiasFunction(mu, sigma, alpha):
    def func(rho):
        x = numpy.log(rho)
        x -= mu
        x /= sigma
        fac = - 0.5 * (1 + numpy.erf(alpha / 2. ** 0.5 * x))
        return numpy.exp(x ** 2 * fac)
    return func
def ASymPolyExpBiasFunction(mu, leftpol, rightpol):
    """Asymmetric polynomial exponential bias fuction is of the
    form:
        leftpol and right pol are rel to mu
        exp(leftpol) if lnrho < mu
        exp(rightpol) if lnrho >= mu
    """
    def func(rho):
        lnrho = numpy.log(rho)
        left = lnrho < mu
        right = lnrho >= mu
        result = numpy.empty_like(lnrho)
        result[left] = numpy.polyval(leftpol, lnrho[left] - mu)
        result[right] = numpy.polyval(rightpol, lnrho[right] - mu)
        return numpy.exp(result)
    return func
def RhoClipFunction(thresh):
    thresh *= 1.0
    def func(rho):
        return thresh ** 2 * rho / (thresh ** 2 + rho)
    return func

def halomock(subsample, Nbar, biasfunction, rhoclipfunction=RhoClipFunction(1000.),
        seed=None, continuous=False):
    """ create qpm halo mocks 
        if seed is None, use /dev/random (numpy semantics)
    """
    prob = sharedmem.empty_like(subsample.rho)
    probsum = [0]
    chunksize = 1024 * 1024
    R = range(0, len(subsample.rho), chunksize)
    rng = numpy.random.RandomState(seed)
    rngs = [numpy.random.RandomState(rng.randint(1000000000)) for i in R]
    def work1(i):
        s = slice(i, i + chunksize)
        rhoclip = rhoclipfunction(subsample.rho[s])
        prob[s] = biasfunction(rhoclip)
        bad = numpy.isnan(prob[s])
        prob[s][bad] = 0
        return prob[s].sum(dtype='f8')
    def reduce1(sum):
        probsum[0] = probsum[0] + sum

    def work2(i, rng):
        s = slice(i, i + chunksize)
        Nexp = Nbar / probsum[0] * prob[s]
        if continuous:
            cut = int(len(prob[s]) * 1.0 * Nbar / len(subsample.pos))
            a = numpy.random.choice(len(prob[s]), cut, replace=False)
            pos = numpy.array(subsample.pos[s][a])
            wt = numpy.array(Nexp[a]) / cut * len(prob[s])
        else:
            Nhalo = rng.poisson(Nexp)
            mask = Nhalo > 0
            pos = numpy.array(subsample.pos[s][mask])
            wt = numpy.array(Nhalo[mask])
        return pos, wt

    with sharedmem.MapReduce() as pool:
        pool.map(work1, R, reduce=reduce1)
        pos, wt = zip(*pool.map(work2, zip(R, rngs), star=True))

    return numpy.concatenate(pos, axis=0), numpy.concatenate(wt, axis=0)

def xi(mocks, rmax=0.1, Nbins=20, mocks2=None):
    """ cross correlation mocks against mocks2, one by one
        returns (r, xi, err) """
    l = []
    for i in range(len(mocks)):
        pos1, wt1 = mocks[i]
        p1 = correlate.points(pos1, wt1, boxsize=1.0)
        if mocks2 is not None:
            pos2, wt2 = mocks2[i]
            p2 = correlate.points(pos2, wt2, boxsize=1.0)
        else:
            wt2 = wt1
            p2 = p1
        r = correlate.paircount(p1, p2, correlate.RBinning(rmax, Nbins))
        DD = r.sum1
        RR = 4 / 3. * numpy.pi * numpy.diff(r.edges ** 3) * wt1.sum() * wt2.sum()
        l.append(DD / RR - 1.0)
       
    return (r.centers, 
            numpy.mean(l, axis=0), 
            numpy.std(l, axis=0) * len(l) ** -0.5
            )
def plotfriendly(xi, boxsize):
    r = xi[0] * boxsize
    return r, r ** 2 * xi[1], r **2 * xi[2]
def MWSim(filename, boxsize):
    """ normalize r in MWSim correlation functions to given boxsize 
    """
    MWSim = numpy.loadtxt(filename, usecols=(0, 1, 2), unpack=True)
    MWSim[0] /= boxsize
    return MWSim

def matchdiff(MWSim, xi):
    match = numpy.interp(xi[0], MWSim[0], MWSim[1])
    return (xi[0], (xi[1] - match), xi[2])

def matchratio(MWSim, xi):
    match = numpy.interp(xi[0], MWSim[0], MWSim[1])
    err = numpy.interp(xi[0], MWSim[0], MWSim[2])
    return (xi[0], xi[1] / match, 1.414 * ((xi[2] / match)+ err * numpy.abs(xi[1]) / match ** 2))

def fitratioslope(MWSim, xi):
    xi = matchratio(MWSim, xi)
    s = slice(2, -5)
    p = numpy.polyfit(xi[0][s], xi[1][s], deg=1, w=xi[2][s]**-1) 
    return p

if __name__ == '__main__':
    import sys
    print QPMFile(sys.argv[1])
