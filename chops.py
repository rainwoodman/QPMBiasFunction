import os
import numpy
import QPM
import glob
import sharedmem
edges = numpy.loadtxt('100bins.txt')

def processfile(file, out, me):
    """compute the cross correlation from bin 0 to me(inc),
    save to out-me"""
    if out.endswith('.npz'):
        out = out[:-4]
    out = out + ('-%03d' % me) + '.npz'
    if os.path.exists(out):
        return
    f = QPM.QPMSubsampleFile(file)
    xis = []
    mybf = QPM.LogWindowBiasFunction(edges[me], edges[me + 1])
    mymock = QPM.halomock(f, 1e6, mybf, continuous=False)
    for i in range(me + 1):
        bf = QPM.LogWindowBiasFunction(edges[i], edges[i + 1])
        mock = QPM.halomock(f, 1e6, bf, continuous=False)
        # the error is nonsense since we have just one sample
        # deal with it in combine
        xi = QPM.xi([mock], mocks2=[mymock])
        r = xi[0]
        xis.append(xi[1])
        print edges[i], edges[i + 1], len(mock[0]), mock[1].sum()
    numpy.savez(out, r=r, xi=xis, edges=edges, N=1e6, me=me)

def combine(prefix, me=None):
    if me is None:
        d = {}
        with sharedmem.MapReduce() as pool:
            def work(i):
                dd = combine(prefix, i)
                return dd
            def reduce(dd):
                d.update(dd)
            pool.map(work, range(0, len(edges) - 1), reduce=reduce)
        return d
    files = glob.glob(prefix + '/*-%03d.npz' % me)
    l = []
    for f in files:
        f = numpy.load(f)
        xi = f['xi']
        m = f['me']
        if m != me: 
            print f, 'is bad'
            continue
        l.append(xi)

    r = f['r']
    mean = numpy.mean(l, axis=0)
    std = numpy.std(l, axis=0) * len(l) ** -0.5
    d = {}
    for i in range(me + 1):
        a = numpy.array([r, mean[i], std[i]])
        d[i, me] = a
        d[me, i] = a
    return d

def combine1(prefix):
    """ combine all xi's measured in prefix """
    files = glob.glob(prefix + '/*.npz')
    files = [numpy.load(f) for f in files]
    edges = files[0]['edges']
    r = []
    m = []
    e = []
    r = files[0]['xi'][0, 0]
    l = []
    for i in range(len(edges) - 1):
        xi = [f['xi'][i, 1, :] for f in files]
        l.append(
                (r, 
                numpy.mean(xi, axis=0),
                numpy.std(xi, axis=0) * len(xi) ** -0.5))
    return numpy.array(l)

def main():
    processfile(sys.argv[1], sys.argv[2], int(sys.argv[3]))

if __name__ == "__main__":
    import sys
    main()


