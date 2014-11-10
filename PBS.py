import subprocess
import numpy
import xml.etree.ElementTree as ET
import re
import time

def submit(string):
    pipe = subprocess.Popen(['qsub'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = pipe.communicate(string)[0]
    match = re.match('([0-9]*)\..*', stdout)
    if pipe.returncode or not match:
        raise Exception("qsub failed: %s", stdout)
    return match.group(1)

def status(jobid):
    """ returns R, Q, E, C, or U(for unknown, eg jobid is not in qstat"""
    try:
        xml = subprocess.check_output(['qstat', '-x', str(jobid)])
        tree = ET.fromstring(xml)
        ele = tree.find('Job/job_state')
        return ele.text
    except subprocess.CalledProcessError:
        return 'U'

def delete(jobid):
    return subprocess.check_call(['qdel', str(jobid)])

def wait(jobid):
    if not isinstance(jobid, (list, tuple, set)):
        while status(jobid) in 'RQ':
            time.sleep(10)
    else:
        for job in jobid:
            wait(job)
