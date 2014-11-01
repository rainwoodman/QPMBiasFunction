import subprocess
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
    xml = subprocess.check_output(['qstat', '-x', str(jobid)])
    tree = ET.fromstring(xml)
    ele = tree.find('Job/job_state')
    return ele.text

def delete(jobid):
    return subprocess.check_call(['qdel', str(jobid)])

def wait(jobid):
    while 'C' != status(jobid):
        time.sleep(10)
    
