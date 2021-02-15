import os
import filecmp
from os.path import dirname

testname = "AlsIfixedf"
cmdline = "-a ls -i fixedf -p 0.1843 -P 1208.8544 -d 0.0011007 -N 50 -S 1.0"

def test_AlsIfixedf_lc1():
    testdir = dirname(dirname(__file__))
    datapath = os.path.abspath(testdir + "/datasets/lc1_python.txt")
    outpath = os.path.abspath(testdir + "/testresults/"+testname)
    goodpath = os.path.abspath(testdir + "/goodresults/"+testname)
    os.system("python -m periodogram "+cmdline+" "+datapath+" "+outpath+"/lc1_python.txt.out")
    match, mismatch, errors = filecmp.cmpfiles(outpath, goodpath, ['lc1_python.txt.out', 'lc1_python.txt.out.top'])
    assert mismatch==[] and errors==[]
