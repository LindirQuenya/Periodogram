import os
import filecmp
from os.path import dirname

def basictest(testname, cmdline):
    testdir = dirname(dirname(__file__))
    datapath = os.path.abspath(testdir + "/datasets/lc1_python.txt")
    outpath = os.path.abspath(testdir + "/testresults/"+testname)
    goodpath = os.path.abspath(testdir + "/goodresults/"+testname)
    os.system("python -m periodogram "+cmdline+" "+datapath+" "+outpath+"/lc1_python.txt.out")
    match, mismatch, errors = filecmp.cmpfiles(outpath, goodpath, ['lc1_python.txt.out', 'lc1_python.txt.out.top'])
    if "Als" in testname:
        fuzzy_cmpfiles(goodpath+'/lc1_python.txt.out', outpath+'/lc1_python.txt.out')
        fuzzy_cmpfiles(goodpath+'/lc1_python.txt.out.top', outpath+'/lc1_python.txt.out.top')
    else:
        assert mismatch==[] and errors==[]


# Assumes unlabeled output.
def fuzzy_cmpfiles(expected, results):
    exp = open(expected, 'r')
    res = open(results, 'r')
    e_cont = exp.readlines()
    r_cont = res.readlines()
    exp.close()
    res.close()
    e_per, e_pow = [0]*len(e_cont), [0]*len(e_cont)
    r_per, r_pow = [0]*len(r_cont), [0]*len(r_cont)
    assert len(e_cont) == len(r_cont)
    for i in range(len(e_cont)):
        e_per[i], e_pow[i] = [float(k) for k in e_cont[i].split(",")]
        r_per[i], r_pow[i] = [float(k) for k in r_cont[i].split(",")]
    diff_per = max([abs(e_i - r_i) for e_i, r_i in zip(e_per, r_per)])
    diff_pow = max([abs(e_i - r_i) for e_i, r_i in zip(e_pow, r_pow)])
    assert diff_per<0.0000000001 and diff_pow<0.0000000001

