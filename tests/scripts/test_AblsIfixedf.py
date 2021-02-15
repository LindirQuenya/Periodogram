from testshared import basictest

testname = "AblsIfixedf"
cmdline = "-a bls -b 150 -q 0.01 -Q 0.1 -i fixedf -p 0.1843 -P 1208.8544 -d 0.0011007 -N 50 -S 1.0"

def test_AblsIfixedf_lc1():
    basictest(testname, cmdline)
