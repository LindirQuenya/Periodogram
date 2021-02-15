from testshared import basictest

testname = "AblsIplav"
cmdline = "-a bls -b 150 -q 0.01 -Q 0.1 -i plav -p 0.1843 -P 1208.8544 -u 0.1 -N 50 -S 1.0"

def test_AblsIplav_lc1():
    basictest(testname, cmdline)
