from testshared import basictest

testname = "AblsIfixedp"
cmdline = "-a bls -b 150 -q 0.01 -Q 0.1 -i fixedp -p 0.1843 -P 1208.8544 -d 0.2452161 -N 50 -S 1.0"

def test_AblsIfixedp_lc1():
    basictest(testname, cmdline)
