from testshared import basictest

testname = "AblsIexp"
cmdline = "-a bls -b 150 -q 0.01 -Q 0.1 -i exp -p 0.1843 -P 1208.8544 -o 10 -N 50 -S 1.0"

def test_AblsIexp_lc1():
    basictest(testname, cmdline)
