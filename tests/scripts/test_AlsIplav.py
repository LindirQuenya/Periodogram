from testshared import basictest

testname = "AlsIplav"
cmdline = "-a ls -i plav -p 0.1843 -P 1208.8544 -u 0.1 -N 50 -S 1.0"

def test_AlsIplav_lc1():
    basictest(testname, cmdline)
