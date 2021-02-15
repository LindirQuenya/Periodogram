from testshared import basictest

testname = "AplavIplav"
cmdline = "-a plav -n 500 -s 0.06 -i plav -p 0.1843 -P 1208.8544 -u 0.1 -N 50 -S 1.0"

def test_AplavIplav_lc1():
    basictest(testname, cmdline)
