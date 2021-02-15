from testshared import basictest

testname = "AplavIfixedp"
cmdline = "-a plav -n 500 -s 0.06 -i fixedp -p 0.1843 -P 1208.8544 -d 0.2452161 -N 50 -S 1.0"

def test_AplavIfixedp_lc1():
    basictest(testname, cmdline)
