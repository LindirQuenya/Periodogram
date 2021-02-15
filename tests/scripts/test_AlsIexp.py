from testshared import basictest

testname = "AlsIexp"
cmdline = "-a ls -i exp -p 0.1843 -P 1208.8544 -o 10 -N 50 -S 1.0"

def test_AlsIexp_lc1():
    basictest(testname, cmdline)
