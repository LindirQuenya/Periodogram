from testshared import basictest

testname = "AplavIexp"
cmdline = "-a plav -n 500 -s 0.06 -i exp -p 0.1843 -P 1208.8544 -o 10 -N 50 -S 1.0"

def test_AplavIexp_lc1():
    basictest(testname, cmdline)
