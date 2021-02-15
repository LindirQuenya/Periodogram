from testshared import basictest

testname = "AlsIfixedf"
cmdline = "-a ls -i fixedf -p 0.1843 -P 1208.8544 -d 0.0011007 -N 50 -S 1.0"

def test_AlsIfixedf_lc1():
    basictest(testname, cmdline)
