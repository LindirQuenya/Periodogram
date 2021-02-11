#!/usr/bin/env python3
# This file contains a class for parsing and holding command-line arguments.

from _periodogram.constants import *
from _periodogram.helptext import *
from getopt import getopt
import os
import sys


# Class to hold command line arguments
class pgramArgs:
    def __init__(this):
        this.debugfile = None
        this.debugfp = None

        # this.numProc=0 #number of processors to split job across
        # this.server=None #remote server on which to launch split jobs
        # this.port=0 #remote port

        # this.serverconfig=None #file containing information about servers and
        # number of processors, etc.
        # this.title=None #title for display on plots, UNUSED

        # input filename
        this.intbl = None
        # output filename
        this.outtbl = None

        this.hdu = DEFAULT_HDU  # if the input file is a fits file, specify the hdu
        # UNUSED

        # output directives
        this.inBase = None
        this.outBase = None
        this.outToStdOut = 0  # allow output to standard out
        this.outDir = None  # output files to this directory

        # input directory
        this.datahome = None

        this.localFile = 0  # UNUSED, but keeping it just so old commands won't
        # throw errors
        this.loadPrecomputed = 0  # UNUSED, same as above

        this.periodFile = None  # UNUSED

        # column names, so we know which data means what
        this.xcol = None  # name of x (time) column in data file
        this.ycol = None  # name of y (data) column in data file
        this.yerrCol = None  # name of uncertainty column (y error) in data file
        this.constraintCol = None  # name of column in data file indicating
        # which measurements to exclude (dictated
        # by constraintMin and constraintMax)

        # exclude all rows that have constraint values that are not
        # between (inclusive) these two
        this.constraintMin = DEFAULT_CONSTRAINT
        this.constraintMax = DEFAULT_CONSTRAINT

        this.minperiod = DEFAULT_MAXPERIOD
        this.maxperiod = DEFAULT_MINPERIOD
        this.asFreq = 0  # For error messages, mostly, to check whether
        # period info was given in frequency or period

        # delimiter separating values in inputfile and outputfile
        this.inputDelimiter = DEFAULT_DELIMITER

        # type of period stepping: std,exp,fixedf,fixedp,plav
        this.pstepType = DEFAULT_PSTEP
        this.oversample = DEFAULT_OVERSAMPLE
        this.substep = DEFAULT_SUBSTEP
        this.dfreq = DEFAULT_DFREQ

        # one of "bls", "plav", or "ls"
        this.algo = DEFAULT_ALGO

        # variables for use with BLS algo
        this.nbins = DEFAULT_NBINS
        this.qmin = DEFAULT_QMIN
        this.qmax = DEFAULT_QMAX

        # variables for use with Plavchan algo
        this.nout = DEFAULT_NOUT

        # smoothing box size for plav. algo or output of smoothed curves
        this.smooth = DEFAULT_SMOOTH

        # number of phased light curves to return
        this.nphased = DEFAULT_NPHASED

        # significance threshold for light curve?
        this.sigThresh = DEFAULT_SIG_THRESH

        # statistical quantities to return (and maybe input?)
        this.powN = DEFAULT_POW_NUM
        this.powMean = DEFAULT_POW_MEAN
        this.powSd = DEFAULT_POW_SD

        # Flag to label output files
        this.outLabeled = False

        # Number of workers in pool
        this.nproc = DEFAULT_NPROC

    # Function to parse arguments passed to the "periodogram"
    # command at the command line
    def populate(self):
        nbSet = 0
        qminSet = 0
        qmaxSet = 0
        nProcSet = 0
        noutSet = 0
        bsmSet = 0
        stepSet = 0
        dfSet = 0
        defAlgo = DEFAULT_ALGO
        defStep = DEFAULT_PSTEP
        if len(sys.argv) == 1:
            print(PGRAM_USAGE_TEXT)
            exit()
        elif '--help' in sys.argv:
            print(PGRAM_HELP_TEXT)
            exit()
        # still available: ABCEgGIJjklmOrtUvzZ
        optList, args = getopt(sys.argv[1:], "a:b:c:d:D:e:f:F:h:H:i:K:LM:n:N:o:p:P:q:Q:R:s:S:T:u:V:w:W:x:X:y:Y:")
        for a in optList:
            if a[0] == '-a':
                # check for invalid algo
                if a[1] not in ['ls', 'bls', 'plav']:
                    raise ValueError(
                        "Inappropriate value for option -a:" + a[1] + "\nargument must be one of ls, bls, or plav!")
                # otherwise, set the algo
                self.algo = a[1]
            elif a[0] == '-b':
                # check for invalid number of bins
                if int(a[1]) < 5 or int(a[1]) > 2000:
                    raise ValueError("Inappropriate value for option -b:" + a[1] + "\nargument must be in [5,2000]!")
                # set the number of bins
                self.nbins = int(a[1])
                nbSet = 1
            elif a[0] == '-c':
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -c:" + a[1] + "\nargument must be greater than 0!")
                if float(a[1]) != int(a[1]):
                    raise ValueError("Inappropriate value for option -c:" + a[1] + "\ncannot have fractional workers!")
                self.nproc = int(a[1])
            elif a[0] == '-d':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -d:" + a[1] + "\nargument must be greater than 0!")
                # set the frequency step
                self.dfreq = float(a[1])
                dfSet = 1
            elif a[0] == '-D':
                # Set the debugfile
                self.debugfile = a[1]
                if a[1] == 'stdout':
                    self.debugfp = sys.stdout
                elif a[1] == 'stderr':
                    self.debugfp = sys.stderr
                else:
                    try:
                        self.debugfp = open(a[1], 'w')
                    except:
                        self.debugfp = sys.stderr
                        self.debugfile = 'stderr'
                        print("Unable to open file " + a[1] + ", using stderr.")  # Maybe print this to stderr?
            elif a[0] == '-e':
                # set the input delimiter
                self.inputDelimiter = a[1]
            elif a[0] == '-f':
                # check for conflict, either in value or with maxperiod
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -f:" + a[1] + "\nargument must be greater than 0!")
                if self.maxperiod != DEFAULT_MAXPERIOD:
                    print("Error. Cannot specify both -f and -P! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set maxperiod and asFreq
                self.maxperiod = 1.0 / float(a[1])
                self.asFreq = 1
            elif a[0] == '-F':
                # check for conflict, either in value or with minperiod
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -F:" + a[1] + "\nargument must be greater than 0!")
                if self.minperiod != DEFAULT_MINPERIOD:
                    print("Error. Cannot specify both -F and -p! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set minperiod and asFreq
                self.minperiod = 1.0 / float(a[1])
                self.asFreq = 1
            elif a[0] == '-h':
                # set datahome
                self.datahome = a[1]
            elif a[0] == '-H':
                # Unused
                self.hdu = int(a[1])
                if self.hdu != DEFAULT_HDU and \
                        self.hdu <= 0:
                    raise ValueError("Inappropriate value for option -H:" + a[1] + "\nargument must begreater than 0!")
            elif a[0] == '-i':
                # check for pstepType validity
                if a[1] not in ["exp", "std", 'plav', 'fixedf', 'fixedp']:
                    raise ValueError("Inappropriate value for option -i:" + a[
                        1] + "\nargument must be one of exp, std, plav, or fixedf!")
                # set pstepType
                self.pstepType = a[1]
            elif a[0] == '-I':
                # Implement later?
                pass
            elif a[0] == '-K':
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -K:" + a[1] + "\nargument must be greater than 0!")
                # set stat number of samples
                self.powN = int(a[1])
            elif a[0] == '-l':
                # For server, not implemented.
                pass
            elif a[0] == '-L':
                # must be boolean flag, check.
                if a[1] not in ["0", "1"]:
                    raise ValueError("Inappropriate value for option -L:" + a[1] + "\n argument must be either 1 or 0!")
                # set outLabeled
                self.outLabeled = bool(int(a[1]))
            elif a[0] == '-M':
                # set statmean
                self.powMean = float(a[1])
            elif a[0] == '-n':
                # check that we don't have a negative number of outliers
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -n:" + a[1] + "\nargument must be greater than 0!")
                # set nout and noutSet
                self.nout = int(a[1])
                noutSet = 1
            elif a[0] == '-N':
                # we can't output a negative number of curves!
                if int(a[1]) < 0:
                    raise ValueError("Inappropriate value for option -N:" + a[1] + "\nargument must be greater than 0!")
                # set number of curves to output
                self.nphased = int(a[1])
            elif a[0] == '-o':
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -o:" + a[1] + "\nargument must be greater than 0!")
                # set the oversample factor
                self.oversample = int(a[1])
            elif a[0] == '-O':
                # Implement later?
                pass
            elif a[0] == '-p':
                # check for conflict, in value or with maxfrequency
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -p:" + a[1] + "\nargument must be greater than 0!")
                if self.minperiod != DEFAULT_MINPERIOD:
                    print("Error. Cannot specify both -F and -p! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set minperiod
                self.minperiod = float(a[1])
            elif a[0] == '-P':
                # check for conflict, in value or with minfrequency
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -P:" + a[1] + "\nargument must be greater than 0!")
                if self.maxperiod != DEFAULT_MAXPERIOD:
                    print("Error. Cannot specify both -f and -P! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set maxperiod
                self.maxperiod = float(a[1])
            elif a[0] == '-q':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -q:" + a[1] + "\nargument must be greater than 0!")
                # save the minimum fraction of a period in transit
                self.qmin = float(a[1])
                qminSet = 1
            elif a[0] == '-Q':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -Q:" + a[1] + "\nargument must be greater than 0!")
                # save the maximum fraction of a period in transit
                self.qmax = float(a[1])
                qmaxSet = 1
            elif a[0] == '-r':
                # remote server stuff. not used.
                pass
            elif a[0] == '-R':
                # set the output directory
                self.outDir = a[1]
            elif a[0] == '-s':
                if float(a[1]) <= 0 or float(a[1]) >= 1:
                    raise ValueError("Inappropriate value for option -s:" + a[1] + "\nargument must be in range (0,1)!")
                # set the phase-smoothing box size
                self.smooth = float(a[1])
                bsmSet = 1
            elif a[0] == '-S':
                if float(a[1]) <= 0 or float(a[1]) > 1:
                    raise ValueError("Inappropriate value for option -S:" + a[1] + "\nargument must be in range (0,1]!")
                # set the significance threshold for output of phased curves
                self.sigThresh = float(a[1])
            elif a[0] == '-t':
                # Remote server stuff. Not implemented
                pass
            elif a[0] == '-T':
                # set the title, but replace ' ' with '_'
                # this.title=a[1].replace(' ','_')
                pass
            elif a[0] == '-u':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -u:" + a[1] + "\nargument must be greater than 0!")
                # set the period step factor
                self.substep = float(a[1])
            elif a[0] == '-V':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -V:" + a[1] + "\nargument must be greater than 0!")
                # set the standard deviation for p-values
                self.powSd = float(a[1])
            elif a[0] == '-w':
                # set the constraintMin
                self.constraintMin = float(a[1])
            elif a[0] == '-W':
                # set the constraintMax
                self.constraintMax = float(a[1])
            elif a[0] == '-x':
                # check that the length of the string it was passed is
                # greater than 0
                if len(a[1]) > 0:
                    # set the xcol name
                    self.xcol = a[1]
            elif a[0] == '-X':
                if len(a[1]) > 0:
                    # set the constraintCol name
                    self.constraintCol = a[1]
            elif a[0] == '-y':
                if len(a[1]) > 0:
                    # set the ycol name
                    self.ycol = a[1]
            elif a[0] == '-Y':
                if len(a[1]) > 0:
                    # set the yerrCol name
                    self.yerrCol = a[1]
            else:
                print('Option not recognised. Somehow you bypassed the GetOptError?')
                print(PGRAM_USAGE_TEXT)

        # check mins/maxes are in correct relative positions (min>max)
        if self.qmin >= self.qmax:
            raise ValueError("Error: FractionOfPeriodInTransitMax must be greater than FractionOfPeriodInTransitMin")
        if self.minperiod != DEFAULT_MINPERIOD and \
                self.maxperiod != DEFAULT_MAXPERIOD and \
                self.minperiod >= self.maxperiod:
            if self.asFreq == 1:
                raise ValueError("Error: FrequencyRangeMax must be greater than FrequencyRangeMin!")
            else:
                raise ValueError("Error: PeriodRangeMax must be greater than PeriodRangeMin!")

        # Check that args haven't been set inconsistently
        if nbSet and self.algo != "bls":
            print("Error. Option -b cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if qminSet and self.algo != "bls":
            print("Error. Option -q cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if qmaxSet and self.algo != "bls":
            print("Error. Option -Q cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if dfSet and "fixedf" != self.pstepType != "fixedp":
            print("Error. Option -d cannot be used with this PeriodStepMethod.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if noutSet and self.algo != "plav":
            print("Error. Option -n cannot be used without -a plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if stepSet and self.pstepType != "plav":
            print("Error. Option -u cannot be used without -i plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if self.pstepType == "plav" and self.oversample != DEFAULT_OVERSAMPLE:
            print("Error. Option -o cannot be used with -i plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if dfSet and self.oversample != DEFAULT_OVERSAMPLE:
            print("Error. Option -d cannot be used with -o.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if self.powMean != DEFAULT_POW_MEAN and \
                self.algo == 'ls':
            print("Error. Option -M cannot be used with -a ls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if self.powSd != DEFAULT_POW_SD and \
                self.algo == 'ls':
            print("Error. Option -V cannot be used with -a ls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if (self.constraintCol is None or self.constraintCol == "none") and \
                (self.constraintMin != DEFAULT_CONSTRAINT or
                 self.constraintMax != DEFAULT_CONSTRAINT):
            print("Error. Cannot set constraint extremes without column.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if (self.constraintCol is not None or self.yerrCol is not None or
            self.ycol is not None or self.xcol is not None) and \
                (self.ycol is None or self.xcol is None):
            print("Error. TimeColumn and DataColumn must be set if any columns are set!")
            print(PGRAM_HELP_TEXT)
            exit()

        # read the required input file and  optional output file
        self.intbl = None
        self.outtbl = None
        for s in args:
            # check that we're not misinterpreting a switch, or trying
            # to open "." or ".."
            if s[0] == '-':
                raise ValueError("Misplaced argument " + s + ": switches must precede arguments!")
            elif s[0] == '.' and (len(s) == 1 or (s[1] == '.' and len(s) == 2)):
                if self.intbl == None:
                    raise ValueError("InputFile: " + s + " is not a valid value.")
                else:
                    raise ValueError("OutputFile: " + s + " is not a valid value.")
            else:
                if self.intbl == None:
                    self.intbl = s
                elif self.outtbl == None:
                    self.outtbl = s
                else:
                    raise ValueError("Argument " + s + " not recognized.")
        if self.intbl == None:
            raise ValueError("InputFile: must supply an input file!")
        if self.intbl[0] != '/':
            if self.datahome != None:
                if self.datahome[len(self.datahome) - 1] == '/':
                    self.datahome = self.datahome[:len(self.datahome) - 1]
                self.intbl = self.datahome + '/' + self.intbl
            else:
                self.intbl = os.getcwd() + '/' + self.intbl
        self.getInBase()
        if len(self.inBase) == 0:
            raise ValueError("InputFile: '" + self.intbl + "' does not contain a file name!")
        self.getOutputFile()
        self.getOutBase()

    # Function to get the pathless input filename
    def getInBase(self):
        if self.inBase == None:
            tmp = self.intbl.split("/")
            tmp = tmp[len(tmp) - 1]
            self.inBase = tmp
        return self.inBase

    # Function to get the pathless output filename, creates
    # one if it is not yet set
    def getOutBase(self):
        if self.outBase == None:
            if self.outtbl == None:
                self.outBase = self.getInBase() + ".out"
            else:
                tmp = self.outtbl.split("/")
                tmp = tmp[len(tmp) - 1]
                self.outBase = tmp
        return self.outBase

    # Function to get the complete absolute output filename,
    # sets it if it is unset
    def getOutputFile(self):
        path = os.getcwd()
        if self.outDir != None:
            if self.outDir[0] != '/':
                path = os.getcwd() + '/' + self.outDir
                self.outDir = path
            else:
                path = self.outDir
        if path[len(path) - 1] == '/':
            path = path[:len(path) - 1]
        if self.outtbl == None:
            tmpFname = path + '/' + self.getOutBase()
            if not tmpFname:
                raise RuntimeError("Error. Outtbl filename empty!")
            self.outtbl = tmpFname
        else:
            if self.outtbl[0] != '/':
                tmpFname = path + '/' + self.outtbl
                self.outtbl = tmpFname
        if self.outtbl == None:
            print("path:" + path)
            print("outdir:" + self.outDir)
        return self.outtbl

    # Function to print the arguments for reproducibility.
    # Goes a bit overboard with the specificity, but it'll protect
    # against default value changes.
    def argsPrint(self, main_file_name):
        arr = ['python3', main_file_name]
        if self.algo:
            arr.append('-a')
            arr.append(str(self.algo))
        if self.nbins != DEFAULT_NBINS:
            arr.append('-b')
            arr.append(str(self.nbins))
        if self.nproc:
            arr.append("-c")
            arr.append(str(self.nproc))
        if self.dfreq != DEFAULT_DFREQ and self.pstepType in ["fixedf", "fixedp"]:
            arr.append('-d')
            arr.append(str(self.dfreq))
        if self.inputDelimiter:
            arr.append('-e')
            arr.append("'" + str(self.inputDelimiter) + "'")
        if self.maxperiod != DEFAULT_MAXPERIOD and self.asFreq:
            arr.append('-f')
            arr.append(str(1.0 / self.maxperiod))
        if self.minperiod != DEFAULT_MINPERIOD and self.asFreq:
            arr.append('-F')
            arr.append(str(1.0 / self.minperiod))
        ##        if this.datahome:
        ##            arr.append('-h')
        ##            arr.append(str(this.datahome))
        if self.pstepType:
            arr.append('-i')
            arr.append(str(self.pstepType))
        if self.powN != DEFAULT_POW_NUM:
            arr.append('-K')
            arr.append(str(self.powN))
        if self.outLabeled:
            arr.append('-L')
            arr.append(str(self.outLabeled))
        if self.powMean != DEFAULT_POW_MEAN:
            arr.append('-M')
            arr.append(str(self.powMean))
        if self.nout and self.algo == "plav":
            arr.append('-n')
            arr.append(str(self.nout))
        if self.nphased:
            arr.append('-N')
            arr.append(str(self.nphased))
        if self.oversample != DEFAULT_OVERSAMPLE:
            arr.append('-o')
            arr.append(str(self.oversample))
        if self.minperiod != DEFAULT_MINPERIOD and not self.asFreq:
            arr.append('-p')
            arr.append(str(self.minperiod))
        if self.maxperiod != DEFAULT_MAXPERIOD and not self.asFreq:
            arr.append('-P')
            arr.append(str(self.maxperiod))
        if self.qmin and self.algo == "bls":
            arr.append('-q')
            arr.append(str(self.qmin))
        if self.qmax and self.algo == "bls":
            arr.append('-Q')
            arr.append(str(self.qmax))
        # removed because outDir is included in outtbl
        #        if this.outDir:
        #            arr.append('-R')
        #            arr.append(str(this.outDir))
        if self.smooth:
            arr.append('-s')
            arr.append(str(self.smooth))
        if self.sigThresh:
            arr.append('-S')
            arr.append(str(self.sigThresh))
        #        if this.title:
        #            arr.append('-T')
        #            arr.append(str(this.title))
        if self.substep and self.pstepType == "plav":
            arr.append('-u')
            arr.append(str(self.substep))
        if self.powSd != DEFAULT_POW_SD:
            arr.append('-V')
            arr.append(str(self.powSd))
        if self.constraintMin != DEFAULT_CONSTRAINT:
            arr.append('-w')
            arr.append(str(self.constraintMin))
        if self.constraintMax != DEFAULT_CONSTRAINT:
            arr.append('-W')
            arr.append(str(self.constraintMax))
        if self.xcol:
            arr.append('-x')
            arr.append(str(self.xcol))
        if self.constraintCol:
            arr.append('-X')
            arr.append(str(self.constraintCol))
        if self.ycol:
            arr.append('-y')
            arr.append(str(self.ycol))
        if self.yerrCol:
            arr.append('-Y')
            arr.append(str(self.yerrCol))
        arr.append(self.intbl)
        arr.append(self.outtbl)

        # len(arr) must be at least two, so it's safe to
        # acess the first element without a check
        string = arr[0]
        for i in range(1, len(arr)):
            string += ' ' + arr[i]
        return string
