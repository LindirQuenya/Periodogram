#!/usr/bin/env python3
# This file contains a class for parsing and holding command-line arguments.
import argparse

from _periodogram.constants import *
from _periodogram.helptext import *
from getopt import getopt
import os
import sys


# TODO: move these somewhere else
def oneOverFloat(val):
    return 1 / float(val)


def boolInt(val):
    return bool(int(val))


# Class to hold command line arguments
class pgramArgs:
    def __init__(self):
        self.debugfile = None
        self.debugfp = None

        # self.numProc=0 #number of processors to split job across
        # self.server=None #remote server on which to launch split jobs
        # self.port=0 #remote port

        # self.serverconfig=None #file containing information about servers and
        # number of processors, etc.
        # self.title=None #title for display on plots, UNUSED

        # input filename
        self.intbl = None
        # output filename
        self.outtbl = None

        self.hdu = DEFAULT_HDU  # if the input file is a fits file, specify the hdu
        # UNUSED

        # output directives
        self.inBase = None
        self.outBase = None
        self.outToStdOut = 0  # allow output to standard out
        self.outDir = None  # output files to this directory

        # input directory
        self.datahome = None

        self.localFile = 0  # UNUSED, but keeping it just so old commands won't
        # throw errors
        self.loadPrecomputed = 0  # UNUSED, same as above

        self.periodFile = None  # UNUSED

        # column names, so we know which data means what
        self.xcol = None  # name of x (time) column in data file
        self.ycol = None  # name of y (data) column in data file
        self.yerrCol = None  # name of uncertainty column (y error) in data file
        self.constraintCol = None  # name of column in data file indicating
        # which measurements to exclude (dictated
        # by constraintMin and constraintMax)

        # exclude all rows that have constraint values that are not
        # between (inclusive) these two
        self.constraintMin = DEFAULT_CONSTRAINT
        self.constraintMax = DEFAULT_CONSTRAINT

        self.minperiod = DEFAULT_MAXPERIOD
        self.maxperiod = DEFAULT_MINPERIOD
        self.asFreq = 0  # For error messages, mostly, to check whether
        # period info was given in frequency or period

        # delimiter separating values in inputfile and outputfile
        self.inputDelimiter = DEFAULT_DELIMITER

        # type of period stepping: std,exp,fixedf,fixedp,plav
        self.pstepType = DEFAULT_PSTEP
        self.oversample = DEFAULT_OVERSAMPLE
        self.substep = DEFAULT_SUBSTEP
        self.dfreq = DEFAULT_DFREQ

        # one of "bls", "plav", or "ls"
        self.algo = DEFAULT_ALGO

        # variables for use with BLS algo
        self.nbins = DEFAULT_NBINS
        self.qmin = DEFAULT_QMIN
        self.qmax = DEFAULT_QMAX

        # variables for use with Plavchan algo
        self.nout = DEFAULT_NOUT

        # smoothing box size for plav. algo or output of smoothed curves
        self.smooth = DEFAULT_SMOOTH

        # number of phased light curves to return
        self.nphased = DEFAULT_NPHASED

        # significance threshold for light curve?
        self.sigThresh = DEFAULT_SIG_THRESH

        # statistical quantities to return (and maybe input?)
        self.powN = DEFAULT_POW_NUM
        self.powMean = DEFAULT_POW_MEAN
        self.powSd = DEFAULT_POW_SD

        # Flag to label output files
        self.outLabeled = False

        # Number of workers in pool
        self.nproc = DEFAULT_NPROC

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
        optList, args = getopt(sys.argv[1:], "a:b:c:d:D:e:f:F:H:i:K:LM:n:N:o:p:P:q:Q:R:s:S:T:u:V:w:W:x:X:y:Y:")
        for a in optList:
            if a[0] == '-a': # Done
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
            elif a[0] == '-H':  # Capitalized to prevent conflicting with help option.
                # set datahome
                self.datahome = a[1]
            #            elif a[0] == '-H':
            #                # Unused
            #                self.hdu = int(a[1])
            #                if self.hdu != DEFAULT_HDU and \
            #                        self.hdu <= 0:
            #                    raise ValueError("Inappropriate value for option -H:" + a[1] + "\nargument must be greater than 0!")
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
                if int(a[1]) <= 1:
                    raise ValueError("Inappropriate value for option -o:" + a[1] + "\nargument must be greater than 1!")
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
                # self.title=a[1].replace(' ','_')
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
                if self.intbl is None:
                    raise ValueError("InputFile: " + s + " is not a valid value.")
                else:
                    raise ValueError("OutputFile: " + s + " is not a valid value.")
            else:
                if self.intbl is None:
                    self.intbl = s
                elif self.outtbl is None:
                    self.outtbl = s
                else:
                    raise ValueError("Argument " + s + " not recognized.")
        if self.intbl is None:
            raise ValueError("InputFile: must supply an input file!")
        if self.intbl[0] != '/':
            if self.datahome is not None:
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

    def populate_new(self, argv):
        nbSet = False

        parser = argparse.ArgumentParser(description='Compute a periodogram, and compute a period-power table.',
                                         prog='python3 -m periodogram', epilog=EPILOG_TEXT, usage=PGRAM_USAGE_TEXT,
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-a', nargs=1, default=DEFAULT_ALGO, type=str, choices=['ls', 'bls', 'plav'], dest='algo',
                            help='Specifies which Periodogram algorithm to run. One of ls, bls, plav.',
                            metavar='PeriodogramType')
        parser.add_argument('-b', nargs=1, type=int, dest='nbins', metavar='NumberOfBins',
                            help='Specifies the number of bins to use in the bls algorithm (only with -a bls).',
                            max=2000, min=5)
        parser.add_argument('-c', nargs=1, default=DEFAULT_NPROC, type=int, dest='nproc', metavar='NumWorkers',
                            help='Specifies how many workers to use in the multiprocessing Pool. Defaults to the '
                                 '\nnumber of cores on the computer.')
        parser.add_argument('-d', nargs=1, default=DEFAULT_DFREQ, type=float, dest='dfreq', metavar='FixedStepSize',
                            help="Specifies the size of the fixed frequency step or period step, depending on which "
                                 "\nof fixedf or fixedp was chosen for -i. (Only with -i fixedp or -i fixedf)")
        parser.add_argument('-e', nargs=1, default=DEFAULT_DELIMITER, type=str, dest='inputDelimiter',
                            help="Delimiter separating data and (optionally) column labels in the InputFile. If not "
                                 "\nsupplied, ',' will be used. This will also be used to separate the output data in "
                                 "\nthe OutputFile. For information on how to use escape characters from the command "
                                 "\nline (e.g. \\t), see \nhttp://www.gnu.org/software/bash/manual/html_node/ANSI_002dC"
                                 "-Quoting.html", metavar='DataDelimiter')
        maxPGroup = parser.add_mutually_exclusive_group()
        maxPGroup.add_argument('-f', nargs=1, default=DEFAULT_MAXPERIOD, type=oneOverFloat, metavar='FrequencyRangeMin',
                               dest='maxperiod', help='The minimum frequency to consider in period selection.')
        maxPGroup.add_argument('-P', nargs=1, default=DEFAULT_MAXPERIOD, type=float, metavar='PeriodRangeMax',
                               dest='maxperiod', help='The maximum period to consider in period selection.')
        minPGroup = parser.add_mutually_exclusive_group()
        minPGroup.add_argument('-F', nargs=1, default=DEFAULT_MINPERIOD, type=oneOverFloat, metavar='FrequencyRangeMax',
                               dest='minperiod', help='The maximum frequency to consider in period selection.')
        minPGroup.add_argument('-p', nargs=1, default=DEFAULT_MINPERIOD, type=float, metavar='PeriodRangeMin',
                               dest='minperiod', help='The minimum period to consider in period selection.')
        parser.add_argument('-H', nargs=1, default=None, type=str, metavar='DataHome', dest='datahome',
                            help='The directory where the data is located. If not supplied, the current working '
                                 '\ndirectory is used.')
        parser.add_argument('-i', nargs=1, default=DEFAULT_PSTEP, choices=['std', 'exp', 'fixedp', 'fixedf', 'plav'],
                            help='Specifies which period-stepping method to use (one of std, exp, fixedf, fixedp, '
                                 '\nplav).', type=str, dest='pstepType', metavar='PeriodStepMethod')
        parser.add_argument('-K', nargs=1, default=DEFAULT_POW_NUM, dest='powN', metavar='StatNumberOfSamples',
                            help='The number of samples to use for computation of p-values for output peaks. If not '
                                 '\nentered, the number of periods for which power is computed will be used.', type=int)
        parser.add_argument('-L', nargs=1, default=False, dest='outLabeled', type=boolInt, metavar='OutFileLabeled',
                            help='Turns on and off column labels for the output file. Can be either 1 (on) or 0 (off). '
                                 '\nIf not entered, it will default to 0, (off).')
        parser.add_argument('-M', nargs=1, default=DEFAULT_POW_MEAN, dest='powMean', type=float, metavar='StatMean',
                            help='Mean to use for computation of p-values for output peaks.  If not entered, the '
                                 '\nobserved mean will be used.')
        parser.add_argument('-n', nargs=1, default=DEFAULT_NOUT, dest='nout', type=int, metavar='NumberOfOutliers',
                            help='Number of outliers to use in power calculation in the Plavchan algo (only with -a '
                                 '\nplav)')
        parser.add_argument('-N', nargs=1, default=DEFAULT_NPHASED, dest='nphased', metavar='NumberOfPeaksToReturn',
                            help='Limit on the number of top peaks to output in table.', type=int)
        parser.add_argument('-o', nargs=1, default=DEFAULT_OVERSAMPLE, dest='oversample', metavar='OversampleFactor',
                            help='Increase number of periods sampled by this factor (not for use with -i plav or -d)',
                            type=int)
        parser.add_argument('-q', nargs=1, default=DEFAULT_QMIN, metavar='FractionOfPeriodInTransitMin', dest='qmin',
                            help='Minimum fraction of period in transit to consider with BLS algo (-a bls only).',
                            type=float)
        parser.add_argument('-Q', nargs=1, default=DEFAULT_QMAX, metavar='FractionOfPeriodInTransitMax', dest='qmax',
                            help='Maximum fraction of period in transit to consider with BLS algo (-a bls only).',
                            type=float)
        parser.add_argument('-R', nargs=1, default=None, type=str, metavar='OutputDirectory', dest='outDir',
                            help='The directory in which to put output files (periodogram, table of top periods). '
                                 '\nThe default is "."')
        parser.add_argument('-s', default=DEFAULT_SMOOTH, type=float, metavar='PhaseSmoothingBoxSize', dest='smooth',
                            help='Size of box over which to average magnitudes for smoothed curve.', nargs=1)
        parser.add_argument('-S', nargs=1, default=DEFAULT_SIG_THRESH, type=float, metavar='PeakSignificanceThreshold',
                            help='Maximum p-value to accept for output peaks in the power spectrum.', dest='sigThresh')
        parser.add_argument('-u', dest='substep', help='Period increment factor for -i plav.', default=DEFAULT_SUBSTEP,
                            type=float, metavar='PeriodStepFactor', nargs=1)
        parser.add_argument('-V', dest='powSd', default=DEFAULT_POW_SD, type=float, metavar='StatStandardDeviation',
                            help='Standard deviation to use for computation of p-values for output peaks. If not '
                                 '\nentered, the observed standard deviation will be used', nargs=1)
        parser.add_argument('-w', nargs=1, type=float, default=DEFAULT_CONSTRAINT, metavar='ConstraintRangeMin',
                            dest='constraintMin', help='Smallest acceptable value of elements in ConstraintColumn.')
        parser.add_argument('-W', nargs=1, type=float, default=DEFAULT_CONSTRAINT, metavar='ConstraintRangeMax',
                            dest='constraintMax', help='Largest acceptable value of elements in ConstraintColumn.')
        parser.add_argument('-x', nargs=1, type=str, dest='xcol', default=None, metavar='TimeColumn',
                            help='Name of column in input file from which to read time info')
        parser.add_argument('-X', nargs=1, type=str, dest='constraintCol', default=None, metavar='ConstraintColumn',
                            help='Name of column in input file from which to read constraint info')
        parser.add_argument('-y', nargs=1, type=str, dest='ycol', default=None, metavar='DataColumn',
                            help='Name of column in input file from which to read measurement values.')
        parser.add_argument('-Y', nargs=1, type=str, dest='yerrCol', default=None, metavar='DataErrorColumn',
                            help='Name of column in input file from which to read measurement errors.')
        parser.add_argument('InputFile', nargs=1, type=argparse.FileType('r'),
                            help="Text file: first row contains either labels or data, separated by"
                                 "\nthe DataDelimiter. If it contains data, do not supply -x, -X, -y,"
                                 "\nor -Y. Otherwise, supply all that you want used. If a labeled"
                                 "\ncolumn exists, but its label is not associated with a data type"
                                 "\nin the arguments, it will be ignored. Subsequent rows contain data,"
                                 "\nwith different data types separated by the DataDelimiter."
                                 "\nIf column labels are not supplied, the following data types will"
                                 "\nbe assumed:"
                                 "\n          2 columns: TIME,DATA"
                                 "\n          3 columns: TIME,DATA,ERROR"
                                 "\n          4 columns: TIME,DATA,ERROR,CONSTRAINT"
                                 "\n          5+ columns: TIME,DATA,ERROR,CONSTRAINT,UNUSED,UNUSED...")
        parser.add_argument('OutputFile', nargs='?', type=argparse.FileType('w'),
                            help="Text file. All content will be overwritten. Columns will"
                                 "\nbe labeled according to OutFileLabeled, with a column format of"
                                 "\nPERIOD,POWER. If no file is specified, one will be constructed in the"
                                 "\noutput directory with the name '<path-free name of input file>.out")
        args = parser.parse_args(argv)
        self.algo = args.algo;
        if args.nbins is None:
            args.nbins = DEFAULT_NBINS
        else:
            nbSet = True
        if not 2000 >= args.nbins >= 5:
            raise ValueError("Inappropriate value for option -b:" + args.nbins + "\nargument must be in [5,2000]!")
        self.nbins = args.nbins;

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


    # Function to get the pathless input filename
    def getInBase(self):
        if self.inBase is None:
            tmp = self.intbl.split("/")
            tmp = tmp[len(tmp) - 1]
            self.inBase = tmp
        return self.inBase

    # Function to get the pathless output filename, creates
    # one if it is not yet set
    def getOutBase(self):
        if self.outBase is None:
            if self.outtbl is None:
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
        if self.outDir is not None:
            if self.outDir[0] != '/':
                path = os.getcwd() + '/' + self.outDir
                self.outDir = path
            else:
                path = self.outDir
        if path[len(path) - 1] == '/':
            path = path[:len(path) - 1]
        if self.outtbl is None:
            tmpFname = path + '/' + self.getOutBase()
            if not tmpFname:
                raise RuntimeError("Error. Outtbl filename empty!")
            self.outtbl = tmpFname
        else:
            if self.outtbl[0] != '/':
                tmpFname = path + '/' + self.outtbl
                self.outtbl = tmpFname
        if self.outtbl is None:
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
        #        if self.datahome:
        #            arr.append('-h')
        #            arr.append(str(self.datahome))
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
        #        if self.outDir:
        #            arr.append('-R')
        #            arr.append(str(self.outDir))
        if self.smooth:
            arr.append('-s')
            arr.append(str(self.smooth))
        if self.sigThresh:
            arr.append('-S')
            arr.append(str(self.sigThresh))
        #        if self.title:
        #            arr.append('-T')
        #            arr.append(str(self.title))
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
