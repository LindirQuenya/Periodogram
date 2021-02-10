#!/usr/bin/env python3
#John Berberian
#This is the preferences file, and it contains flags
#and values that the user can set.

#Flag to enable phase-wrapping in phaseLightCurve()
#Default:True
PHASE_WRAPPING=True

#For some extra print statements, etc.
#Default:False
DEBUG=False

#For BLS, to weight the data points
#according to their uncertainty.
#Default:False
WEIGHT_BY_ERR=False

#To restrict periods to those that we could observe
#in their entirety.
#Default:False
RESTRICT_TO_COMPLETELY_OBSERVABLE=False

#To make minperiod default to the average space
#between data points.
#Default:False
AVG_TIME_STEP=False

#To use the smallest space between two data points
#as minperiod
#Default:False
SMALLEST_TIME_STEP=False

#Flag indicating whether we should offset the periods
#returned by the assorted methods to the range
#specified by minperiod-maxperiod (as opposed to
#omitting periods that lie outside the range)
#Default:False
ADJUST_PERIOD_TO_RANGE=False

#Flag to fit a lognormal distribution to the powers
#returned by BLS
#Default:False
USE_LOGNORMAL_BLS=False

#Flag to ensure that multiple points do not represent
#each peak in computPgramStats()
#Default:False
TRIM_DISTRIBUTION=False

#Flag to trim outliers or not
#Default:False
TRIM_OUTLIERS=False

#Flag to activate print statements that show which
#cycle is currently being completed. Warning:
#generates a lot of output.
#Default:False
TRACK_CYCLES=False

#Flag to print a time estimate before running,
#and the actual elapsed time afterwards. Disable
#if you don't want those messages printing.
#Default:True
PRINT_TIMES = True

#Minimum acceptable internally-calculated minperiod value
#Default:0.1
MIN_MINPERIOD=0.1

#Maximum number of standard deviations away from the mean
#a peice of data can be, and not be labeled an outlier.
#Default:2.5
NORMAL_OUTLIER=2.5

#Maximum fraction of the data that can be labeled as
#an outlier
#Default:0.1
MAX_TRIM_FRACTION=0.1

#Maximum number of iterations for statTrimOutliers()
#Default:10
MAX_ITERATIONS=10