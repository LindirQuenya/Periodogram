#!/usr/bin/env python3
# John Berberian
# This file contains constants and defaults internal to
# the program. Messing with them isn't recommended unless
# you're sure of what you're doing.

import multiprocessing as mp
import math
from enum import IntEnum

# These are the scaling factors for estimating how long it
# will take the program to run. They were derived on a
# cluster running at 2.8 GHz.
# Adjust as needed (or just to taste).
BLS_T0 = 1.39e-6
LS_T0 = 1.964e-6
PLAV_T0 = 7.878398e-7
# This is a function for estimating how time scales
# with more threads. 1/sqrt(nproc) is generally about
# right, given that single-core turbo frequencies are
# generally much higher than multi-core turbo.
SCALING_FUNC = lambda nproc: 1 / math.sqrt(nproc)

# These are constants. They must be set to these values
# to ensure that the program is able to run properly.
# Change them at your own risk.
UNSET_VALUE = -32768
MAXSTR = 32768
TINY_NUM = 0.0000001
UNSET_MEAN = -1e7
# This class determines which index for every column of dataArray
# signifies each type of data
class DATA_FIELD_TYPE(IntEnum):
    DATA_X = 0
    DATA_X_UNCERTAINTY = 1
    DATA_Y = 2
    DATA_Y_UNCERTAINTY = 3
    DATA_CONSTRAINT = 4
    DATA_N_TYPES = 5

# These are some default values for CLI argument parsing.
DEFAULT_HDU = UNSET_VALUE
DEFAULT_CONSTRAINT = UNSET_VALUE
DEFAULT_MAXPERIOD = UNSET_VALUE
DEFAULT_MINPERIOD = UNSET_VALUE
DEFAULT_PSTEP = "exp"
DEFAULT_OVERSAMPLE = 1
DEFAULT_SUBSTEP = 0.10
DEFAULT_DFREQ = UNSET_VALUE
DEFAULT_NOUT = 500
DEFAULT_SMOOTH = 0.06
DEFAULT_NPHASED = 50
DEFAULT_SIG_THRESH = 1.0
DEFAULT_POW_NUM = UNSET_VALUE
DEFAULT_POW_MEAN = UNSET_VALUE
DEFAULT_POW_SD = UNSET_VALUE
DEFAULT_ALGO = "bls"
DEFAULT_NBINS = -32768
DEFAULT_QMIN = 0.05
DEFAULT_QMAX = 0.10
DEFAULT_DELIMITER = ','
DEFAULT_NPROC = mp.cpu_count()

MIN_NDATA = 2
