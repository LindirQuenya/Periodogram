#!/usr/bin/env python3

# Python multiprocessing periodogram code
# Converted from C++ to Python by John Berberian, Jr.
# Original C++ code by Peter Plavchan
# This version uses multiprocessing, and is the recommended
# multiprocessing version.
# Meant to be run from the terminal with the command:
# user@computer:~$ python3 driver.py [options] <InputFile>
# For more information, run this with the option --help.

import time
from _periodogram.algo_ls import *
from _periodogram.algo_bls import *
from _periodogram.algo_plav import *

# "start" the timer
t_i = time.time()


# computePeriodogram()
# Wrapper routine to compute the periodogram based
# on the input data. FuncArgs is already populated,
# so the messy stuff getting data into the proper format
# for each algo is already done. All that's left is
# calling the algo functions
def computePeriodogram(args, data, fargs, pool):
    if not isinstance(args, pgramArgs):
        raise TypeError("Error: args must be of type pgramArgs!")
    if not isinstance(data, dataTbl):
        raise TypeError("Error: data must be of type dataTbl!")
    if not isinstance(fargs, funcArgs):
        raise TypeError("Error: fargs must be of type funcArgs!")
    if not isinstance(pool, mp.pool.Pool):
        raise TypeError("Error: fargs must be of type funcArgs!")

    # print out a time estimate, calculated earlier
    # by fargs.populate()
    if fargs.timeEst > 0 and PRINT_TIMES:
        print("Estimated time for processing " + str(fargs.nsamp) + " periods: " + format(fargs.timeEst,
                                                                                          '.4f') + ' seconds (' + format(
            fargs.timeEst / 60, '.4f') + ' minutes)')

    # compute the periodogram (results will go into
    # fargs.power)
    if args.algo == "ls":
        computeLombScargle(data, args, fargs, pool)
    elif args.algo == "bls":
        computeBLS(data, args, fargs, pool)
    elif args.algo == "plav":
        computePlavchan(data, args, fargs, pool)
    else:
        print("Error: invalid algo! Somehow you bypassed the argument-population check on this...")
        print(PGRAM_HELP_TEXT)
        exit()
    findPeaks(args, fargs)





