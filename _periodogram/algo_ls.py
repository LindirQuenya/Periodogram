#!/usr/bin/env python3
# This file contains the core functions for the Lomb-Scargle algorithm.

from _periodogram.func_args import *
from _periodogram.cli_args import *


# computeLombScargle()
# Function to compute the Lomb-Scargle Periodogram for an input light curve
#
# ref: [Scargle, J.D., "Studies in Astronomical Time Series Analysis II.
#     Statistical Aspects of Spectral Analysis of Unevenly Spaced Data."
#     Astrophysical Journal 263:835-853 (1982)];
#     http://adsabs.harvard.edu/full/1982ApJ...263..835S *\/
#
# Periods are sampled according to the time period covered, or based on
# the input values of minperiod and maxperiod.
#
# The coefficients of the transform are selected so the statistical
# distribution of powers for the unevenly spaced power spectrum is the
# same as that of the evenly spaced one.
#
# At each period, a time offset is calculated to diagonalize the
# least-squares fit to sinusoids in the transform.
#
# Power at period p is the magnitude of the transform at p.
#
# Arguments:
#  data = populated dataTbl object
#  args = populated pgramArgs object
#  fargs = funcArgs object with ndata, time, mag
#               nsamp, and period set. power will
#               be populated by this function
#  pool = the pool object to give tasks to
def computeLombScargle(data, args, fargs, pool):
    if not isinstance(pool, mp.pool.Pool):
        raise TypeError('Inappropriate argument type:\npool must be of type Pool!')
    if not isinstance(data, dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args, pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    time = fargs.getTime()
    mag = fargs.getMag()
    period = fargs.getPeriod()
    nsamp = len(period)

    # Compute stats on magnitude
    sdMag = data.getDev(DATA_FIELD_TYPE.DATA_Y)
#    print(sdMag)
    if sdMag == 0:
        raise ValueError("Error in InputFile: Zero deviation in data values!")
    myArgList = []
    startArr, endArr = getSplitNums(nsamp, args.nproc)
    for i in range(len(startArr)):
        myArgList.append([time, mag, sdMag, period, startArr[i], endArr[i]])
    results = pool.starmap(splitLS, myArgList)
    #    results=pool.starmap(doLS,transpose([[time]*nsamp,[mag]*nsamp,[sdMag]*nsamp,period]))
#    print(results)
    fargs.power = merge_arrs(results)


# This function runs a sub-loop of the LS algo.
# It operates on a small slice of the data, and loops over that.
def splitLS(time, mag, sdMag, period, startNum, stopNum):
    ret = []
    for i in range(startNum, stopNum):
        ret.append(doLS(time, mag, sdMag, period[i]))
#    print(ret)
    return ret


# This actually runs the LS algorithm on one period.
def doLS(time, mag, sdMag, p):
    ndata = len(time)
    w = 2 * math.pi / p
    tnum = 0
    tdenom = 0
    for j in range(ndata):
        tnum += math.sin(2.0 * w * time[j])
        tdenom += math.cos(2.0 * w * time[j])
    t = (1 / (2 * w)) * math.atan2(tnum, tdenom)
    lnum = 0
    ldenom = 0
    rnum = 0
    rdenom = 0
    for j in range(ndata):
        s = math.sin(w * (time[j] - t))
        c = math.cos(w * (time[j] - t))
        rnum += mag[j] * s
        lnum += mag[j] * c
        rdenom += s * s
        ldenom += c * c
#    print((1 / (2 * sdMag * sdMag)) * ((lnum * lnum) / ldenom + (rnum * rnum / rdenom)))
    return (1 / (2 * sdMag * sdMag)) * ((lnum * lnum) / ldenom + (rnum * rnum / rdenom))
