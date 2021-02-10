# Python multiprocessing periodogram code
# Converted from C++ to Python by John Berberian, Jr.
# Original C++ code by Peter Plavchan
# This version uses multiprocessing, but is very inefficient.
# mp_periodogram.py is recommended over this program.
# Meant to be run from the terminal with the command:
# user@computer:~$ python3 mp_alt.py [options] <InputFile>
# For more information, run this with the option --help.

import multiprocessing as mp
import os
import sys
import time
from enum import IntEnum
from getopt import getopt
from constants import *
from preferences import *

# "start" the timer
t_i = time.time()


#######################
##
##  Algorithm Functions
##
#######################

# assumes rectangular matrix
def transpose(l):
    n = []
    for i in range(len(l)):
        for j in range(len(l[0])):
            if j == len(n):
                n.append([l[i][j]])
            else:
                n[j].append(l[i][j])
    return n


def makeArrOfCopies(inList, num):
    bigArr = [[]] * num
    for k in inList:
        for i in range(len(bigArr)):
            bigArr[i].append(k)
    return bigArr


def computeBLS(data, args, fargs, pool):
    # Check for type errors
    if not isinstance(data, dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args, pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    if not isinstance(pool, mp.pool.Pool):
        raise TypeError('Inappropriate argument type:\npool must be of type Pool!')
    [ndata, time] = funcArgsGetTime(fargs)
    [ndata, mag] = funcArgsGetMag(fargs)
    [nsamp, period] = funcArgsGetPeriods(fargs)
    [nsamp, power] = funcArgsGetPower(fargs)

    # Initializes variables with nonsense values
    blsR = [0.0] * nsamp
    blsS = [0.0] * nsamp
    lowBin0 = [0] * nsamp
    lowBin1 = [0] * nsamp

    wt = []
    # In case we want weight as a function of uncertainty
    if WEIGHT_BY_ERR:
        err = dtGetFilteredArray(data, DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY)
    totalWt = 0
    for j in range(ndata):
        if WEIGHT_BY_ERR:
            wt.append(err[j])
        else:
            # If we don't, weight everything the same
            wt.append(1)
        totalWt += wt[j]

    # If unset, make a value for nbins, based off of ndata
    if args.nbins == DEFAULT_NBINS:
        if ndata <= 500:
            args.nbins = 50
        elif ndata <= 20000:
            args.nbins = int(ndata / 10)
        else:
            args.nbins = 2000
    nbins = args.nbins

    # Write nbins value to debugfp
    if args.debugfp != None:
        args.debugfp.write("IN COMPUTE BLS: nbins = " + str(nbins) + "\n")

    # Run checks on FractionOfPeriodInTransitMin/Max before
    # running the algorithm, in case it somehow escaped our
    # pgramArgs.populate() checks
    qmin = args.qmin
    qmax = args.qmax
    if qmin <= 0 or qmax <= 0 or qmax < qmin:
        raise ValueError("Error: invalid values for qmin/qmax!")

    minBins = int(qmin * nbins)
    if minBins < 1:
        minBins = 1

    minWt = totalWt * qmin
    # I don't love that this is fixed at "5" -- if we convert to
    # weighting with errors, it'll have to change
    if minWt < 5:
        minWt = 5  # min weight over "low" set of bins

    # maximum number of abins over which a "low" phase can extend:
    # (this is also the amount by which we want to pad the bin array)
    binExt = int(qmax * nbins) + 1
    binMax = int(nbins + binExt)

    # Initialize arrays with null values, because we refer to
    # specific elements of these arrays by index later on:
    # just appending as needed won't work here.
    binMag = [0.0] * binMax
    binWt = [0.0] * binMax

    # Compute periodogram
    results = pool.starmap(doBLS, transpose([[time] * nsamp, [mag] * nsamp, [wt] * nsamp, \
                                             makeArrOfCopies(binWt, nsamp), \
                                             makeArrOfCopies(binMag, nsamp), \
                                             period, [nbins] * nsamp, [nsamp] * nsamp, \
                                             [binExt] * nsamp, [minBins] * nsamp, \
                                             [minWt] * nsamp, [totalWt] * nsamp]))
    results = transpose(results)
    fargs.power = results[0]
    fargs.blsR = results[1]
    fargs.blsS = results[2]
    fargs.lowBin0 = results[3]
    fargs.lowBin1 = results[4]


# time and mag are arrs. binWt and binMag are working vars,
# so they can't be shared. copies will have to be made for
# each run of the function.
def doBLS(time, mag, wt, binWt, binMag, period, nbins, nsamp, binExt, minBins, minWt, \
          totalWt):
    maxPwr = 0.0
    for b in range(nbins):
        binMag[b] = 0.0
        binWt[b] = 0.0
    ndata = len(time)
    binMax = len(binWt)
    for j in range(ndata):
        phase = (time[j] / period) % 1
        b = int(math.floor(nbins * phase))
        binWt[b] += wt[j]
        binMag[b] += wt[j] * mag[j]
    for b in range(nbins, binMax):
        binWt[b] = binWt[b - nbins]
        binMag[b] = binMag[b - nbins]
    maxPwr = 0.0
    for b in range(nbins):
        binCt = 0
        sumWt = 0.0
        sumMag = 0.0
        for k in range(b, b + binExt + 1):
            binCt += 1
            sumWt += binWt[k]
            sumMag += binMag[k]
            if binCt >= minBins and sumWt >= minWt and \
                    sumWt < totalWt:
                pwr = (sumMag ** 2) / (sumWt * (totalWt - sumWt))
                if pwr >= maxPwr:
                    maxPwr = pwr
                    lowStart = b
                    lowEnd = k
                    lowWt = sumWt
                    lowMag = sumMag
    maxPwr = math.sqrt(maxPwr)
    if maxPwr > 0:
        return (maxPwr, lowWt / totalWt, lowMag, lowStart, lowEnd)
    return (0, 0, 0, 0, 0)


def computeLombScargle(data, args, fargs, pool):
    if not isinstance(pool, mp.pool.Pool):
        raise TypeError('Inappropriate argument type:\npool must be of type Pool!')
    if not isinstance(data, dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args, pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    [ndata, time] = funcArgsGetTime(fargs)
    [ndata, mag] = funcArgsGetMag(fargs)
    [nsamp, period] = funcArgsGetPeriods(fargs)
    [nsamp, power] = funcArgsGetPower(fargs)

    # Compute stats on magnitude
    sdMag = dtGetDev(data, DATA_FIELD_TYPE.DATA_Y)
    if sdMag == 0:
        raise ValueError("Error in InputFile: Zero deviation in data values!")
    results = pool.starmap(doLS, transpose([[time] * nsamp, [mag] * nsamp, [sdMag] * nsamp, period]))
    fargs.power = results


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
    return (1 / (2 * sdMag * sdMag)) * ((lnum * lnum) / ldenom + (rnum * rnum / rdenom))


def phaseLightCurve(time, mag, mySmooth, boxSize, p):
    ndata = len(time)
    mySort = []
    myPhase = [0] * ndata
    myMag = [0] * ndata
    myChi = [0] * ndata
    for j in range(ndata):
        mySort.append([mod(time[j], p) / p, j])
    mySort.sort()
    for j in range(ndata):
        myIdx = mySort[j][1]
        myPhase[j] = mySort[j][0]
        myMag[j] = mag[myIdx]
        if mySmooth:
            mySmooth[j] = -1
        if j > 0 and myPhase[j] < myPhase[j - 1]:
            print("Sort error?")
    if mySmooth:
        # Now smooth the phases:
        #
        # To reduce computation time, remember values from the
        # last round. Specifically:

        # Given that phase[j]>=phase[j-1] (data is sorted)

        # If
        # (phase[j-1] - phase[prevLo-1]) > smooth/2 then
        # (phase[j] - phase[prevLo-1]) > smooth/2 and, generally,
        # (phase[j] - phase[k]) > smooth/2 for any k<prevLo

        # -> there's no need to consider any k<prevLo as the lower
        # edge of the box.

        # If
        # (phase[prevHi] - phase[j-1]) <= smooth/2 then
        # (phase[prevHi] - phase[j]) <= smooth/2 and, generally,
        # (phase[k] - phase[j]) <= smooth/2 for any k <= prevHi

        # -> there's no need to consider any k<prevHi as the upper
        # edge of the box

        ###

        # "Wrapping":
        # Phase-smoothing can still take place for values at the beginning
        # or end of a period by "wrapping" around to the other end of the
        # array. For example, we could smmoth values for phase 0 with
        # those at phase 1-s, or values at phase 1 with those at
        # phase 0+s. Since the assumption in phase-folding is that the
        # signal is periodic, looking at the other end of the array is like
        # wrapping around to the "next" or "previous" period
        bLo = 0;
        bHi = 0;
        prevLo = 0;
        prevHi = 0;
        count = 0
        boxSum = 0.0;
        halfbox = boxSize / 2.0

        for j in range(ndata):
            if PHASE_WRAPPING:
                # Determine value for bHi. If j=0, prevHi=0,
                # thereafter updated to the largest index such
                # that phase[bHi] - phase[j-1]<halfBox
                for bHi in range(prevHi, 2 * ndata - 1):
                    if bHi < (ndata - 1):
                        if (myPhase[bHi + 1] - myPhase[j]) > halfbox:
                            break
                    else:
                        # adjust index and phase range if we're off the
                        # end of the array
                        if (myPhase[bHi - ndata + 1] - myPhase[j] + 1) > halfbox:
                            break
                if (j == 0):
                    # Get initial value for bLo for j=0: step back until
                    # ndata + (bLo -1) is out of the range
                    for bLo in range(0, -ndata, -1):  # +1 removed because of >=
                        if (myPhase[j] - myPhase[ndata + bLo - 1] + 1) > halfbox:
                            break
                else:
                    # once prevLo hase ben properly initialized (i.e. j>0),
                    # start from prevLo and shift box as needed
                    for bLo in range(prevLo, j):
                        if bLo < 0:
                            if (myPhase[j] - myPhase[ndata + bLo] + 1) < halfbox:
                                break
                        else:
                            if (myPhase[j] - myPhase[bLo]) < halfbox:
                                break
            else:
                # Find edges of box (no phase-wrapping)
                for bLo in range(prevLo, j):
                    if (myPhase[j] - myPhase[bLo]) < halfbox:
                        break
                for bHi in range(prevHi, ndata - 1):
                    if (myPhase[bHi + 1] - myPhase[j]) > halfbox:
                        break
                if DEBUG:
                    # check that the values for bLo and bHi satisfy
                    # the conditions we were trying to meet
                    if ((myPhase[bHi] - myPhase[j]) > halfbox) or \
                            (((bHi + 1) < ndata) and ((myPhase[bHi + 1] - myPhase[j]) <= halfbox)) \
                            or ((myPhase[j] - myPhase[bLo]) > halfbox) or \
                            ((bLo > 0) and ((myPhase[j] - myPhase[bLo - 1]) <= halfbox)):
                        print("BOX ERROR!")

            # Initialize the sum of magnitudes in our box. We will
            # start fro scratch if j = 0 or if our new low edge is
            # above our previous high edge
            if (j == 0) or (bLo >= prevHi):
                boxSum = 0
                count = 0
                for k in range(bLo, bHi + 1):  # shifted +1 due to <=
                    if k < 0:
                        myIdx = ndata + k
                    elif k >= ndata:
                        myIdx = k - ndata
                    else:
                        myIdx = k
                    boxSum += myMag[myIdx]
                    count += 1
            else:
                # if there is overlap between this box and the
                # previous one, subtract off the left edge and
                # add the right
                for k in range(prevLo, bLo):
                    if k < 0:
                        myIdx = ndata + k
                    else:
                        myIdx = k
                    boxSum -= myMag[myIdx]
                    count -= 1
                for k in range(prevHi + 1, bHi + 1):  # shifted+= due to <=
                    if k >= ndata:
                        myIdx = k - ndata
                    else:
                        myIdx = k
                    boxSum += myMag[myIdx]
                    count += 1

            # save the values of bLo and bHi for the next value of j
            prevLo = bLo
            prevHi = bHi
            if count > 0:
                mySmooth[j] = boxSum / count
                if myChi:
                    myChi[j] = (myMag[j] - mySmooth[j]) ** 2
    return myChi


def computePlavchan(data, args, fargs, pool):
    # Check for type errors
    if not isinstance(data, dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args, pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')

    [ndata, time] = funcArgsGetTime(fargs)
    [ndata, mag] = funcArgsGetMag(fargs)
    [ndata, smooth] = funcArgsGetSmoothedMag(fargs)
    [nsamp, period] = funcArgsGetPeriods(fargs)
    [nsamp, power] = funcArgsGetPower(fargs)
    noutliers = args.nout

    # array to hold the deviation from the smoothed curve for each
    # data point
    [ndata, tmpChi] = funcArgsGetChi(fargs)

    # make sure we don't have more outliers than we have data points
    if noutliers > ndata:
        noutliers = ndata

    # determine reference deviations (recycle "tmpChi" array)
    meanMag = dtGetMean(data, DATA_FIELD_TYPE.DATA_Y)
    for j in range(ndata):
        tmpChi[j] = (mag[j] - meanMag) ** 2
    tmpChi.sort()

    # sum the values most _poorly_ fit by the model mag=meanMag
    maxStd = 0.0;
    maxChi = 0.0  # maxChi is the analogous var for each pd
    for j in range(ndata - 1, ndata - noutliers - 1, -1):
        maxStd += tmpChi[j]
    maxStd /= noutliers

    boxSize = fargs.boxSize

    # Compute periodogram
    fargs.power = pool.starmap(doPlav, transpose([[time] * nsamp, [mag] * nsamp, \
                                                  makeArrOfCopies(smooth, nsamp), \
                                                  [boxSize] * nsamp, [maxStd] * nsamp, \
                                                  [noutliers] * nsamp, period]))


def doPlav(time, mag, mySmooth, boxSize, maxStd, noutliers, p):
    ndata = len(time)
    errval = 0
    tmpChi = phaseLightCurve(time, mag, mySmooth, boxSize, p)
    tmpChi.sort()
    count = 0
    maxChi = 0
    for j in range(ndata - 1, -1, -1):
        if tmpChi[j] != errval:
            maxChi += tmpChi[j]
            count += 1
            if count >= noutliers:
                break
    maxChi /= count
    if maxChi > 0:
        return maxStd / maxChi
    return errval


#######################
##
##  Fargs Getter functions
##
#######################
# Returns [fargs.ndata,fargs.time]
def funcArgsGetTime(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.time]


# Returns [fargs.ndata,fargs.mag]
def funcArgsGetMag(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.mag]


# Returns [fargs.ndata,fargs.period]
def funcArgsGetPeriods(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.nsamp, fargs.period]


# Returns [fargs.ndata,fargs.power]
def funcArgsGetPower(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.nsamp, fargs.power]


# Returns [fargs.ndata,fargs.chi]
def funcArgsGetChi(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.chi]


# Returns [fargs.ndata,fargs.phase]
def funcArgsGetPhase(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.phase]


# Returns [fargs.ndata,fargs.phasedMag]
def funcArgsGetPhasedMag(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.phasedMag]


# Returns [fargs.ndata,fargs.smoothedMag]
def funcArgsGetSmoothedMag(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.smoothedMag]


# Returns [fargs.ndata,fargs.sortable]
def funcArgsGetSortable(fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    return [fargs.ndata, fargs.sortable]


# funcArgsGetPeakWidth()
# Function to compute and return the start and end+1 indices of each "peak."
# width[i][0] is the starting index of the current peak and width[i][1] is
# its end index+1 (thus, width[i][0]<=i<width[i][1]).
# The start and end indices of a peak are as follows:
# For each i, width[i][0]=i
# - if power[i] is not above the mean, width[i][1] = i+1.
# - if it is above the mean (or above it by "siglevel", currently 0) find
#   j such that power[i]...power[j-1] are all above the mean and power[j]
#   is below it. Then set width[i][1]=j.
# Continue the above loop until i=j
def funcArgsGetPeakWidth(fargs, useLog):
    sigLevel = 0

    # retrieve power array
    [nsamp, power] = funcArgsGetPower(fargs)

    if useLog:
        tmp = []
        for i in range(nsamp):
            if not power[i]:
                tmp.append(math.log(1.0e-7))
            else:
                tmp.append(math.log(power[i]))
        power = tmp

    # compute mean and dev for power
    meanPow = statMean(power, None)
    sdPow = statStdDev(power, None)

    # make width be a 2d array filled with nonsense
    # so that we can refer to indecies later on
    width = [[0.0, 0.0]] * nsamp

    for i in range(nsamp):
        width[i] = [i, i + 1]

        # If stats failed for this set of powers, make width 1 for all
        if sdPow > 0:
            j = i
            while j < nsamp and (power[j] - meanPow) / sdPow > sigLevel:
                j += 1
                width[i] = [width[i][0], j]

            # catch up: everything between i and j is the same peak
            for k in range(i + 1, j):
                width[k] = [i, j]
            if j > i:
                i = j - 1
    fargs.width = width
    return width


#######################
##
##  Data table functions
##
#######################

# Returns the mean for data type f in data table t.
# Calculates and sets the mean if unset.
def dtGetMean(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')

    # If the mean hasn't been calculated yet, do so and set it.
    if isMeanSet(t, f) < 1:
        myArray = dtGetArray(t, f)
        m = statMean(myArray, t.omitRow)
        setMean(t, f, m)
    return t.mean[f]


# Returns the data array for data type f given by data table t
def dtGetArray(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    return t.dataArray[f]


# Returns the data array for data type f given by data table t,
# filtered by t.omitRow: omitRow is 0 or None for good rows, and
# any other value for ones that should be ignored.
def dtGetFilteredArray(t, f):
    # No need to run checks, dtGetArray should do that.
    refArray = dtGetArray(t, f)
    myArray = []
    for i in range(t.ndata):
        if t.omitRow == None or not t.omitRow[i]:
            myArray.append(refArray[i])
    return myArray


# Returns the mininum value of the filtered array for data type f
# given by data table t
def dtGetMin(t, f):
    if t.ndata <= 0:
        raise ValueError("Inappropriate value:\nt must be populated before dtGetMin can be called!")
    arr = dtGetFilteredArray(t, f)
    return min(arr)


# Returns the maximum value of the filtered array for data type f
# given by data table t
def dtGetMax(t, f):
    if t.ndata <= 0:
        raise ValueError("Inappropriate value:\nt must be populated before dtGetMin can be called!")
    arr = dtGetFilteredArray(t, f)
    return max(arr)


# Returns the standard deviation for data type f in data table t.
# Calculates and sets it if it is unset
def dtGetDev(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of type dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    if isDevSet(t, f) < 1:
        arr = dtGetArray(t, f)
        m = statStdDev(arr, t.omitRow)
        setDev(t, f, m)
    return t.sd[f]


# Returns the minimum difference between two adjacent elements
# in the array t.dataArray[f].
def dtGetMinDiff(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of type dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    ndata = t.ndata
    myArray = dtGetArray(t, f)
    minDiff = myArray[1] - myArray[0]
    if minDiff == 0:
        raise ValueError('InputFile: two measurements cannot exist at the same time!')
    diff = 0
    for i in range(1, ndata):
        diff = myArray[i] - myArray[i - 1]
        if diff < 0:
            raise ValueError("Array must be sorted before calling dtGetMinDiff!")
        if diff > 0 and diff < minDiff:
            minDiff = diff
    return minDiff


# Returns the median difference between two adjacent elements
# in the array t.dataArray[f].
def dtGetMedianDiff(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of type dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    ndata = t.ndata
    myArray = dtGetArray(t, f)
    diff = []
    for i in range(1, ndata):
        diff.append(myArray[i] - myArray[i - 1])
        if diff[i - 1] < 0:
            raise ValueError("Array must be sorted before calling dtGetMinDiff!")
    return statMedian(diff, None)


#######################
##
##  Data table value checkers
##
#######################

# Returns a value signifying whether the mean in dataTbl t
# for data type signified by f is set. -1:error,0:unset,1:set
def isMeanSet(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        return -1
    if t.mean[f] == UNSET_MEAN:
        return 0
    return 1


# Returns a value representing whether the standard deviation has
# been set for data type f in data table t
def isDevSet(t, f):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of type dataTbl!')
    if t == None or f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        return -1
    if t.sd[f] == UNSET_MEAN or t.sd[f] == 0:
        return 0
    return 1


#######################
##
##  Data table value setters
##
#######################

# Sets the mean for data type f in data table t to m
def setMean(t, f, m):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    t.mean[f] = m


# Sets the standard deviation for data type f in data table t to d
def setDev(t, f, d):
    if not isinstance(t, dataTbl):
        raise TypeError('Inappropriate argument type:\nt must be of dataTbl!')
    if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nf is greater than maximum value!')
    t.sd[f] = d


#######################
##
##  Statistical functions
##
#######################

# Returns the mean of arr, filtering unusable data points specified by omitRow:
# if omitRow[n] is nonzero, row n should be ignored.
def statMedian(arr, omitRow):
    if not isinstance(arr, list):
        raise TypeError("Inappropriate argument type:\narr must be a list!")
    if len(arr) <= 0:
        raise IndexError("Inappropriate index:\ncannot find the median of an empty list!")
    tmp = []
    count = 0
    for i in range(len(arr)):
        if omitRow == None or omitRow[i] == 0:
            tmp.append(arr[i])
            count += 1
    if count <= 0:
        raise ValueError("No usable points in array!")
    tmp.sort()
    if count % 2 == 1:
        md = tmp[int(count - 1) / 2]
    else:
        md = (tmp[int(count / 2) - 1] + tmp[int(count / 2)]) / 2.0
    return md


# Returns the mean of arr. Filters bad data specified by omitRow
def statMean(arr, omitRow):
    if not isinstance(arr, list):
        raise TypeError("Inappropriate argument type:\narr must be a list!")
    if len(arr) <= 0:
        raise IndexError("Inappropriate index:\ncannot find the mean of an empty list!")
    if omitRow:
        arr1 = [arr[i] for i in range(len(arr)) if not omitRow[i]]
    else:
        arr1 = arr
    return sum(arr1) / len(arr1)


# Returns the standard deviation of an array
# Filters bad data specified by omitRow
def statStdDev(inArr, omitRow):
    if omitRow:
        arr = [inArr[i] for i in range(len(inArr)) if not omitRow[i]]
    else:
        arr = inArr
    s = sum(arr)
    sumsq = 0
    for i in arr:
        sumsq += i ** 2
    if len(arr) <= 0:
        raise ValueError("Cannot find standard deviation of empty array!")
    # "divide by n" method
    ##return math.sqrt(len(arr)*sumsq-s**2)/len(arr)
    if len(arr) > 1:
        diff = sumsq - (s * s / len(arr))
        if diff < 0:
            if diff > TINY_NUM:
                raise RuntimeError("Calculation error in std dev!")
            else:
                diff = 0
        return math.sqrt(diff / (len(arr) - 1))
    else:
        # Non-fatal error, note it.
        return -1


# statTrimOuliers()
# remove outliers from the input array, where an outlier is defined
# as a value more than "threshold" standard deviations from the mean.
# if more than maxFract samples are removed, terminate
def statTrimOutliers(length, x, threshold, maxFract):
    nout = 100
    minCount = int(math.ceil((1.0 - maxFract) / length))

    x1 = [0.0] * length
    x2 = [0.0] * length
    its = 0
    for i in range(length):
        x1[i] = x[i]
    while (nout > 1) and (its < MAX_ITERATIONS):
        its += 1
        nout = 0
        nin = 0
        mean = statMean(x1, None)
        sd = statStdDev(x1, None)
        if sd <= 0:
            # something went wrong
            raise RuntimeError("Calculation error in std dev!")
        for i in range(length):
            if abs((x1[i] - mean) / sd) <= threshold:
                x2[nin] = x1[i]
                nin += 1
        # accept this iteration if we have at least minCount samples
        if nin > minCount:
            out = length - nin
            length = nin
            for i in range(length):
                x1[i] = x2[i]
        else:
            break
    x1 = x1[:length]
    return length, x1, mean, sd


#######################
##
##  Misc. Functions
##
#######################

# Function to read the data from a file, and construct dataArray from it,
# using user-defined column names, or the default ones (see help message)
# Ignores any rows with non-numerical data.
def readDoubleColumns(fname, arrayNames, delim):
    # open the file
    f = open(fname, 'r')
    data = f.readlines()
    f.close()

    # strip the newline character from each line and splits
    # it according to the delimiter
    for i in range(len(data)):
        data[i] = data[i].strip().split(delim)

    # checks that file is not empty
    if len(data) == 0:
        raise ValueError("InputFile cannot be empty!")

    # sets the default value for numCol, the number of
    # columns "deep" that we will read into a row
    numCol = len(data[0])
    if arrayNames[DATA_FIELD_TYPE.DATA_X] == None:
        # No labels, first row is data.
        startRow = 0
        if len(data[0]) < 2:
            raise ValueError("InputFile: too few columns for meaningful data!")
        # Make an index map. this is the default one detailed in
        # the help message.
        indMap = [DATA_FIELD_TYPE.DATA_X, DATA_FIELD_TYPE.DATA_Y, DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY,
                  DATA_FIELD_TYPE.DATA_CONSTRAINT]
        # Ensure that we do not read deeper into either array
        # than its length
        numCol = min(len(data[0]), len(indMap))
    else:
        # Labels in first row, second row+ is data.
        startRow = 1
        # In this case, we make an index map the length of numCol, and
        # construct non-None elements where nessesary
        indMap = [None] * numCol
        # for every index in the label row (first row)
        for i in range(len(data[0])):
            try:
                # get the label there
                label = data[0][i]
                # find its index in arrayNames (this is where an exception
                # would be thrown, if the label does not exist in arrayNames)
                ind = DATA_FIELD_TYPE(arrayNames.index(label))
                # change that index map element
                indMap[i] = ind
            except:
                # Evidently that label is not one that we care about
                pass

    # Initialize columns array with an array for each data type
    columns = [[]] * DATA_FIELD_TYPE.DATA_N_TYPES

    # Iterate through the rows
    for rownum in range(startRow, len(data)):
        row = data[rownum]
        # Iterate through the elements of that row
        for i in range(numCol):
            # If that index has been mapped to a data type
            if indMap[i] != None:
                try:
                    # Try casting to float. This is where an error would be
                    # thrown.
                    fl = float(row[i])
                    # Add it on to columns at the proper index
                    columns[indMap[i]] = columns[indMap[i]] + [fl]
                except:
                    # It must be non-numerical. Let's pretend that row
                    # never existed.
                    pass
    # It's possible that some types of data just weren't in the file,
    # so change all of the columns that are still empty to be filled
    # with None. Side note: python list comprehesion is awsome!
    columns = [[None] * len(columns[0]) if n == [] else n for n in columns]

    # Finally, return our value
    return columns


def sortColumns(arr2d, sortIdx):
    # Check that the sortIdx is valid
    if sortIdx >= DATA_FIELD_TYPE.DATA_N_TYPES:
        raise ValueError('Inappropriate value:\nsortIdx is greater than maximum value!')

    # Make a temporary array for sorting
    sortArr = []
    for i in range(len(arr2d[sortIdx])):
        # first element of each element in sortArr is the value, the
        # second is its initial index, to keep track of where each element
        # went.
        sortArr.append([arr2d[sortIdx][i], i])

    # Sort sortArr. This will sort by the first element (the value),
    # and only use the second element (the index) as a tie-breaker,
    # which is fine for our purposes
    sortArr.sort()

    # Now make all of the elements of the columns in arr2d rearrange
    # themselves to match the order of sortArr
    for i in range(len(arr2d)):
        # make a temporary column variable
        col = []
        for k in sortArr:
            # fill it with the elements in the right order
            col.append(arr2d[i][k[1]])
        # and set it
        arr2d[i] = col

    # Return the now-sorted array
    return arr2d


# Basic algo to save period/power columns, unlabeled
def saveOutput(outstream, period, power, delim):
    if len(period) != len(power):
        raise RuntimeError("Length of period data not equal to length of power data!")
    for i in range(len(period)):
        # write each value to the outfile
        outstream.write(format(period[i], '.11f') + delim + format(power[i], '.11f') + '\n')


# Basic algo to save period/power columns, labeled
def saveLabeledOutput(outstream, period, power, periodLabel, powerLabel, delim):
    # run saveOutput, but append the labels to the beginning of each column
    outstream.write(periodLabel + delim + powerLabel)
    saveOutput(outstream, period, power, delim)


# estimateProcessingTime()
# Function to estimate the time required to execute the command with
# the input values of nsamp, ndata, nbins, and qmax.
# The relative time spent in the different loops of BLS is empirically
# estimated. For really reliable estimates, might do linear regression.
def estimateProcessingTime(nsamp, ndata, args, qmax, algo):
    if algo == "bls":
        if args.nbins == DEFAULT_NBINS:
            if ndata <= 500:
                args.nbins = 50
            elif ndata <= 20000:
                args.nbins = int(ndata / 10)
            else:
                args.nbins = 2000
        scaledNsamp = BLS_T0 * nsamp
        timeEst = (scaledNsamp * (ndata + (2.096 * args.nbins) + ((qmax * args.nbins) * (-2.6 + 0.38 * args.nbins))))
    elif algo == "ls":
        timeEst = ((LS_T0 * nsamp) * ndata)
    elif algo == "plav":
        timeEst = (((PLAV_T0 * nsamp) * ndata) * math.log2(ndata))
    else:
        raise ValueError("PeriodogramType: invalid value " + algo)
    if timeEst == 0:
        # this must have been a non-fata error: note it
        timeEst = -1
    else:
        timeEst *= SCALING_FUNC(args.nproc)
    return timeEst


def findPeaks(args, fargs):
    if not isinstance(fargs, funcArgs):
        raise TypeError("Error: fargs must be of type funcArgs!")
    if not isinstance(args, pgramArgs):
        raise TypeError("Error: args must be of type pgramArgs!")

    # Make lists to hold significant periods and their powers
    sigPeriod = []
    sigPower = []

    isBls = 0;
    isLs = 0;
    isPlav = 0
    if args.algo == "bls":
        isBls = 1
    elif args.algo == "ls":
        isLs = 1
    elif args.algo == "plav":
        isPlav = 1
    else:
        raise ValueError("Invalid periodogram type: " + args.algo)

    # Retrieve data from fargs
    [ndata, time] = funcArgsGetTime(fargs)
    [ndata, mag] = funcArgsGetMag(fargs)
    [nsamp, period] = funcArgsGetPeriods(fargs)
    [nsamp, power] = funcArgsGetPower(fargs)

    # determine the width of each peak
    width = funcArgsGetPeakWidth(fargs, (isPlav or (USE_LOGNORMAL_BLS and isBls)))

    # skipme will be used to identify periods that are part of previously
    # identified peaks
    skipme = [0] * nsamp

    # sort by power and compute initial mean and deviation
    sortedPower = [0.0] * nsamp
    sortable = [[0.0, 0.0]] * nsamp

    for i in range(nsamp):
        p = 0
        if isLs:
            p = power[i]
        elif isBls:
            if USE_LOGNORMAL_BLS:
                if power[i] <= 0:
                    p = math.log(TINY_NUM)
                else:
                    p = math.log(power[i])
            else:
                p = power[i]
        elif isPlav:
            if power[i] <= 0:
                p = math.log(TINY_NUM)
            else:
                p = math.log(power[i])
        ##        print(p)
        ##        print(i)
        sortable[i] = [p, i]
    ##    print(sortable)
    sortable.sort()
    ##    print(sortable)
    for i in range(nsamp):
        sortedPower[i] = sortable[i][0]
    statNsamp, meanPwr, sdPwr = computePgramStats(args, nsamp, sortedPower, width)

    # bls-specific stuff
    H, L = 0.0, 0.0

    SDE = 0.0

    sig = 0
    j = 0
    nph = 0
    # print(sortable)
    while nph < args.nphased and sig <= args.sigThresh and j < nsamp:
        myIdxJ = sortable[nsamp - j - 1][1]
        # print(myIdxJ)
        if not isLs:
            # Compute normalized variate SDE and the probability of seeing
            # SDE by chance in the nsamp samples:
            if sdPwr > 0:
                SDE = (sortedPower[nsamp - j - 1] - meanPwr) / sdPwr
            if isBls:

                # prob that any given objs is < SDE
                prob = (1.0 + math.erf(SDE / math.sqrt(2.0))) / 2.0

                # prob that at least one obs in nsamp is >= SDE
                sig = 1.0 - pow(prob, statNsamp)
            else:
                # for plav assume gaussian distribution until
                # something better:
                # now using log normal (7/16/09)
                prob = (1.0 + math.erf(SDE / math.sqrt(2.0))) / 2.0
                sig = 1.0 - pow(prob, statNsamp)
        else:
            # the distribution for ls is exponential
            SDE = sortedPower[nsamp - j - 1]
            prob = 1.0 - math.exp(-SDE)
            sig = 1.0 - pow(prob, statNsamp)

        # If this is the bls algo, compute and save the values for H/L
        if (isBls and fargs.blsR and fargs.blsR[myIdxJ]):
            L = fargs.blsS[myIdxJ] / fargs.blsR[myIdxJ]
            H = fargs.blsS[myIdxJ] / (1 - fargs.blsR[myIdxJ])

        # Don't output the phased curve unless the significance is better than
        # the threshold
        if sig <= args.sigThresh:
            if skipme[myIdxJ]:
                pass
                # print("SKIPPING "+str(period[myIdxJ]))
            if not skipme[myIdxJ]:
                # print(str(period[myIdxJ])+" made it in!")
                # block out the other points associated with this peak
                for k in range(width[myIdxJ][0], width[myIdxJ][1]):
                    skipme[k] = 1
                if not 0 in skipme:
                    pass
                    # print("All periods skipped!")
                sigPeriod.append(period[myIdxJ])
                sigPower.append(power[myIdxJ])
            nph += 1
        else:
            args.nphased = nph
            # print("Stopping with nphased="+str(nph))
        j += 1

    # Get the output filename for the peaks
    outName = args.getOutputFile() + ".top"

    # Open the file
    f = open(outName, 'w')
    if args.outLabeled:
        saveLabeledOutput(f, sigPeriod, sigPower, "PERIOD", "POWER", args.inputDelimiter)
    else:
        saveOutput(f, sigPeriod, sigPower, args.inputDelimiter)
    f.close()


# computePgramStats()
# Function to compute the statistics to use for p-values
# Several options:
# - TRIM_DISTRIBUTION: remove "peaks" that are
#   part of the same peak. define being in the same peak as failing to
#   pass below the mean power between periods
# - TRIM_OUTLIERS: trim powers that are "outliers" before computing stats
def computePgramStats(args, nsamp, sortedPower, peakWidth):
    i = 0
    sumPwr = 0.0
    sumSqPwr = 0.0
    myP = sortedPower
    numP = nsamp

    # Make return variables
    statNsamp = None
    meanPwr = None
    sdPwr = None

    # if we have the input values for number of samples, mean and dev
    # save them in the local variables
    if args.powN != DEFAULT_POW_NUM:
        statNsamp = args.powN
    else:
        args.powN = nsamp
        statNsamp = nsamp
    if args.powMean != DEFAULT_POW_MEAN:
        meanPwr = args.powMean
    if args.powSd != DEFAULT_POW_SD:
        sdPwr = args.powSd

    # If all the values were input, just return
    if args.powN != DEFAULT_POW_NUM and \
            args.powMean != DEFAULT_POW_MEAN and \
            args.powSd != DEFAULT_POW_SD:
        return statNsamp, meanPwr, sdPwr

    # If this is lomb-scargle, we don't want to comput stats
    # and we set statNsamp above so just return
    if args.algo == "ls":
        return statNsamp, meanPwr, sdPwr

    # otherwise, compute stats
    if TRIM_DISTRIBUTION:
        # filter so we don't have multiple points representing
        # each peak

        trimmedP = [0.0] * nsamp
        numP = 0
        while i < nsamp:
            p = myP[i]
            trimmedP[numP] = p
            numP += 1
            sumPwr += p
            sumSqPwr += (p * p)
            i = peakWidth[i][1]  # index of the start of the next peak
        mn = sumPwr / numP
        if numP > 1:
            diff = sumSqPwr - (sumPwr * sumPwr) / numP
            if diff < 0:
                if diff > TINY_NUM:
                    raise RuntimeError("Calculation error in std dev")
                diff = 0
            sd = math.sqrt(diff / (numP - 1))
        else:
            sd = DEFAULT_POW_SD
        myP = trimmedP
    else:
        for i in range(numP):
            p = myP[i]
            sumPwr += p
            sumSqPwr += (p * p)
        mn = sumPwr / numP

        if numP > 1:
            diff = sumSqPwr - (sumPwr ** 2) / numP
            if diff < 0:
                if diff > TINY_NUM:
                    raise RuntimeError("Calculation error in std dev")
                diff = 0
            sd = math.sqrt(diff / (numP - 1))
        else:
            sd = DEFAULT_POW_SD
        if TRIM_OUTLIERS:
            # Trim outliers and re-compute mean and deviation
            # trimming gaussian outliers is not appropriate for the ls algo
            # because the distribution is exponential
            tnsamp, tpower, tmean, tsd = statTrimOutliers(numP, myP, NORMAL_OUTLIER, MAX_TRIM_FRACTION)
            mn = tmean
            sd = tsd
        # save computed values for output
        if args.powN == DEFAULT_POW_NUM:
            args.powN = nsamp
        if args.powMean == DEFAULT_POW_MEAN:
            args.powMean = mn
        if args.powSd == DEFAULT_POW_SD:
            args.powSd = sd
        statNsamp = args.powN
        meanPwr = args.powMean
        sdPwr = args.powSd
        return statNsamp, meanPwr, sdPwr


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


# New version of mod function, once I remembered
# that python has a builtin
def mod(n, m):
    return n % m


# Finds the smallest nonnegative
# number that is a multiple of m away from n.
# Recursive, for fun.
def mod_OLD(n, m):
    if m <= 0:
        raise ValueError("Inappropriate value: m must be positive!")
    if 0 <= n < m:
        return n
    if n < 0:
        return mod(n + m, m)
    else:
        return mod(n - m, m)


#######################
##
##  Classes
##
#######################

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

    # Function to parse arguments passed to the "periodogram"
    # command at the command line
    def populate(this):
        nbSet = 0;
        qminSet = 0;
        qmaxSet = 0;
        nProcSet = 0
        noutSet = 0;
        bsmSet = 0;
        stepSet = 0;
        dfSet = 0;
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
                this.algo = a[1]
            elif a[0] == '-b':
                # check for invalid number of bins
                if int(a[1]) < 5 or int(a[1]) > 2000:
                    raise ValueError("Inappropriate value for option -b:" + a[1] + "\nargument must be in [5,2000]!")
                # set the number of bins
                this.nbins = int(a[1])
                nbSet = 1
            elif a[0] == '-c':
                # not doing paralell processing right now.
                pass
            elif a[0] == '-d':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -d:" + a[1] + "\nargument must be greater than 0!")
                # set the frequency step
                this.dfreq = float(a[1])
                dfSet = 1
            elif a[0] == '-D':
                # Set the debugfile
                this.debugfile = a[1]
                if a[1] == 'stdout':
                    this.debugfp = sys.stdout
                elif a[1] == 'stderr':
                    this.debugfp = sys.stderr
                else:
                    try:
                        this.debugfp = open(a[1], 'w')
                    except:
                        this.debugfp = sys.stderr
                        this.debugfile = 'stderr'
                        print("Unable to open file " + a[1] + ", using stderr.")  # Maybe print this to stderr?
            elif a[0] == '-e':
                # set the input delimiter
                this.inputDelimiter = a[1]
            elif a[0] == '-f':
                # check for conflict, either in value or with maxperiod
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -f:" + a[1] + "\nargument must be greater than 0!")
                if this.maxperiod != DEFAULT_MAXPERIOD:
                    print("Error. Cannot specify both -f and -P! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set maxperiod and asFreq
                this.maxperiod = 1.0 / float(a[1])
                this.asFreq = 1
            elif a[0] == '-F':
                # check for conflict, either in value or with minperiod
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -F:" + a[1] + "\nargument must be greater than 0!")
                if this.minperiod != DEFAULT_MINPERIOD:
                    print("Error. Cannot specify both -F and -p! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set minperiod and asFreq
                this.minperiod = 1.0 / float(a[1])
                this.asFreq = 1
            elif a[0] == '-h':
                # set datahome
                this.datahome = a[1]
            elif a[0] == '-H':
                # Unused
                this.hdu = int(a[1])
                if this.hdu != DEFAULT_HDU and \
                        this.hdu <= 0:
                    raise ValueError("Inappropriate value for option -H:" + a[1] + "\nargument must begreater than 0!")
            elif a[0] == '-i':
                # check for pstepType validity
                if a[1] not in ["exp", "std", 'plav', 'fixedf', 'fixedp']:
                    raise ValueError("Inappropriate value for option -i:" + a[
                        1] + "\nargument must be one of exp, std, plav, or fixedf!")
                # set pstepType
                this.pstepType = a[1]
            elif a[0] == '-I':
                # Implement later?
                pass
            elif a[0] == '-K':
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -K:" + a[1] + "\nargument must be greater than 0!")
                # set stat number of samples
                this.powN = int(a[1])
            elif a[0] == '-l':
                # For server, not implemented.
                pass
            elif a[0] == '-L':
                # must be boolean flag, check.
                if a[1] not in ["0", "1"]:
                    raise ValueError("Inappropriate value for option -L:" + a[1] + "\n argument must be either 1 or 0!")
                # set outLabeled
                this.outLabeled = bool(int(a[1]))
            elif a[0] == '-M':
                # set statmean
                this.powMean = float(a[1])
            elif a[0] == '-n':
                # check that we don't have a negative number of outliers
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -n:" + a[1] + "\nargument must be greater than 0!")
                # set nout and noutSet
                this.nout = int(a[1])
                noutSet = 1
            elif a[0] == '-N':
                # we can't output a negative number of curves!
                if int(a[1]) < 0:
                    raise ValueError("Inappropriate value for option -N:" + a[1] + "\nargument must be greater than 0!")
                # set number of curves to output
                this.nphased = int(a[1])
            elif a[0] == '-o':
                if int(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -o:" + a[1] + "\nargument must be greater than 0!")
                # set the oversample factor
                this.oversample = int(a[1])
            elif a[0] == '-O':
                # Implement later?
                pass
            elif a[0] == '-p':
                # check for conflict, in value or with maxfrequency
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -p:" + a[1] + "\nargument must be greater than 0!")
                if this.minperiod != DEFAULT_MINPERIOD:
                    print("Error. Cannot specify both -F and -p! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set minperiod
                this.minperiod = float(a[1])
            elif a[0] == '-P':
                # check for conflict, in value or with minfrequency
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -P:" + a[1] + "\nargument must be greater than 0!")
                if this.maxperiod != DEFAULT_MAXPERIOD:
                    print("Error. Cannot specify both -f and -P! See usage below.")
                    print(PGRAM_USAGE_TEXT)
                    exit()
                # set maxperiod
                this.maxperiod = float(a[1])
            elif a[0] == '-q':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -q:" + a[1] + "\nargument must be greater than 0!")
                # save the minimum fraction of a period in transit
                this.qmin = float(a[1])
                qminSet = 1
            elif a[0] == '-Q':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -Q:" + a[1] + "\nargument must be greater than 0!")
                # save the maximum fraction of a period in transit
                this.qmax = float(a[1])
                qmaxSet = 1
            elif a[0] == '-r':
                # remote server stuff. not used.
                pass
            elif a[0] == '-R':
                # set the output directory
                this.outDir = a[1]
            elif a[0] == '-s':
                if float(a[1]) <= 0 or float(a[1]) >= 1:
                    raise ValueError("Inappropriate value for option -s:" + a[1] + "\nargument must be in range (0,1)!")
                # set the phase-smoothing box size
                this.smooth = float(a[1])
                bsmSet = 1
            elif a[0] == '-S':
                if float(a[1]) <= 0 or float(a[1]) > 1:
                    raise ValueError("Inappropriate value for option -S:" + a[1] + "\nargument must be in range (0,1]!")
                # set the significance threshold for output of phased curves
                this.sigThresh = float(a[1])
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
                this.substep = float(a[1])
            elif a[0] == '-V':
                if float(a[1]) <= 0:
                    raise ValueError("Inappropriate value for option -V:" + a[1] + "\nargument must be greater than 0!")
                # set the standard deviation for p-values
                this.powSd = float(a[1])
            elif a[0] == '-w':
                # set the constraintMin
                this.constraintMin = float(a[1])
            elif a[0] == '-W':
                # set the constraintMax
                this.constraintMax = float(a[1])
            elif a[0] == '-x':
                # check that the length of the string it was passed is
                # greater than 0
                if len(a[1]) > 0:
                    # set the xcol name
                    this.xcol = a[1]
            elif a[0] == '-X':
                if len(a[1]) > 0:
                    # set the constraintCol name
                    this.constraintCol = a[1]
            elif a[0] == '-y':
                if len(a[1]) > 0:
                    # set the ycol name
                    this.ycol = a[1]
            elif a[0] == '-Y':
                if len(a[1]) > 0:
                    # set the yerrCol name
                    this.yerrCol = a[1]
            else:
                print('Option not recognised. Somehow you bypassed the GetOptError?')
                print(PGRAM_USAGE_TEXT)

        # check mins/maxes are in correct relative positions (min>max)
        if this.qmin >= this.qmax:
            raise ValueError("Error: FractionOfPeriodInTransitMax must be greater than FractionOfPeriodInTransitMin")
        if this.minperiod != DEFAULT_MINPERIOD and \
                this.maxperiod != DEFAULT_MAXPERIOD and \
                this.minperiod >= this.maxperiod:
            if this.asFreq == 1:
                raise ValueError("Error: FrequencyRangeMax must be greater than FrequencyRangeMin!")
            else:
                raise ValueError("Error: PeriodRangeMax must be greater than PeriodRangeMin!")

        # Check that args haven't been set inconsistently
        if nbSet and this.algo != "bls":
            print("Error. Option -b cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if qminSet and this.algo != "bls":
            print("Error. Option -q cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if qmaxSet and this.algo != "bls":
            print("Error. Option -Q cannot be used without -a bls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if dfSet and "fixedf" != this.pstepType != "fixedp":
            print("Error. Option -d cannot be used with this PeriodStepMethod.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if noutSet and this.algo != "plav":
            print("Error. Option -n cannot be used without -a plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if stepSet and this.pstepType != "plav":
            print("Error. Option -u cannot be used without -i plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if this.pstepType == "plav" and this.oversample != DEFAULT_OVERSAMPLE:
            print("Error. Option -o cannot be used with -i plav.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if dfSet and this.oversample != DEFAULT_OVERSAMPLE:
            print("Error. Option -d cannot be used with -o.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if this.powMean != DEFAULT_POW_MEAN and \
                this.algo == 'ls':
            print("Error. Option -M cannot be used with -a ls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if this.powSd != DEFAULT_POW_SD and \
                this.algo == 'ls':
            print("Error. Option -V cannot be used with -a ls.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if (this.constraintCol == None or this.constraintCol == "none") and \
                (this.constraintMin != DEFAULT_CONSTRAINT or \
                 this.constraintMax != DEFAULT_CONSTRAINT):
            print("Error. Cannot set constraint extremes without column.")
            print(PGRAM_USAGE_TEXT)
            exit()
        if (this.constraintCol != None or this.yerrCol != None or \
            this.ycol != None or this.xcol != None) and \
                (this.ycol == None or this.xcol == None):
            print("Error. TimeColumn and DataColumn must be set if any columns are set!")
            print(PGRAM_HELP_TEXT)
            exit()

        # read the required input file and  optional output file
        this.intbl = None
        this.outtbl = None
        for s in args:
            # check that we're not misinterpreting a switch, or trying
            # to open "." or ".."
            if s[0] == '-':
                raise ValueError("Misplaced argument " + s + ": switches must precede arguments!")
            elif s[0] == '.' and (len(s) == 1 or (s[1] == '.' and len(s) == 2)):
                if this.intbl == None:
                    raise ValueError("InputFile: " + s + " is not a valid value.")
                else:
                    raise ValueError("OutputFile: " + s + " is not a valid value.")
            else:
                if this.intbl == None:
                    this.intbl = s
                elif this.outtbl == None:
                    this.outtbl = s
                else:
                    raise ValueError("Argument " + s + " not recognized.")
        if this.intbl == None:
            raise ValueError("InputFile: must supply an input file!")
        if this.intbl[0] != '/':
            if this.datahome != None:
                if this.datahome[len(this.datahome) - 1] == '/':
                    this.datahome = this.datahome[:len(this.datahome) - 1]
                this.intbl = this.datahome + '/' + this.intbl
            else:
                this.intbl = os.getcwd() + '/' + this.intbl
        this.getInBase()
        if len(this.inBase) == 0:
            raise ValueError("InputFile: '" + this.intbl + "' does not contain a file name!")
        this.getOutputFile()
        this.getOutBase()

    # Function to get the pathless input filename
    def getInBase(this):
        if this.inBase == None:
            tmp = this.intbl.split("/")
            tmp = tmp[len(tmp) - 1]
            this.inBase = tmp
        return this.inBase

    # Function to get the pathless output filename, creates
    # one if it is not yet set
    def getOutBase(this):
        if this.outBase == None:
            if this.outtbl == None:
                this.outBase = this.getInBase() + ".out"
            else:
                tmp = this.outtbl.split("/")
                tmp = tmp[len(tmp) - 1]
                this.outBase = tmp
        return this.outBase

    # Function to get the complete absolute output filename,
    # sets it if it is unset
    def getOutputFile(this):
        path = os.getcwd()
        if this.outDir != None:
            if this.outDir[0] != '/':
                path = os.getcwd() + '/' + this.outDir
                this.outDir = path
            else:
                path = this.outDir
        if path[len(path) - 1] == '/':
            path = path[:len(path) - 1]
        if this.outtbl == None:
            tmpFname = path + '/' + this.getOutBase()
            if not tmpFname:
                raise RuntimeError("Error. Outtbl filename empty!")
            this.outtbl = tmpFname
        else:
            if this.outtbl[0] != '/':
                tmpFname = path + '/' + this.outtbl
                this.outtbl = tmpFname
        if this.outtbl == None:
            print("path:" + path)
            print("outdir:" + this.outDir)
        return this.outtbl

    def argsPrint(this):
        arr = ['python3', 'mp_alt.py']
        if this.algo:
            arr.append('-a')
            arr.append(str(this.algo))
        if this.nbins != DEFAULT_NBINS:
            arr.append('-b')
            arr.append(str(this.nbins))
        if this.dfreq != DEFAULT_DFREQ and this.pstepType in ["fixedf", "fixedp"]:
            arr.append('-d')
            arr.append(str(this.dfreq))
        if this.inputDelimiter:
            arr.append('-e')
            arr.append("'" + str(this.inputDelimiter) + "'")
        if this.maxperiod != DEFAULT_MAXPERIOD and this.asFreq:
            arr.append('-f')
            arr.append(str(1.0 / this.maxperiod))
        if this.minperiod != DEFAULT_MINPERIOD and this.asFreq:
            arr.append('-F')
            arr.append(str(1.0 / this.minperiod))
        ##        if this.datahome:
        ##            arr.append('-h')
        ##            arr.append(str(this.datahome))
        if this.pstepType:
            arr.append('-i')
            arr.append(str(this.pstepType))
        if this.powN != DEFAULT_POW_NUM:
            arr.append('-K')
            arr.append(str(this.powN))
        if this.outLabeled:
            arr.append('-L')
            arr.append(str(this.outLabeled))
        if this.powMean != DEFAULT_POW_MEAN:
            arr.append('-M')
            arr.append(str(this.powMean))
        if this.nout and this.algo == "plav":
            arr.append('-n')
            arr.append(str(this.nout))
        if this.nphased:
            arr.append('-N')
            arr.append(str(this.nphased))
        if this.oversample != DEFAULT_OVERSAMPLE:
            arr.append('-o')
            arr.append(str(this.oversample))
        if this.minperiod != DEFAULT_MINPERIOD and not this.asFreq:
            arr.append('-p')
            arr.append(str(this.minperiod))
        if this.maxperiod != DEFAULT_MAXPERIOD and not this.asFreq:
            arr.append('-P')
            arr.append(str(this.maxperiod))
        if this.qmin and this.algo == "bls":
            arr.append('-q')
            arr.append(str(this.qmin))
        if this.qmax and this.algo == "bls":
            arr.append('-Q')
            arr.append(str(this.qmax))
        # removed because outDir is included in outtbl
        ##        if this.outDir:
        ##            arr.append('-R')
        ##            arr.append(str(this.outDir))
        if this.smooth:
            arr.append('-s')
            arr.append(str(this.smooth))
        if this.sigThresh:
            arr.append('-S')
            arr.append(str(this.sigThresh))
        #        if this.title:
        #            arr.append('-T')
        #            arr.append(str(this.title))
        if this.substep and this.pstepType == "plav":
            arr.append('-u')
            arr.append(str(this.substep))
        if this.powSd != DEFAULT_POW_SD:
            arr.append('-V')
            arr.append(str(this.powSd))
        if this.constraintMin != DEFAULT_CONSTRAINT:
            arr.append('-w')
            arr.append(str(this.constraintMin))
        if this.constraintMax != DEFAULT_CONSTRAINT:
            arr.append('-W')
            arr.append(str(this.constraintMax))
        if this.xcol:
            arr.append('-x')
            arr.append(str(this.xcol))
        if this.constraintCol:
            arr.append('-X')
            arr.append(str(this.constraintCol))
        if this.ycol:
            arr.append('-y')
            arr.append(str(this.ycol))
        if this.yerrCol:
            arr.append('-Y')
            arr.append(str(this.yerrCol))
        arr.append(this.intbl)
        arr.append(this.outtbl)

        # len(arr) must be at least two, so it's safe to
        # acess the first element without a check
        string = arr[0]
        for i in range(1, len(arr)):
            string += ' ' + arr[i]
        return string


# Class to hold the arguments that will be passed in
# and out of a periodogram function call. Will also be
# used to store utility arrays
class funcArgs:
    def __init__(this):
        # populated by populateLite() before calling
        # specific algorithms
        this.ndata = 0
        this.time = None
        this.mag = None

        # these arguments are set based on pstepType and
        # associated command-line arguments
        this.nsamp = 0
        this.period = None

        this.timeEst = -1

        # populated by the specific algorithm
        this.power = None

        # start/end index of each peak
        this.width = None

        # BLS results... incase anyone wants to know where the phase
        # bounds of the transit were
        this.blsS = None
        this.blsR = None
        this.lowBin0 = None
        this.lowBin1 = None

        # From input arguments
        this.boxSize = None

        # reusable arrays, recomputed during
        # each call to phaseLightCurve based
        # on the value of p
        this.p = 0

        # the following are sorted by phase
        this.phase = None
        this.phasedMag = None
        this.smoothedMag = None
        this.chi = None  # for use by plav

        # sortable will hold the original order of the last array sorted
        this.sortable = None

    # populateLite()
    # Function to import time and mag (possibly adjusted by minDay and
    # meanMag) into the funcArgs object. This will be called by populate()
    def populateLite(this, algo, data):

        # Retrieve data for processing: number of points, day and mag
        time = dtGetFilteredArray(data, DATA_FIELD_TYPE.DATA_X)
        mag = dtGetFilteredArray(data, DATA_FIELD_TYPE.DATA_Y)
        if len(time) != len(mag):
            raise RuntimeError("Houston, we have a problem. len(mag)!=len(time).")
        ndata = len(time)

        # Adjust the local copy of the data as appropirate for the different
        # algorithms
        meanMag = 0.0;
        minDay = 0.0
        adjustByMean = 0
        adjustByMinDay = 0

        if algo == "ls":
            adjustByMean = 1
        elif algo == "bls":
            adjustByMean = 1
            adjustByMinDay = 1
        elif algo == "plav":
            pass
        else:
            raise ValueError("Inappropriate value for algo:\n" + algo)

        if adjustByMean or adjustByMinDay:
            if adjustByMean:
                meanMag = dtGetMean(data, DATA_FIELD_TYPE.DATA_Y)
            if adjustByMinDay:
                minDay = dtGetMin(data, DATA_FIELD_TYPE.DATA_X)
            for i in range(ndata):
                if adjustByMean:
                    mag[i] -= meanMag
                if adjustByMinDay:
                    time[i] -= minDay

        this.ndata = ndata
        this.time = time
        this.mag = mag

    # populate()

    # Function to import time and mag (possibly adjusted by minDay and
    # meanMag) into the funcArgs object. In addition, input arguments
    # will be used to set the periods at which we will compute power in
    # computePeriodogram(). Calls populateLite()
    def populate(this, args, data):
        # set algo for populateLite (could just feed it args.algo,
        # but it's nice to have a local version.)
        algo = args.algo

        # populate and adjust time and mag
        this.populateLite(algo, data)
        ndata = this.ndata
        time = this.time
        mag = this.mag

        dell = TINY_NUM

        # We have set an extremely mild restriction: require at least 2
        # points in the file
        if ndata < MIN_NDATA:
            raise ValueError("InputFile: Not enough data in file to process: " + str(ndata))

        # Get minimum and maximum time values, along with the smallest
        # time-difference between two data points
        minDay = dtGetMin(data, DATA_FIELD_TYPE.DATA_X)
        maxDay = dtGetMax(data, DATA_FIELD_TYPE.DATA_X)
        minDt = dtGetMinDiff(data, DATA_FIELD_TYPE.DATA_X)

        # run a sanity check on the data
        if minDay >= maxDay:
            raise ValueError("InputFile: minimum time is >= maximum time!")

        # compute min and max periods from data if not set from cmdline
        minperiod = args.minperiod
        maxperiod = args.maxperiod
        if maxperiod == DEFAULT_MAXPERIOD:
            maxperiod = maxDay - minDay
        # If we want to restrict to periods that we could actually observe
        # in their entirety
        elif RESTRICT_TO_COMPLETELY_OBSERVABLE:
            if maxperiod > (maxDay - minDay):
                raise ValueError("PeriodRangeMax: Periods should not exceed time spanned " + str(maxDay - minDay))

        if minperiod == DEFAULT_MINPERIOD:
            # median time step
            minperiod = dtGetMedianDiff(data, DATA_FIELD_TYPE.DATA_X)

            if AVG_TIME_STEP:
                # average time step
                minperiod = (maxDay - minDay) / ndata

            if SMALLEST_TIME_STEP:
                # smallest time step
                if minDt == 0:
                    print("ERROR: Min DT = 0... dtGetMinDiff error bypass?")
                    minperiod = maxperiod / ndata
                else:
                    minperiod = minDt

            if minperiod < MIN_MINPERIOD:
                minperiod = MIN_MINPERIOD
        # check that max is greater than min
        if minperiod > maxperiod or minperiod <= 0 or maxperiod <= 0:
            if args.asFreq:
                raise ValueError("FrequencyRangeMin: frequency range error: from " + str(1 / maxperiod) + " to " + str(
                    1 / minperiod))
            else:
                raise ValueError("PeriodRangeMin: period range error: from " + str(minperiod) + " to " + str(maxperiod))

        # Finally, set minperiod and maxperiod
        args.minperiod = minperiod
        args.maxperiod = maxperiod

        # Period stepping

        # What type of steppin are we doing and how many periods
        # will we be considering?
        ptype = args.pstepType
        oversample = args.oversample

        period = []
        count = 0

        # If we're not reading the periods from a file, generate them
        if args.periodFile == None:

            if ptype == "std":

                # Standard period stepping, the array length is
                # <=ndata*oversample
                nsamp = ndata * oversample

                if ADJUST_PERIOD_TO_RANGE:
                    offset = minperiod - 1.0 / ((1.0 + (nsamp - 1)) / nsamp)
                else:
                    offset = 0

                for i in range(nsamp):
                    w = 2.0 * math.pi * ((1.0 + (nsamp - 1 - i)) / nsamp)
                    p = 2.0 * math.pi / w + offset
                    if p > maxperiod:
                        break
                    elif (p >= minperiod):
                        period.append(p)
                        count += 1
                    elif ADJUST_PERIOD_TO_RANGE:
                        raise RuntimeError("PeriodStepMethod: Still skipping period " + str(p))
                nsamp = count

            elif ptype == "exp":

                # Peter's exponential period stepping: the array length
                # is <=ndata*oversample

                nsamp = ndata * oversample

                pdMagSpan = math.ceil(math.log10(maxperiod / minperiod))
                if ADJUST_PERIOD_TO_RANGE:
                    offset = minperiod - (10 ** (-pdMagSpan)) * maxperiod
                else:
                    offset = 0

                # this should be done with now: output of version that checked
                # out without this now saved for future reference
                if False:
                    # for comparison with Peter' output
                    if algo == "bls":
                        pdMagSpan = 4
                    elif algo == "ls":
                        pdMagSpan = 6

                for i in range(nsamp):
                    p = (10 ** (pdMagSpan * ((1.0 * i / nsamp) - 1))) * maxperiod + offset
                    if p > maxperiod:
                        break
                    elif p >= minperiod:
                        period.append(p)
                        count += 1
                    elif ADJUST_PERIOD_TO_RANGE:
                        raise RuntimeError("PeriodStepMethod: Still skipping period " + str(p))
                    if count != len(period):
                        print("Error: period length!=count!")
                        print("Period length: " + str(len(period)))
                        print("Count: " + str(count))
                        exit()

                nsamp = count

            elif ptype == "fixedf":

                # fixed frequency stepping: this is what the next method
                # alleges to be...?

                minf = 1.0 / maxperiod
                maxf = 1.0 / minperiod

                if args.dfreq == DEFAULT_DFREQ:
                    nsamp = int(ndata * oversample)
                    df = (maxf - minf) / (nsamp - 1)
                else:
                    df = args.dfreq
                    nsamp = int((maxf - minf) / df) + 1

                args.dfreq = df
                f = maxf

                for i in range(nsamp):
                    period.append(1.0 / f)
                    if (period[i] > (maxperiod + dell)) or \
                            (period[i] < (minperiod - dell)):
                        raise RuntimeError("PeriodStepMethod: bad period in fixedf: " + str(period[i]) + " (" + str(
                            minperiod) + " to " + str(maxperiod) + ")")
                    f -= df

            elif ptype == "fixedp":

                # fixed period stepping

                if args.dfreq == DEFAULT_DFREQ:
                    nsamp = int(ndata * oversample)
                    dp = (maxperiod - minperiod) / (nsamp - 1)
                else:
                    dp = args.dfreq
                    nsamp = int((maxperiod - minperiod) / dp) + 1

                p = minperiod
                for i in range(nsamp):
                    period.append(p)
                    if (period[i] > maxperiod + dell) or \
                            (period[i] < minperiod - dell):
                        raise RuntimeError("PeriodStepMethod: bad period in fixedp: " + str(period[i]) + " (" + str(
                            minperiod) + " to " + str(maxperiod) + ")")
                    p += dp

            elif ptype == "plav":

                # size of the period step is based on the current period.
                # Peter says it was suggested by someone on his committee
                # whose opinion he trusts.

                pstep = args.substep / (maxDay - minDay)
                p = minperiod

                while p < maxperiod:
                    period.append(p)
                    count += 1
                    p += pstep * p * p

                nsamp = count

            else:
                raise ValueError("PeriodStepMethod: Invalid period stepping method!")
            if nsamp == 0:
                raise ValueError(
                    "PeriodRangeMin: No periods satisfy min/max period constraints (" + str(minperiod) + "-" + str(
                        maxperiod) + ") with this PeriodStepMethod!")
        else:
            # read periods to process from a file
            fp = open(args.periodFile, 'r')
            lines = fp.readlines()
            fp.close()
            for l in lines:
                period.append(float(l.strip()))
            nsamp = len(period)

        # save the data into the fargs structure
        this.boxSize = args.smooth
        this.ndata = ndata
        this.time = time
        this.mag = mag
        this.nsamp = nsamp
        this.period = period
        this.power = [0.0] * len(period)
        # Initialize chi for use in computePlavchan
        if algo == 'plav':
            this.chi = [0.0] * this.ndata

        # Initialize other arrays
        this.phase = [0.0] * this.ndata
        this.phasedMag = [0.0] * this.ndata
        this.smoothedMag = [0.0] * this.ndata

        # compute time estimate
        this.timeEst = estimateProcessingTime(nsamp, ndata, args, args.qmax, algo)
        if nsamp != len(this.period):
            raise RuntimeError("Error: nsamp!=len(period)")
        # print(period)


# Class used to hold data from files, column names (labels),
# filtering info, and various statistics about each set of
# data. It is organized according to DATA_FIELD_TYPE, an enum
# that dictates which indicies denote which type of data (xdata,
# ydata,yUncertainty, etc)
class dataTbl:
    def __init__(this):
        this.description = None  # UNUSED
        this.origin = None  # UNUSED

        this.ndata = 0  # number of data points
        this.ndataUnfilt = 0  # number of data points not filtered
        this.dataArray = []  # each index contains a dataArray of
        # length ndata corresponding to the
        # DATA_FIELD_TYPE of that index
        this.colKeyNames = []  # UNUSED
        this.arrayNames = []  # contains the column names associated
        # with the DATA_FIELD_TYPE of that index
        this.units = []  # UNUSED
        this.median = []  # the median for each array
        this.mean = []  # the mean for each array
        this.sd = []  # the std dev for each array
        for i in range(DATA_FIELD_TYPE.DATA_N_TYPES):
            this.colKeyNames.append(None)
            this.arrayNames.append(None)
            this.dataArray.append(None)
            this.units.append(None)
            this.median.append(UNSET_MEAN)
            this.mean.append(UNSET_MEAN)
            this.sd.append(UNSET_MEAN)

        this.omitRow = None  # array containing flags to skip a row.
        # length:ndata
        this.tinfo = None  # UNUSED
        this.nhead = 0  # UNUSED
        this.header = 0  # UNUSED

    # Function to populate the dataTbl. Calls readDoubleColumns to
    # get data from file
    def populate(this, intbl, xcol, ycol, yerrCol, constraintCol, \
                 constraintMin, constraintMax, delim):
        # set the column names
        this.setColName(DATA_FIELD_TYPE.DATA_X, xcol)
        this.setColName(DATA_FIELD_TYPE.DATA_Y, ycol)
        this.setColName(DATA_FIELD_TYPE.DATA_CONSTRAINT, constraintCol)
        this.setColName(DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY, yerrCol)

        # Even though we don't use it, let's set it anyway.
        this.origin = intbl

        # Index by which we sort data will be the xdata
        sortIdx = DATA_FIELD_TYPE.DATA_X

        # Read the data, omitting rows that contain non-numerical data
        unsortedData = readDoubleColumns(intbl, this.arrayNames, delim)
        # Sort the data
        this.dataArray = sortColumns(unsortedData, sortIdx)
        this.ndata = len(this.dataArray[DATA_FIELD_TYPE.DATA_X])
        if this.ndata <= 0:
            raise ValueError("InputFile: No data found in file!")

        # Set an initial value for omitRow
        this.omitRow = [None] * this.ndata

        # print(this.omitRow==None)
        # Filter bad data
        this.setFilter(constraintMin, constraintMax)

    # Function to set the arrayNames entry to the specified value
    def setColName(this, f, s):
        if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
            raise ValueError('Inappropriate value:\nf is greater than maximum value!')
        if s:
            this.arrayNames[f] = s

    # Function to set omitRow flags based on values in the
    # constraint (aka filter) column
    def setFilter(this, minVal, maxVal):
        # print(this.omitRow==None)
        # Initially, nothing filtered
        this.ndataUnfilt = this.ndata

        # check that filtering is applicable here
        if this.arrayNames[DATA_FIELD_TYPE.DATA_CONSTRAINT] != None and \
                this.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT] != None:

            # have we set lower and upper limits?
            if minVal == UNSET_VALUE or maxVal == UNSET_VALUE:
                tmp = this.arrayNames[DATA_FIELD_TYPE.DATA_CONSTRAINT].lower()
                # status!=0 signifies error - do not allow anything else
                if "status" in tmp:
                    minVal = 0
                    maxVal = 0
                else:
                    # identify the minimum and maximum values -- not constrained
                    # but informational
                    myArray = this.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT]
                    minVal = min(myArray)
                    maxVal = max(myArray)
                    return None  # no restrictions set -> nothing to do
            # if we get here, we have a constraint field and min and max values
            nfilt = 0
            for i in range(this.ndata):
                # if anything violates the max/min values
                if not this.omitRow[i] and \
                        (this.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT][i] < minVal or \
                         this.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT][i] > maxVal):
                    # exclude it
                    this.omitRow[i] = 1
                    nfilt += 1
            if nfilt == this.ndata:
                raise ValueError("InputFile: all data filtered by constraints!")
            this.ndataUnfilt = this.ndata - nfilt


# This class determines which index for every column of dataArray
# signifies each type of data
class DATA_FIELD_TYPE(IntEnum):
    DATA_X = 0
    DATA_X_UNCERTAINTY = 1
    DATA_Y = 2
    DATA_Y_UNCERTAINTY = 3
    DATA_CONSTRAINT = 4
    DATA_N_TYPES = 5


#######################
##
##  Help texts
##
#######################


PGRAM_USAGE_TEXT = " Usage: periodogram \n\t[-a <PeriodogramType (algorithm): one of ls, bls, plav>]" \
                   "\n\t[-b <NumberOfBins (-a bls only)>]" \
                   "\n\t[-d <FixedStepSize (-i fixedf or -i fixedp only)>]" \
                   "\n\t[-D <output file name for debugging (may say 'stdout' or 'stderr')>]" \
                   "\n\t[-e <DataDelimiter>]" \
                   "\n\t[-f <FrequencyRangeMin> | -P <PeriodRangeMax>]" \
                   "\n\t[-F <FrequencyRangeMax> | -p <PeriodRangeMin>]" \
                   "\n\t[-h <DataHome directory: location of the input file>]" \
                   "\n\t[-i <PeriodStepMethod: std, exp, fixedf, fixedp, plav>]" \
                   "\n\t[-K <StatNumberOfSamples>]" \
                   "\n\t[-L <OutFileLabeled>]" \
                   "\n\t[-M <StatMean> (not with -a ls)]" \
                   "\n\t[-n <NumberOfOutliers (-a plav only)>]" \
                   "\n\t[-N <NumberOfPeaksToReturn>]" \
                   "\n\t[-o <OversampleFactor (not with -i plav)>]" \
                   "\n\t[-q <FractionOfPeriodInTransitMin (-a bls only)>]" \
                   "\n\t[-Q <FractionOfPeriodInTransitMax (-a bls only)>]" \
                   "\n\t[-R <OutputDirectory>]" \
                   "\n\t[-s <PhaseSmoothingBoxSize>]" \
                   "\n\t[-S <PeakSignificanceThreshold (on power for output)>]" \
                   "\n\t[-u <PeriodStepFactor (-i plav only)>]" \
                   "\n\t[-V <StatStandardDeviation> (not with -a ls)]" \
                   "\n\t[-w <ConstraintRangeMin>]" \
                   "\n\t[-W <ConstraintRangeMax>]" \
                   "\n\t[-x <TimeColumn>]" \
                   "\n\t[-X <ConstraintColumn>]" \
                   "\n\t[-y <DataColumn>]" \
                   "\n\t[-Y <DataErrorColumn>]" \
                   "\n\t<InputFile>" \
                   "\n\t[<OutputFile>]"
#                  "\n\t[-c <Number of processors to split job across>]"\
#                  "\n\t[-g <Remote server configuration file>]"\
#                  "\n\t[-H <FitsHeaderDataUnit to use (for input files in FITS format)>]"\
#                  "\n\t[-r <Remote server name>]"\
#                  "\n\t[-t <Remote server port number>]"\
#                  "\n\t[-T <Title (name of star)>]"\
#    "\n -T <Title>"                                                     \
#    "\n    The name of the star.  This will be used primarily for graphics" \

PGRAM_HELP_TEXT = "\n\n Description:  " \
                  "\n " \
                  "\n Compute a periodogram using one of three algorithms: Lomb-Scargle, " \
                  "\n BLS, or Plavchan." \
                  "\n " \
                  "\n Syntax: " \
                  "\n" + PGRAM_USAGE_TEXT + \
                  "\n Switches: " \
                  "\n " \
                  "\n --help" \
                  "\n    Returns this message.\n" \
                  "\n -a <PeriodogramType (algorithm)> " \
                  "\n    Specifies which algorithm to run (one of ls, bls, plav). " \
                  "\n -b <NumberOfBins>" \
                  "\n    Specifies the number of bins to use in the bls algorithm" \
                  "\n -d <FixedStepSize (-i fixedf or -i fixedp only)>" \
                  "\n    Specifies the size of the fixed frequency step or period step," \
                  "\n    depending on which of fixedf or fixedp was chosen for -i." \
                  "\n -e <DataDelimiter>" \
                  "\n    Delimiter separating data and (optionally) column labels in the" \
                  "\n    InputFile. If not supplied, ',' will be used. This will also be" \
                  "\n    used to separate the output data in the OutputFile. For information" \
                  "\n    on how to use escape characters from the command line (e.g. \\t), see" \
                  "\n    http://www.gnu.org/software/bash/manual/html_node/ANSI_002dC-Quoting.html" \
                  "\n -f <FrequencyRangeMin> | -P <PeriodRangeMax>" \
                  "\n    Maximum period to consider (may be optionally specified as minimum freq)" \
                  "\n -F <FrequencyRangeMax> | -p <PeriodRangeMin>" \
                  "\n    Minimum period to consider (may be specified as maximum freq)" \
                  "\n -h <DataHome>" \
                  "\n    The directory where the data is located. If not supplied," \
                  "\n    the location from which this is being run will be used." \
                  "\n -i <PeriodStepMethod>" \
                  "\n    Specifies which type of period stepping to use " \
                  "\n    (one of std, exp, fixedf, plav)" \
                  "\n -K <StatNumberOfSamples>" \
                  "\n    The number of samples to use for computation of p-values for output peaks." \
                  "\n    If not entered, the number of periods for which power is computed will " \
                  "\n    be used." \
                  "\n -L <OutFileLabeled>" \
                  "\n    Turns on and off column labels for the output file. Can be either" \
                  "\n    1 (on) or 0 (off). If not entered, it will default to 0, or off." \
                  "\n -M <StatMean>" \
                  "\n    Mean to use for computation of p-values for output peaks.  If not entered," \
                  "\n    the observed mean will be used." \
                  "\n -n <NumberOfOutliers>" \
                  "\n    Number of outliers to use in power calculation in the Plavchan algo" \
                  "\n -N <NumberOfPeaksToReturn>" \
                  "\n    Limit on the number of top peaks to output in table." \
                  "\n -o <OversampleFactor>]" \
                  "\n    Increase number of periods sampled by this factor (not for use" \
                  "\n    with -i plav or -d)" \
                  "\n -q <FractionOfPeriodInTransitMin>" \
                  "\n    Minimum fraction of period in transit to consider with BLS algo" \
                  "\n -Q <FractionOfPeriodInTransitMax>" \
                  "\n    Maximum fraction of period in transit to consider with BLS algo" \
                  "\n -R <OutputDirectory>" \
                  "\n    The directory in which to put output files (periodogram, " \
                  "\n    table of top periods).  The default is '.'" \
                  "\n -s <PhaseSmoothingBoxSize>" \
                  "\n    Size of box over which to average magnitudes for smoothed curve" \
                  "\n -S <PeakSignificanceThreshold>" \
                  "\n    Maximum p-value to accept for output peaks in the power spectrum" \
                  "\n -u <PeriodStepFactor>" \
                  "\n    Period increment factor for -i plav" \
                  "\n -V <StatStandardDeviation>" \
                  "\n    Standard deviation to use for computation of p-values for output peaks" \
                  "\n    If not entered, the observed standard deviation will be used" \
                  "\n -w <ConstraintRangeMin>" \
                  "\n    Smallest acceptable value of elements in ConstraintColumn" \
                  "\n -W <ConstraintRangeMax>" \
                  "\n    Largest acceptable value of elements in ConstraintColumn" \
                  "\n -x <TimeColumn>" \
                  "\n    Name of column in input file from which to read time info" \
                  "\n -X <ConstraintColumn>" \
                  "\n    Name of column in input file from which to read constraint info" \
                  "\n -y <DataColumn>" \
                  "\n    Name of column in input file from which to read measurement values" \
                  "\n -Y <DataErrorColumn>" \
                  "\n    Name of column in input file from which to read measurement errors" \
                  "\n " \
                  "\n Arguments: " \
                  "\n " \
                  "\n <input file> " \
                  "\n    Text file: first row contains either labels or data, separated by" \
                  "\n    the DataDelimiter. If it contains data, do not supply -x, -X, -y," \
                  "\n    or -Y. Otherwise, supply all that you want used. labeled column" \
                  "\n    exists, but its label is not associated with If a a data type" \
                  "\n    in the arguments, it will be ignored. Subsequent rows contain data," \
                  "\n    with different data types separated by the DataDelimiter." \
                  "\n    If column labels are not supplied, the following data types will" \
                  "\n    be assumed:" \
                  "\n              2 columns: TIME,DATA" \
                  "\n              3 columns: TIME,DATA,ERROR" \
                  "\n              4 columns: TIME,DATA,ERROR,CONSTRAINT" \
                  "\n              5+ columns: TIME,DATA,ERROR,CONSTRAINT,UNUSED,UNUSED..." \
                  "\n" \
                  "\n <output file> " \
                  "\n    [Optional] Text file. All content will be overwritten. Columns will" \
                  "\n    be labeled according to OutFileLabeled, with a column format of" \
                  "\n    PERIOD,POWER. If no file is specified, one will be constructed in the" \
                  "\n    output directory with the name '<path-free name of input file>.out" \
                  "\n " \
                  "\n Results: " \
                  "\n " \
                  "\n If successful, periodogram creates an output table file containing " \
                  "\n period and power, prints \"[struct stat=\"OK\", msg=\"<msg>\"]\" to stdout, " \
                  "\n and exits with 0.  The output message contains the command line arguments" \
                  "\n needed to replicate the exact results, including derived quantities if any." \
                  "\n " \
                  "\n Examples: " \
                  "\n " \
                  "\n The following example runs periodogram on a table file with the period " \
                  "\n range from .5 days to 1000 days and saves the output to out.tbl: " \
                  "\n " \
                  "\n $ periodogram test/test.txt -p .5 -P 1000 out" \
                  "\n "

#######################
##
##  Code that is run
##
#######################

if __name__ == '__main__':
    with mp.Pool(processes=mp.cpu_count()) as pool:
        # Construct and populate the pgramArgs object
        args = pgramArgs()
        args.populate()

        # Get the output file stream
        if args.outToStdOut:
            out = sys.stdout
        else:
            fname = args.getOutputFile()
            if fname == None:
                print(args.getInBase())
                print(args.getOutBase())
                print(args.getOutputFile())
            out = open(fname, 'w+')

        # Construct and populate the dataTbl object
        data = dataTbl()
        data.populate(args.intbl, args.xcol, args.ycol, args.yerrCol, \
                      args.constraintCol, args.constraintMin, args.constraintMax, args.inputDelimiter)

        # Construct and populate the funcArgs object
        fargs = funcArgs()
        fargs.populate(args, data)

        # Compute the periodogram
        computePeriodogram(args, data, fargs, pool)

        # Save the output
        if args.outLabeled == 1:
            saveLabeledOutput(out, fargs.period, fargs.power, "PERIOD", 'POWER', args.inputDelimiter)
        else:
            saveOutput(out, fargs.period, fargs.power, args.inputDelimiter)
        if not args.outToStdOut:
            out.close()
        print("All done!")
        print("Here is the command to reconstruct this:")
        print(args.argsPrint())

        # "stop" the timer
        t_f = time.time()
        t_d = t_f - t_i
        if PRINT_TIMES:
            # Print the result from our timer
            print("Time elapsed: " + format(t_d, '.4f') + " seconds (" + format(t_d / 60, '.4f') + " minutes)")

##def multiFunc(a,b,c):
##    return a/b,c
##if __name__=='__main__':
##    with mp.Pool(4) as pool:
##        results=pool.starmap(multiFunc,transpose([range(10),range(1,11),range(2,12)]))
