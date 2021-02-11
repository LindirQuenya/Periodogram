#!/usr/bin/env python3
# TODO: add comment

from func_args import *
from cli_args import *


# computePlavchan()
# Function to compute periodogram based on Plavchan 2008 algo
#
# ref: Peter Plavchan, M. Jura, J. Davy Kirkpatrick, Roc M. Cutri,
#     and S. C. Gallagher, "NEAR-INFRARED VARIABILITY IN THE 2MASS
#     CALIBRATION FIELDS: A SEARCH FOR PLANETARY TRANSIT CANDIDATES."
#     ApJS 175:191Y228 (2008)
#
# For each of a set of candidate periods, this algorithm folds a light
# curve to that period and then computes a "smoothed" curve by averaging
# the curve over a box spanning a certain phase range to either side (defined
# by the parameter "smooth").  The ratio of the sum of squared deviations
# from the mean (over the "nout" worst-fitting points) is divided by
# the sum of squared deviations from the smoothed values (again, over nout).
# The smaller the deviation from the smoothed curve, the larger this ratio
# will be, indicating that the smooth curve is a substantially better fit
# than the straight line "mag = mean mag".  This ratio is interpreted as the
# "power" at that period.
#
# Arguments:
#  data = populated dataTbl object
#  args = populated pgramArgs object
#  fargs = funcArgs object with ndata, time, mag
#               nsamp, and period set. power will
#               be populated by this function
#  pool = the pool object to give tasks to
def computePlavchan(data, args, fargs, pool):
    # Check for type errors
    if not isinstance(data, dataTbl):
        raise TypeError('Inappropriate argument type:\ndata must be of type dataTbl!')
    if not isinstance(fargs, funcArgs):
        raise TypeError('Inappropriate argument type:\nfargs must be of type funcArgs!')
    if not isinstance(args, pgramArgs):
        raise TypeError('Inappropriate argument type:\nargs must be of type pgramArgs!')
    time = fargs.getTime()
    mag = fargs.getMag()
    period = fargs.getPeriod()
    smooth = fargs.getSmoothedMag()
    ndata = len(time)
    nsamp = len(period)
    noutliers = args.nout

    # array to hold the deviation from the smoothed curve for each
    # data point
    tmpChi = fargs.getChi()

    # make sure we don't have more outliers than we have data points
    if noutliers > ndata:
        noutliers = ndata

    # determine reference deviations (recycle "tmpChi" array)
    meanMag = data.getMean(DATA_FIELD_TYPE.DATA_Y)
    for j in range(ndata):
        tmpChi[j] = (mag[j] - meanMag) ** 2
    tmpChi.sort()

    # sum the values most _poorly_ fit by the model mag=meanMag
    maxStd = 0.0
    maxChi = 0.0  # maxChi is the analogous var for each pd
    for j in range(ndata - 1, ndata - noutliers - 1, -1):
        maxStd += tmpChi[j]
    maxStd /= noutliers

    boxSize = fargs.boxSize
    myArgList = []
    startArr, endArr = getSplitNums(nsamp, args.nproc)
    for i in range(len(startArr)):
        myArgList.append([time, mag, makeCopy(smooth), boxSize, maxStd,
                          noutliers, period, startArr[i], endArr[i]])
    # Compute periodogram
    fargs.power = merge_arrs(pool.starmap(splitPlav, myArgList))


# A function to split up the tasks
def splitPlav(time, mag, mySmooth, boxSize, maxStd, noutliers, period,
              startNum, endNum):
    ret = []
    for i in range(startNum, endNum):
        ret.append(doPlav(time, mag, mySmooth, boxSize, maxStd,
                          noutliers, period[i]))
    return ret

# Actually does the algorithm computation.
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

# Function to fold/phase a light curve to an input period.
# Currently does smoothing based on args->smooth
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
        bLo = 0
        bHi = 0
        prevLo = 0
        prevHi = 0
        count = 0
        boxSum = 0.0
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
