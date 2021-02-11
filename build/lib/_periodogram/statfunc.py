#!/usr/bin/env python3
# TODO: add a comment

from _periodogram.preferences import *
from _periodogram.fileutils import *


# Returns the mean of arr, filtering unusable data points specified by omitRow:
# if omitRow[n] is nonzero, row n should be ignored.
def statMedian(arr, omitrow):
    if not isinstance(arr, list):
        raise TypeError("Inappropriate argument type:\narr must be a list!")
    if len(arr) <= 0:
        raise IndexError("Inappropriate index:\ncannot find the median of an empty list!")
    tmp = []
    count = 0
    for i in range(len(arr)):
        if omitrow == None or omitrow[i] == 0:
            tmp.append(arr[i])
            count += 1
    if count <= 0:
        raise ValueError("No usable points in array!")
    tmp.sort()
    if count % 2 == 1:
        md = tmp[int(int(count - 1) / 2)]
    else:
        md = (tmp[int(count / 2) - 1] + tmp[int(count / 2)]) / 2.0
    return md


# Returns the mean of arr. Filters bad data specified by omitRow
def statMean(arr, omitrow):
    if not isinstance(arr, list):
        raise TypeError("Inappropriate argument type:\narr must be a list!")
    if len(arr) <= 0:
        raise IndexError("Inappropriate index:\ncannot find the mean of an empty list!")
    if omitrow:
        arr1 = [arr[i] for i in range(len(arr)) if not omitrow[i]]
    else:
        arr1 = arr
    return sum(arr1) / len(arr1)


# Returns the standard deviation of an array
# Filters bad data specified by omitRow
def statStdDev(inarr, omitrow):
    if omitrow:
        arr = [inarr[i] for i in range(len(inarr)) if not omitrow[i]]
    else:
        arr = inarr
    s = sum(arr)
    sumsq = 0
    for i in arr:
        sumsq += i ** 2
    if len(arr) <= 0:
        raise ValueError("Cannot find standard deviation of empty array!")
    # "divide by n" method
    # return math.sqrt(len(arr)*sumsq-s**2)/len(arr)
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
def statTrimOutliers(length, inputarr, threshold, maxfract):
    nout = 100
    minCount = int(math.ceil((1.0 - maxfract) / length))

    x1 = [0.0] * length
    x2 = [0.0] * length
    its = 0
    # Set to error values, so that we know if the loop
    # never runs.
    mean = -1
    sd = -1
    for i in range(length):
        x1[i] = inputarr[i]
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


def findPeaks(args, fargs):

    # Make lists to hold significant periods and their powers
    sigPeriod = []
    sigPower = []

    isBls = 0
    isLs = 0
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
    period = fargs.getPeriod()
    power = fargs.getPower()
    nsamp = len(period)

    # determine the width of each peak
    width = fargs.getPeakWidth((isPlav or (USE_LOGNORMAL_BLS and isBls)))

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
        #        print(p)
        #        print(i)
        sortable[i] = [p, i]
    #    print(sortable)
    sortable.sort()
    #    print(sortable)
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
        myIdxJ = int(sortable[nsamp - j - 1][1])
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
                if 0 not in skipme:
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
