#!/usr/bin/env python3
# TODO:add comment

from func_args import *
from cli_args import *

# computeBLS()
# Function to compute the BLS "periodogram"

# BLS = Box-fitting Least Squares
# ref: Kovacs, G., Zucker, S. and Mazeh, T. "A box-fitting algorithm
#     in the search for periodic transits." A&A 391:369-377 (2002).
# http://adsabs.harvard.edu/abs/2002A%26A...391..369K
# The BLS algorithm starts from the premise that for a specific fraction
# of the period of an orbiting planet, the planet will transit in front
# of its star.  This time during which the star's light is
# obstructed ranges from qmin to qmax, expressed as a fraction of the
# total period.

# For each candidate period p, the number of bins (nbins) is considered
# to span one period: each bin corresponds to a time span of p/nbins.

# The observed data is "folded" to match the period: observations
# at time t = p + dt are placed into the bin corresponding to dt.

# A model in which the mean signal level in the occluded phase is L and
# the level in the un-occluded phase is H is considered for each
# candidate length of the L phase (qmin * nbins to qmax * nbins).  The
# least squares fit is given by maximizing s**2/(r*(1-r)) where
# s is the weighted sum of magnitudes in the low period and r
# the sum of the weights in the low period.

# Arguments:
#  data = populated dataTbl object
#  args = populated pgramArgs object
#  fargs = funcArgs object with ndata, time, mag
#               nsamp, and period set. power will
#               be populated by this function
#  pool = the pool object to give tasks to
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
    time = fargs.getTime()
    mag = fargs.getMag()
    period = fargs.getPeriod()
    ndata = len(time)
    nsamp = len(period)

    # Initializes variables with nonsense values
    blsR = [0.0] * nsamp
    blsS = [0.0] * nsamp
    lowBin0 = [0] * nsamp
    lowBin1 = [0] * nsamp

    wt, err = [], []
    # In case we want weight as a function of uncertainty
    if WEIGHT_BY_ERR:
        err = data.getFilteredArray(DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY)
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
    myArgList = []
    startArr, endArr = getSplitNums(nsamp, args.nproc)
    for i in range(len(startArr)):
        myArgList.append([time, mag, wt, period, nbins, binExt, minBins, minWt,
                          totalWt, binMax, startArr[i], endArr[i]])

    # Compute periodogram
    results = pool.starmap(splitBLS, myArgList)
    ##    results=pool.starmap(doBLS,transpose([[time]*nsamp,[mag]*nsamp,[wt]*nsamp,\
    ##                                      makeArrOfCopies(binWt,nsamp),\
    ##                                      makeArrOfCopies(binMag,nsamp),\
    ##                                      period,[nbins]*nsamp,[nsamp]*nsamp,\
    ##                                          [binExt]*nsamp,[minBins]*nsamp,\
    ##                                          [minWt]*nsamp,[totalWt]*nsamp]))
    results = transpose(merge_arrs(results))
    fargs.power = results[0]
    fargs.blsR = results[1]
    fargs.blsS = results[2]

    # Commented out because they are unused
    # fargs.lowBin0=results[3]
    # fargs.lowBin1=results[4]


# stopNum is exclusive
def splitBLS(time, mag, wt, period, nbins, binExt, minBins, minWt, totalWt, binMax,
             startNum, stopNum):
    nsamp = len(period)
    res = []
    for i in range(startNum, stopNum):
        binWt = [0.0] * binMax
        binMag = [0.0] * binMax
        res.append(doBLS(time, mag, wt, binWt, binMag, period[i], nbins, nsamp, binExt,
                         minBins, minWt, totalWt))
    return res


# time and mag are arrs. binWt and binMag are working vars,
# so they can't be shared. copies will have to be made for
# each run of the function.
def doBLS(time, mag, wt, binWt, binMag, period, nbins, nsamp, binExt, minBins, minWt,
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
            if binCt >= minBins and minWt <= sumWt < totalWt:
                pwr = (sumMag ** 2) / (sumWt * (totalWt - sumWt))
                if pwr >= maxPwr:
                    maxPwr = pwr
                    lowStart = b
                    lowEnd = k
                    lowWt = sumWt
                    lowMag = sumMag
    maxPwr = math.sqrt(maxPwr)
    if maxPwr > 0:
        return maxPwr, lowWt / totalWt, lowMag  # ,lowStart,lowEnd)
    return 0, 0, 0  # ,0,0)