#!/usr/bin/env python3
# TODO: Make a comment for this.

from data_table import *
from preferences import *


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
        time = data.getFilteredArray(DATA_FIELD_TYPE.DATA_X)
        mag = data.getFilteredArray(DATA_FIELD_TYPE.DATA_Y)
        if len(time) != len(mag):
            raise RuntimeError("Houston, we have a problem. len(mag)!=len(time).")
        ndata = len(time)

        # Adjust the local copy of the data as appropirate for the different
        # algorithms
        meanMag = 0.0
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
                meanMag = data.getMean(DATA_FIELD_TYPE.DATA_Y)
            if adjustByMinDay:
                minDay = data.getMin(DATA_FIELD_TYPE.DATA_X)
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
    def populate(self, args, data):
        # set algo for populateLite (could just feed it args.algo,
        # but it's nice to have a local version.)
        algo = args.algo

        # populate and adjust time and mag
        self.populateLite(algo, data)
        ndata = self.ndata
        time = self.time
        mag = self.mag

        dell = TINY_NUM

        # We have set an extremely mild restriction: require at least 2
        # points in the file
        if ndata < MIN_NDATA:
            raise ValueError("InputFile: Not enough data in file to process: " + str(ndata))

        # Get minimum and maximum time values, along with the smallest
        # time-difference between two data points
        minDay = data.getMin(DATA_FIELD_TYPE.DATA_X)
        maxDay = data.getMax(DATA_FIELD_TYPE.DATA_X)
        minDt = data.getMinDiff(DATA_FIELD_TYPE.DATA_X)

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
            minperiod = data.getMedianDiff(DATA_FIELD_TYPE.DATA_X)

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
        self.boxSize = args.smooth
        self.ndata = ndata
        self.time = time
        self.mag = mag
        self.nsamp = nsamp
        self.period = period
        self.power = [0.0] * len(period)
        # Initialize chi for use in computePlavchan
        if algo == 'plav':
            self.chi = [0.0] * self.ndata

        # Initialize other arrays
        self.phase = [0.0] * self.ndata
        self.phasedMag = [0.0] * self.ndata
        self.smoothedMag = [0.0] * self.ndata

        # compute time estimate
        self.timeEst = estimateProcessingTime(nsamp, ndata, args, args.qmax, algo)
        if nsamp != len(self.period):
            raise RuntimeError("Error: nsamp!=len(period)")
        # print(period)

    def getTime(self):
        return self.time

    def getMag(self):
        return self.mag

    def getPeriod(self):
        return self.period

    def getPower(self):
        return self.power

    def getChi(self):
        return self.chi

    def getPhase(self):
        return self.phase

    def getPhasedMag(self):
        return self.phasedMag

    def getSmoothedMag(self):
        return self.smoothedMag

    def getSortable(self):
        return self.sortable

    # getPeakWidth()
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
    def getPeakWidth(self, useLog):
        sigLevel = 0

        # retrieve power array
        power = self.getPower()
        nsamp = len(power)

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
        self.width = width
        return width
