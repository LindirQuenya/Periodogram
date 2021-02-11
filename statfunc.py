#!/usr/bin/env python3
# TODO: add a comment

from constants import *
from preferences import *

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
def statTrimOutliers(inputarr, threshold, maxfract):
    length = len(inputarr)
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
