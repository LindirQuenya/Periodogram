#!/usr/bin/env python3
# TODO: Make a comment for this.

from statfunc import *
from utilfunc import *


# Class used to hold data from files, column names (labels),
# filtering info, and various statistics about each set of
# data. It is organized according to DATA_FIELD_TYPE, an enum
# that dictates which indicies denote which type of data (xdata,
# ydata,yUncertainty, etc)
class dataTbl:
    def __init__(self):
        self.description = None  # UNUSED
        self.origin = None  # UNUSED

        self.ndata = 0  # number of data points
        self.ndataUnfilt = 0  # number of data points not filtered
        self.dataArray = []  # each index contains a dataArray of
        # length ndata corresponding to the
        # DATA_FIELD_TYPE of that index
        self.colKeyNames = []  # UNUSED
        self.arrayNames = []  # contains the column names associated
        # with the DATA_FIELD_TYPE of that index
        self.units = []  # UNUSED
        self.median = []  # the median for each array
        self.mean = []  # the mean for each array
        self.sd = []  # the std dev for each array
        for i in range(DATA_FIELD_TYPE.DATA_N_TYPES):
            self.colKeyNames.append(None)
            self.arrayNames.append(None)
            self.dataArray.append(None)
            self.units.append(None)
            self.median.append(UNSET_MEAN)
            self.mean.append(UNSET_MEAN)
            self.sd.append(UNSET_MEAN)

        self.omitRow = None  # array containing flags to skip a row.
        # length:ndata
        self.tinfo = None  # UNUSED
        self.nhead = 0  # UNUSED
        self.header = 0  # UNUSED

    # Function to populate the dataTbl. Calls readDoubleColumns to
    # get data from file
    def populate(self, intbl, xcol, ycol, yerrcol, constraintcol,
                 constraintmin, constraintmax, delim):
        # set the column names
        self.setColName(DATA_FIELD_TYPE.DATA_X, xcol)
        self.setColName(DATA_FIELD_TYPE.DATA_Y, ycol)
        self.setColName(DATA_FIELD_TYPE.DATA_CONSTRAINT, constraintcol)
        self.setColName(DATA_FIELD_TYPE.DATA_Y_UNCERTAINTY, yerrcol)

        # Even though we don't use it, let's set it anyway.
        self.origin = intbl

        # Index by which we sort data will be the xdata
        sortIdx = DATA_FIELD_TYPE.DATA_X

        # Read the data, omitting rows that contain non-numerical data
        unsortedData = readDoubleColumns(intbl, self.arrayNames, delim)
        # Sort the data
        self.dataArray = sortColumns(unsortedData, sortIdx)
        self.ndata = len(self.dataArray[DATA_FIELD_TYPE.DATA_X])
        if self.ndata <= 0:
            raise ValueError("InputFile: No data found in file!")

        # Set an initial value for omitRow
        self.omitRow = [None] * self.ndata

        # print(this.omitRow==None)
        # Filter bad data
        self.setFilter(constraintmin, constraintmax)

    # Function to set the arrayNames entry to the specified value
    def setColName(self, f, s):
        if f >= DATA_FIELD_TYPE.DATA_N_TYPES:
            raise ValueError('Inappropriate value:\nf is greater than maximum value!')
        if s:
            self.arrayNames[f] = s

    # Function to set omitRow flags based on values in the
    # constraint (aka filter) column
    def setFilter(self, minval, maxval):
        # print(this.omitRow==None)
        # Initially, nothing filtered
        self.ndataUnfilt = self.ndata

        # check that filtering is applicable here
        if self.arrayNames[DATA_FIELD_TYPE.DATA_CONSTRAINT] is not None and \
                self.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT] is not None:

            # have we set lower and upper limits?
            if minval == UNSET_VALUE or maxval == UNSET_VALUE:
                tmp = self.arrayNames[DATA_FIELD_TYPE.DATA_CONSTRAINT].lower()
                # status!=0 signifies error - do not allow anything else
                if "status" in tmp:
                    minval = 0
                    maxval = 0
                else:
                    # identify the minimum and maximum values -- not constrained
                    # but informational
                    myArray = self.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT]
                    minval = min(myArray)
                    maxval = max(myArray)
                    return None  # no restrictions set -> nothing to do
            # if we get here, we have a constraint field and min and max values
            nfilt = 0
            for i in range(self.ndata):
                # if anything violates the max/min values
                if not self.omitRow[i] and \
                        (self.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT][i] < minval or
                         self.dataArray[DATA_FIELD_TYPE.DATA_CONSTRAINT][i] > maxval):
                    # exclude it
                    self.omitRow[i] = 1
                    nfilt += 1
            if nfilt == self.ndata:
                raise ValueError("InputFile: all data filtered by constraints!")
            self.ndataUnfilt = self.ndata - nfilt

    # Returns a value signifying whether the mean for the data type
    # specified by dataType is set. -1:error,0:unset,1:set
    def isMeanSet(self, datatype):
        try:
            if self.mean[datatype] == UNSET_MEAN:
                return 0
            return 1
        except IndexError:
            return -1
        except TypeError:
            return -1

    # Returns a value signifying whether the standard deviation for
    # the data type specified by dataType is set. -1:error,0:unset,1:set
    def isDevSet(self, datatype):
        try:
            if self.sd[datatype] == UNSET_MEAN or self.sd[datatype] == 0:
                return 0
            return 1
        except IndexError:
            return -1
        except TypeError:
            return -1

    # Sets the mean for datatype to mean.
    def setMean(self, datatype, mean):
        try:
            self.mean[datatype] = mean
        except IndexError:
            raise ValueError('Inappropriate value:\ndatatype is greater than maximum value!')

    # Sets the standard deviation for datatype to deviation.
    def setDev(self, datatype, deviation):
        try:
            self.sd[datatype] = deviation
        except IndexError:
            raise ValueError('Inappropriate value:\ndatatype is greater than maximum value!')

    # Returns the data array for datatype.
    def getArray(self, datatype):
        try:
            return self.dataArray[datatype]
        except IndexError:
            raise ValueError('Inappropriate value:\ndatatype is greater than maximum value!')

    # Returns the mean for datatype.
    # Calculates and returns the mean if it's unset.
    def getMean(self, datatype):
        try:
            if self.isMeanSet(datatype) < 1:
                myArray = self.getArray(datatype)
                m = statMean(myArray, self.omitRow)
                self.setMean(datatype, m)
            return self.mean[datatype]
        except IndexError:
            raise ValueError('Inappropriate value:\ndatatype is greater than maximum value!')

    # Returns the standard deviation for datatype.
    # Calculates and returns the stddev if it's unset.
    def getDev(self, datatype):
        try:
            if self.isMeanSet(datatype) < 1:
                myArray = self.getArray(datatype)
                d = statStdDev(myArray, self.omitRow)
                self.setDev(datatype, d)
            return self.mean[datatype]
        except IndexError:
            raise ValueError('Inappropriate value:\ndatatype is greater than maximum value!')

    # Returns the data array for datatype, filtered by omitRow:
    # omitRow is 0 or None for good rows, and any other value
    # for rows that should be ignored.
    def getFilteredArray(self, datatype):
        # No need to run checks, getArray should do that.
        refArray = self.getArray(datatype)
        myArray = []
        for i in range(self.ndata):
            if self.omitRow is None or not self.omitRow[i]:
                myArray.append(refArray[i])
        return myArray

    # Returns the minimum value of the filtered array for datatype
    def getMin(self, datatype):
        try:
            arr = self.getFilteredArray(datatype)
            return min(arr)
        except ValueError:
            raise ValueError("Inappropriate value:\nt must be populated before getMin can be called!")

    # Returns the maximum value of the filtered array for datatype
    def getMax(self, datatype):
        try:
            arr = self.getFilteredArray(datatype)
            return max(arr)
        except ValueError:
            raise ValueError("Inappropriate value:\nt must be populated before getMax can be called!")

    # Returns the minimum difference between two adjacent elements
    # in the array for datatype
    def getMinDiff(self, datatype):
        # getArray will run checks
        myArray = self.getArray(datatype)
        minDiff = myArray[1] - myArray[0]
        if minDiff == 0:
            raise ValueError('InputFile: two measurements cannot exist at the same time!')
        diff = 0
        for i in range(1, self.ndata):
            diff = myArray[i] - myArray[i - 1]
            if diff < 0:
                raise ValueError("Array must be sorted before calling dtGetMinDiff!")
            if 0 < diff < minDiff:
                minDiff = diff
        return minDiff

    # Returns the median difference between two adjacent elements
    # in the array t.dataArray[f].
    def getMedianDiff(self, datatype):
        # getArray will run checks
        myArray = self.getArray(datatype)
        diff = []
        for i in range(1, self.ndata):
            diff.append(myArray[i] - myArray[i - 1])
            if diff[i - 1] < 0:
                raise ValueError("Array must be sorted before calling dtGetMinDiff!")
        return statMedian(diff, None)
