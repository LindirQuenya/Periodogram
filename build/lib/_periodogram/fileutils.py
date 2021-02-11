#!/usr/bin/env python3
# TODO: make a comment for this

from _periodogram.constants import *


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
    if arrayNames[DATA_FIELD_TYPE.DATA_X] is None:
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
            if indMap[i] is not None:
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
