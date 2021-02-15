#!/usr/bin/env python3
# This file runs the overall program. It is meant to be run from the
# command-line, with the command "python3 -m periodogram ....". However,
# the module can also be imported and the functions can be used in other
# Python programs.

from _periodogram.driver import *


# This runs the program.
if __name__ == '__main__':
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
    data.populate(args.intbl, args.xcol, args.ycol, args.yerrCol,
                  args.constraintCol, args.constraintMin, args.constraintMax, args.inputDelimiter)

    # Construct and populate the funcArgs object
    fargs = funcArgs()
    fargs.populate(args, data)
    with mp.Pool(processes=args.nproc) as pool:
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
    print(args.argsPrint("-m periodogram"))

    # "stop" the timer
    t_f = time.time()
    t_d = t_f - t_i
    if PRINT_TIMES:
        # Print the result from our timer
        print("Time elapsed: " + format(t_d, '.4f') + " seconds (" + format(t_d / 60, '.4f') + " minutes)")