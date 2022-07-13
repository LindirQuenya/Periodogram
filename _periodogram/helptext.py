#!/usr/bin/env python3
# This file contains the help texts, printed when the program is run
# with invalid arguments.

PGRAM_USAGE_TEXT = " Usage: python3 -m periodogram \n\t[-a <PeriodogramType (algorithm): one of ls, bls, plav>]" \
                   "\n\t[-b <NumberOfBins (-a bls only)>]" \
                   "\n\t[-c <NumWorkers>]" \
                   "\n\t[-d <FixedStepSize (-i fixedf or -i fixedp only)>]" \
                   "\n\t[-e <DataDelimiter>]" \
                   "\n\t[-f <FrequencyRangeMin> | -P <PeriodRangeMax>]" \
                   "\n\t[-F <FrequencyRangeMax> | -p <PeriodRangeMin>]" \
                   "\n\t[-H <DataHome directory: location of the input file>]" \
                   "\n\t[-i <PeriodStepMethod: one of std, exp, fixedf, fixedp, plav]" \
                   "\n\t[-K <StatNumberOfSamples>]" \
                   "\n\t[-L <OutFileLabeled>]" \
                   "\n\t[-M <StatMean> (not with -a ls)]" \
                   "\n\t[-n <NumberOfOutliers (-a plav only)>]" \
                   "\n\t[-N <NumberOfPeaksToReturn>]" \
                   "\n\t[-o <OversampleFactor (not with -i plav or -d)>]" \
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
#                  "\n\t[-D <output file name for debugging (may say 'stdout' or 'stderr')>]" \
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
                  "\n    Specifies which algorithm to run (one of ls, bls, plav). Defaults" \
                  "\n    to bls." \
                  "\n -b <NumberOfBins>" \
                  "\n    Specifies the number of bins to use in the bls algorithm" \
                  "\n -c <NumWorkers>" \
                  "\n    Specifies how many workers to use in the multiprocessing Pool." \
                  "\n    Defaults to the number of processers on the computer." \
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
                  "\n -H <DataHome>" \
                  "\n    The directory where the data is located. If not supplied," \
                  "\n    the location from which this is being run will be used." \
                  "\n -i <PeriodStepMethod>" \
                  "\n    Specifies which type of period stepping to use " \
                  "\n    (one of std, exp, fixedf, fixedp, plav)" \
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
                  "\n    table of top periods). The default is '.'" \
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
                  "\n    or -Y. Otherwise, supply all that you want used. If a labeled" \
                  "\n    column exists, but its label is not associated with a data type" \
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
                  "\n period and power, prints the internal argument state to stdout, " \
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

EPILOG_TEXT = "Results: " \
              "\n " \
              "\n If successful, periodogram creates an output table file containing " \
              "\n period and power, prints the internal argument state to stdout, " \
              "\n and exits with 0.  The output message contains the command line arguments" \
              "\n needed to replicate the exact results, including derived quantities, if any." \
              "\n " \
              "\n Examples: " \
              "\n " \
              "\n The following example runs periodogram on a table file with the period " \
              "\n range from .5 days to 1000 days and saves the output to outfile: " \
              "\n " \
              "\n $ python3 -m periodogram -p 0.5 -P 1000 test/test.csv test/outfile" \
              "\n "
