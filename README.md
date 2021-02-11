For all code that is run in this project, it is meant to be run from the command line with the command `python3 program.py args`, where `program.py` changes depending on which file is being run. 

# Periodogram
This is the python code for a periodogram (converted from C++).
There are three periodogram algorithms availiable: Box Least Squares, Lomb-Scargle, and Plavchan.
There are also five period-stepping types: standard, exponential, fixed period, fixed frequency, and plavchan.
To run without multiprocessing, use the commandline option "-c 1" to create only one thread.

## Multiprocessing
By default, the program will spawn as many threads as there are cores on the host machine.
Values for the time estimates are all based off of running the program on a dual-E5462 server.