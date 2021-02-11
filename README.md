#UNSTABLE
This code is for developers only. Do not use it for anything you care about, as it can and will fail.

Currently, I'm hoping that this will work if you run `python -m periodogram` and then tack on some options.

# Periodogram
This is the python code for a periodogram (converted from C++).
There are three periodogram algorithms availiable: Box Least Squares, Lomb-Scargle, and Plavchan.
There are also five period-stepping types: standard, exponential, fixed period, fixed frequency, and plavchan.
To run without multiprocessing, use the commandline option "-c 1" to create only one thread.

## Multiprocessing
By default, the program will spawn as many threads as there are cores on the host machine.
Values for the time estimates are all based off of running the program on a dual-E5462 server.
