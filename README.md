For all code that is run in this project, it is meant to be run from the command line with the command `python3 program.py args`, where `program.py` changes depending on which file is being run. 

# Periodogram
This is the python code for a periodogram (converted from C++).
There are three periodogram algorithms availiable: Box Least Squares, Lomb-Scargle, and Plavchan.
There are also five period-stepping types: standard, exponential, fixed period, fixed frequency, and plavchan.
The file [`periodogram.py`](periodogram.py) does not use multiprocessing. It uses a loop (iterating through the list of periods) to find the power of each one. This is generally slower than multiprocessing.

## Multiprocessing
Multiprocessing programs are found in the [`mp`](mp) directory. There are two programs: [`mp_periodogram.py`](mp_periodogram.py) and [`mp_alt.py.`](mp_alt.py) `mp_periodogram` is generally faster. It splits up the task into a number of (close to) equal peices corresponding to the number of processors on the machine, and then adds each peice as a task that the pool of workers needs to complete. `mp_alt` adds each period seperately as a task to the pool. This requires much more overhead than `mp_periodogram,` and also requires making huge variables (lists that contain ndata*nsamp elements), so it slightly slower for small (~1000) numbers of periods, and much slower for large (~100000) numbers of periods. `mp_alt` was only left in for completeness. (Its use is very much not recommended!)

Time estimates are all based off of running the program on Peter's cluster. Only `periodogram` and `mp_periodogram` have usable time estimates. The time estimate for `mp_alt` is close to useless.
