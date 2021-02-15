# UNSTABLE
This is the dev branch. It is only for people who want to contribute or bug-test. Do not use it for anything you care about, as it can and will fail.

If you find any problems/bugs/etc in it, please open an [issue](https://github.com/LindirQuenya/Periodogram/issues). I very much appreciate such contributions! See the [bugs](#bugs) section for more info.

# Periodogram
This is the python code for a periodogram (converted from C++).
It is currently availiable on [testing PyPI](https://test.pypi.org/project/periodogram), and I have plans to migrate it to normal PyPI sometime soon.
There are three periodogram algorithms availiable: Box Least Squares, Lomb-Scargle, and Plavchan.
There are also five period-stepping types: standard, exponential, fixed period, fixed frequency, and plavchan.

## Using the code
Given the instability of the code, I would recommend that you [create a dedicated venv](https://docs.python.org/3/library/venv.html) for the periodogram.
You should be able to install the periodogram with `pip install -i https://test.pypi.org/simple/ periodogram`.
With any luck, it will "just work" if you run `python -m periodogram [options] InputFile [OutputFile]`.
If you want to see the options, run `python -m periodogram --help`.

## Multiprocessing
By default, the program will spawn as many threads as there are cores on the host machine.
Values for the time estimates are all based off of running the program on a dual-E5462 server.
To run without multiprocessing, use the command-line option "-c 1" to create only one thread.

## Bugs
When filing a bug report, please remember to include:
 - The error message.
 - The command-line arguments that triggered this bug.
 - The output of `pip freeze`.
 - The dataset that you're running the periodogram on.
 - The version of the code that you're using.
 - Anything else that might be helpful to replicate the bug.
Thank you!
