#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

filename = os.path.abspath(sys.argv[1])

with open(filename, 'r') as f:
    k = f.readlines()

period=[]
power=[]
for i in range(len(k)):
    splitrow = k[i].split(",")
    period.append(float(splitrow[0]))
    power.append(float(splitrow[1]))

plt.semilogx(period,power)
plt.xlabel("Period")
plt.ylabel("Power")
plt.show()
