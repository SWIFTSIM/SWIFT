#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
B = np.array(the_statistics[34])
B = B/B[0]
divB = np.abs(np.array(the_statistics[35]))

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(10, 5))
ax[0].plot(Time, B)
ax[1].plot(Time, divB)
ax[0].set_xlabel("Time [s]")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("$B$/ $B(t=0)$")
ax[0].set_yscale('log')
ax[1].set_ylabel("divB*h/B")
ax[1].set_yscale('log')
plt.savefig("B_divB_vs_t.png", dpi=100)
