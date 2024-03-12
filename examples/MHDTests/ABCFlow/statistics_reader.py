#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
E_mag = the_statistics[34]

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(10, 5))
ax.plot(Time, E_mag/E_mag[0])
ax.set_xlabel("Time [s]")
ax.set_ylabel("$E_{mag}$/ $E_{mag}(t=0)$")
ax.set_yscale('log')
plt.savefig("E_change.png", dpi=100)
