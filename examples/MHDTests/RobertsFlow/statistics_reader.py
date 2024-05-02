#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

take_last = int(10 / (5e-2))

the_statistics = np.transpose(np.loadtxt("statistics.txt"))

timestamps = the_statistics[1]
B = np.array(the_statistics[34])

dBdt = np.diff(B)/np.diff(timestamps)
cut_B = B[1:]
local_growth_rate = dBdt[-take_last:]/cut_B[-take_last:]
cut_timestamps = timestamps[-take_last:]

B = B / B[0]
divB = np.abs(np.array(the_statistics[35]))

growth_rate = round(np.mean(local_growth_rate),5)
growth_rate_err_2sigma = round(2*np.std(local_growth_rate),5)

name = "<s>="+str(growth_rate)+"$\pm$"+str(growth_rate_err_2sigma)

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
ax[0].plot(timestamps, B)
ax[1].plot(timestamps, divB)
ax[2].plot(cut_timestamps,local_growth_rate,label=name)
ax[0].set_xlabel("Time [s]")
ax[1].set_xlabel("Time [s]")
ax[2].set_xlabel("Time [s]")
ax[0].set_ylabel("$B$/ $B(t=0)$")
ax[0].set_yscale("log")
ax[1].set_ylabel("divB*h/B")
ax[1].set_yscale("log")
ax[2].legend(loc="best", fontsize=8)
ax[2].set_ylabel("d(LnB)/dt$")
plt.savefig("B_divB_vs_t.png", dpi=100)
