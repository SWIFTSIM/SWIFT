#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

take_last = int(10 / (5e-2))


def find_growth_rate(the_time, B_field, nlast=take_last):
    l = len(B_field)
    B_field_cut = B_field[l - 1 - nlast : -1]
    time_cut = the_time[l - 1 - nlast : -1]
    # print(time_cut,B_field_cut)
    res = np.polyfit(time_cut, B_field_cut, 1)[0]
    return np.round(res, 4)


the_statistics = np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
B = np.array(the_statistics[34])
B = B / B[0]
divB = np.abs(np.array(the_statistics[35]))

print("The growth rate is: " + str(find_growth_rate(Time, np.log(B))))

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(10, 5))
ax[0].plot(Time, B)
ax[1].plot(Time, divB)
ax[0].set_xlabel("Time [s]")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("$B$/ $B(t=0)$")
ax[0].set_yscale("log")
ax[1].set_ylabel("divB*h/B")
ax[1].set_yscale("log")
plt.savefig("B_divB_vs_t.png", dpi=100)
