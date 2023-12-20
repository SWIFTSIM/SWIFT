#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
B = the_statistics[38]
Ekin = the_statistics[13]
Eint = the_statistics[14]
Emag = the_statistics[34] 
Etot = Ekin+Eint+Emag
Beq = np.sqrt(2*Ekin)
divB = np.abs(the_statistics[35])

fig, ax = plt.subplots(1, 3, sharex=False, figsize=(15, 5))
ax[0].plot(Time, B/B[0], label = 'B/$B_0$', color = 'red')
ax[0].plot(Time, Beq/B[0], label = '$B_{eq}$/$B_0$', color = 'black')
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("B / $B_0$")
ax[0].set_yscale("log")
ax[0].legend()
#ax[0].set_xscale("log")
ax[0].grid()

ax[1].plot(Time, Ekin, label = '$E_{kin}$',color = 'blue')
ax[1].plot(Time, Eint, label = '$E_{int}$',color = 'green')
ax[1].plot(Time, Emag, label = '$E_{mag}$',color = 'red')
#ax[1].plot(Time, Etot, label = '$E_{tot}$',color = 'black')

ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel("$E$")
ax[1].set_yscale("log")
#ax[1].set_ylim([1,1e2])
ax[1].grid()
ax[1].legend()
#ax[1].set_xscale("log")

ax[2].plot(Time, divB, label = '$\langle$ divB * h/B $\rangle$', color = 'red')
ax[2].set_yscale("log")

plt.tight_layout()
plt.savefig("B_change.png", dpi=100)
