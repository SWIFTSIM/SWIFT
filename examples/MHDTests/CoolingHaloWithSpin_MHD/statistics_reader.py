#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
from glob import glob
from swiftsimio import load_statistics
import unyt as u

filename = sys.argv[1]
outputname = sys.argv[2]

data = load_statistics(filename)

time = data.time

Ekin = data.kin_energy.to(u.erg)
Eint = data.int_energy.to(u.erg)
Emag = data.mag_energy.to(u.erg)
Epot = np.abs(data.pot_energy.to(u.erg))

Etot = np.abs(Ekin + Eint + Emag - Epot)

name = ""

# print(max(Epot))

fig, ax = plt.subplots(1, 1, sharex=True, figsize=(7, 5))
ax.plot(time, Emag, label="Emag" + name, color="green", linestyle="dashed")
ax.plot(time, Ekin, label="Ekin" + name, color="blue", linestyle="dashed")
ax.plot(time, Eint, label="Eint" + name, color="red", linestyle="dashed")
ax.plot(time, Epot, label="|Epot|" + name, color="brown", linestyle="dotted")
# ax.plot(time, Epot,label='Etot'+name,color='black', linestyle='dotted')
ax.set_xlabel("time [Gyr]")
ax.set_ylabel(f"Energy [erg]")
ax.set_ylim([1e48, 1e60])
ax.set_yscale("log")
ax.grid()
ax.legend()
plt.savefig(outputname, dpi=220)
