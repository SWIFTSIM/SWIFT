###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
#                    Loic Hausammann (loic.hausammann@id.ethz.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
import matplotlib

import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import h5py

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")


##########################################
# read specific energy from the solution
##########################################

# Read in units.
resultfile = "CoolingHeatingBox_results.txt"

f = open(resultfile, "r")
firstline = f.readline()
massline = f.readline()
lengthline = f.readline()
velline = f.readline()
f.close()
units = []
for l in [massline, lengthline, velline]:
    before, after = l.split("used:")
    val, unit = after.split("[")
    val = val.strip()
    units.append(float(val))

mass_units = units[0]
length_units = units[1]
velocity_units = units[2]
time_units = velocity_units / length_units
density_units = mass_units / length_units ** 3

# Read in all other data
data = np.loadtxt(resultfile)

Time = data[:, 1]
Time_Myr = Time * 1e-6
Energy = data[:, 12]  # internal energy IU


##########################################
# compute specific energy from the models
##########################################


# First snapshot
snap_filename = "snap/snapshot_0000.hdf5"

# Read the initial state of the gas
f = h5py.File(snap_filename, "r")

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]


# Read snapshots
files = glob("snap/snapshot_*.hdf5")
files.sort()
N = len(files)
energy_snap = np.zeros(N)
time_snap_cgs = np.zeros(N)
for i in range(N):
    snap = h5py.File(files[i], "r")
    masses = snap["/PartType0/Masses"][:]
    u = snap["/PartType0/InternalEnergies"][:] * masses
    u = sum(u) / masses.sum()
    energy_snap[i] = u
    time_snap_cgs[i] = snap["/Header"].attrs["Time"] * unit_time


##########################################
# plot
##########################################

plt.figure()

Myr_in_s = 1e6 * 365.25 * 24.0 * 60.0 * 60.0
yr_in_s = 365.25 * 24.0 * 60.0 * 60.0

plt.plot(Time_Myr, Energy, "r", ms=3, label="Expected solution")

plt.plot(time_snap_cgs / Myr_in_s, energy_snap, "k", ms=2)


plt.legend(loc="right", fontsize=8, frameon=False, handlelength=3, ncol=1)
plt.xlabel("${\\rm{Time~[Myr]}}$", labelpad=0)
plt.ylabel(r"$\rm{Internal Energy\,\,[IU]}$")


plt.tight_layout()

plt.show()
