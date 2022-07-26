###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import h5py

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

# Some constants in cgs units
k_b_cgs = 1.38e-16  # boltzmann
m_h_cgs = 1.67e-24  # proton mass

# File containing the total energy
stats_filename = "./statistics.txt"

# First snapshot
snap_filename = "coolingBox_0000.hdf5"

# Read the initial state of the gas
f = h5py.File(snap_filename, "r")

# Read the units parameters from the snapshot
units = f["InternalCodeUnits"]
unit_mass = units.attrs["Unit mass in cgs (U_M)"]
unit_length = units.attrs["Unit length in cgs (U_L)"]
unit_time = units.attrs["Unit time in cgs (U_t)"]

# Read the adiabatic index
gamma = float(f["HydroScheme"].attrs["Adiabatic index"])


def energyUnits(u):
    """ Compute the temperature from the internal energy. """
    u *= (unit_length / unit_time) ** 2
    return u * m_h_cgs / k_b_cgs


# Read energy and time arrays
array = np.genfromtxt(stats_filename, skip_header=1)
time = array[:, 1] * unit_time
total_mass = array[:, 4]
kinetic_energy = array[:, 13]
internal_energy = array[:, 14]
potential_energy = array[:, 15]
radiated_energy = array[:, 16]
total_energy = kinetic_energy + internal_energy + potential_energy
initial_energy = total_energy[0]

# Conversions to cgs
total_energy_cgs = total_energy / total_mass[0]
total_energy_cgs = energyUnits(total_energy_cgs)

kinetic_energy_cgs = kinetic_energy / total_mass[0]
kinetic_energy_cgs = energyUnits(kinetic_energy_cgs)

internal_energy_cgs = internal_energy / total_mass[0]
internal_energy_cgs = energyUnits(internal_energy_cgs)

radiated_energy_cgs = radiated_energy / total_mass[0]
radiated_energy_cgs = energyUnits(radiated_energy_cgs)

# Read snapshots
files = glob("coolingBox_*.hdf5")
N = len(files)
temp_snap = np.zeros(N)
time_snap_cgs = np.zeros(N)
for i in range(N):
    snap = h5py.File(files[i], "r")
    u = snap["/PartType0/InternalEnergies"][:] * snap["/PartType0/Masses"][:]
    u = sum(u) / total_mass[0]
    temp_snap[i] = energyUnits(u)
    time_snap_cgs[i] = snap["/Header"].attrs["Time"] * unit_time


plt.figure()

Myr_in_s = 1e6 * 365.25 * 24.0 * 60.0 * 60.0
plt.plot(time / Myr_in_s, total_energy_cgs, "r-", lw=1.6, label="Gas total energy")
# statistics and snapshots may not be at same timestep and frequency
plt.plot(time_snap_cgs / Myr_in_s, temp_snap, "rD", ms=3)
plt.plot(time / Myr_in_s, radiated_energy_cgs, "g-", lw=1.6, label="Radiated energy")
plt.plot(
    time / Myr_in_s,
    total_energy_cgs + radiated_energy_cgs,
    "b-",
    lw=0.6,
    label="Gas total + radiated",
)

plt.legend(loc="right", fontsize=8, frameon=False, handlelength=3, ncol=1)
plt.xlabel("${\\rm{Time~[Myr]}}$", labelpad=0)
plt.ylabel("${\\rm{Internal ~Energy ~(u ~m_H / k_B) ~[K]}}$")

plt.tight_layout()

plt.savefig("energy.png", dpi=200)
