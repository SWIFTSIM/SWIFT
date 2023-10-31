###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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

MyrInSec = 31557600000000.0
gcm3InAcc = 1 / 5.978637406556783e23


def ComputeDensity(snap):

    # Read the initial state of the gas
    f = h5py.File(snap, "r")

    # Read the units parameters from the snapshot
    units = f["InternalCodeUnits"]
    unit_mass = units.attrs["Unit mass in cgs (U_M)"]
    unit_length = units.attrs["Unit length in cgs (U_L)"]
    unit_time = units.attrs["Unit time in cgs (U_t)"]
    unit_density = unit_mass / unit_length ** 3

    # Header
    header = f["Header"]
    BoxSize = header.attrs["BoxSize"]
    Time = header.attrs["Time"]

    # Read data
    pos = f["/PartType0/Coordinates"][:]
    rho = f["/PartType0/Densities"][:]
    mass = f["/PartType0/Masses"][:]
    ids = f["/PartType0/ParticleIDs"][:]

    # Center the model and compute particle radius
    pos = pos - BoxSize / 2
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    r = np.sqrt(x * x + y * y + z * z)
    n = len(r)

    # sort particles according to their distance
    idx = np.argsort(r)
    r = r[idx]
    mass = mass[idx]
    ids = ids[idx]

    nparts_per_bin = 100

    # loop over bins containing nparts_per_bin particles
    idx = np.arange(n)

    # bins radius and mass
    rs_beg = np.array([])
    rs_end = np.array([])
    ms = np.array([])

    i = 0

    while i + nparts_per_bin < n:

        rs_beg = np.concatenate((rs_beg, [r[i]]))
        rs_end = np.concatenate((rs_end, [r[i + nparts_per_bin]]))
        m_this_bin = np.sum(mass[i : i + nparts_per_bin])
        ms = np.concatenate((ms, [np.sum(m_this_bin)]))

        # shift
        i = i + nparts_per_bin

    # compute density
    vol = 4 / 3 * np.pi * (rs_end ** 3 - rs_beg ** 3)
    rho = ms / vol

    # compute radius, we use the mean
    rs = 0.5 * (rs_beg + rs_end)

    # convert rho to acc
    rho = rho * unit_density / gcm3InAcc

    # convert time to Myr
    Time = Time * unit_time / MyrInSec

    return rs, rho, Time


# Do the plot

plt.figure()

rs, rho, Time = ComputeDensity("snapshot_0000.hdf5")
plt.plot(rs, rho, c="b", label=r"$t=%5.1f\,\rm{[Myr]}$" % Time, lw=1)

rs, rho, Time = ComputeDensity("snapshot_0050.hdf5")
plt.plot(rs, rho, c="r", label=r"$t=%5.1f\,\rm{[Myr]}$" % Time, lw=1)

plt.loglog()

plt.legend()
plt.xlabel("${\\rm{Radius~[kpc]}}$", labelpad=0)
plt.ylabel("${\\rm{Density~[atom/cm^3]}}$", labelpad=0)


plt.savefig("GasDensityProfile.png", dpi=200)
