###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import pylab as pl
import sys

if len(sys.argv) < 2:
    print("No final snapshot number provided!")
    exit()
lastsnap = int(sys.argv[1])

# Read the simulation data
t = np.zeros(lastsnap + 1)
ey = np.zeros(lastsnap + 1)
b2 = np.zeros(lastsnap + 1)
ei = np.zeros(lastsnap + 1)
ek = np.zeros(lastsnap + 1)

for i in range(lastsnap + 1):
    file = h5py.File("KelvinHelmholtzMHD_{0:04d}.hdf5".format(i), "r")
    t_snap = float(file["/Header"].attrs["Time"])
    vy = file["/PartType0/Velocities"][:, 1]
    vx = file["/PartType0/Velocities"][:, 0]
    vz = file["/PartType0/Velocities"][:, 2]
    by = file["/PartType0/MagneticFluxDensities"][:, 1]
    bx = file["/PartType0/MagneticFluxDensities"][:, 0]
    bz = file["/PartType0/MagneticFluxDensities"][:, 2]
    m = file["/PartType0/Masses"][:]
    ee = file["/PartType0/InternalEnergies"][:]

    ey_snap = 0.5 * m * vy ** 2
    ek_snap = 0.5 * m * (vy ** 2 + vx ** 2 + vz ** 2)
    b2_snap = 0.5 * (by ** 2 + bx ** 2 + bz ** 2)

    t[i] = t_snap
    ey[i] = ey_snap.sum()
    b2[i] = b2_snap.sum()
    ek[i] = ek_snap.sum()
    ei[i] = ee.sum()


fi, ax = pl.subplots(2, 2)
# pl.semilogy(t, ey, "k.")
ax[0, 0].semilogy(t, ey, "k.")
ax[0, 0].set_title("Kinetik Y")
ax[1, 0].semilogy(t, b2, "k.")
ax[1, 0].set_title("b2")
ax[0, 1].semilogy(t, ek, "k.")
ax[0, 1].set_title("Kinnetik")
ax[1, 1].semilogy(t, ei, "k.")
ax[1, 1].set_title("Internal")
fi.tight_layout()

fi.savefig("kelvinHelmholtzMHD.png")
