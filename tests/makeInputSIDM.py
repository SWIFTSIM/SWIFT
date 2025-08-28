###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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

import h5py
from numpy import *

# Generates a swift IC file containing a cartesian distribution of SIDM particles
# in a cubic box

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1.0
L = 4  # Number of particles along one axis
fileName = "input_SIDM.hdf5"

# ---------------------------------------------------
numPart = L ** 3
mass = 1

# Generate particles
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
m = zeros((numPart, 1))
ids = zeros((numPart, 1), dtype="L")

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i * L * L + j * L + k
            x = i * boxSize / L + boxSize / (2 * L)
            y = j * boxSize / L + boxSize / (2 * L)
            z = k * boxSize / L + boxSize / (2 * L)
            coords[index, 0] = x
            coords[index, 1] = y
            coords[index, 2] = z
            v[index, 0] = 0.0
            v[index, 1] = 0.0
            v[index, 2] = 0.0
            m[index] = mass
            ids[index] = index


# --------------------------------------------------

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [0, 0, 0, 0, 0, 0, 0, numPart]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, 0, 0, 0, 0, 0, 0, numPart]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0, 0, 0]


# Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0


# Particle group
grp = file.create_group("/PartType7")
ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
ds[()] = coords
ds = grp.create_dataset("Velocities", (numPart, 3), "f")
ds[()] = v
ds = grp.create_dataset("Masses", (numPart, 1), "f")
ds[()] = m
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids

file.close()
