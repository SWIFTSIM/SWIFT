################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import h5py
from numpy import *

# Generates a swift IC file for the Gresho-Chan vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1  # Gas density
P0 = 0.0  # Constant additional pressure (should have no impact on the dynamics)
fileOutputName = "greshoVortex.hdf5"
fileGlass = "glassCube_64.hdf5"
# ---------------------------------------------------

# Get position and smoothing lengths from the glass
fileInput = h5py.File(fileGlass, "r")
coords = fileInput["/PartType0/Coordinates"][:, :]
h = fileInput["/PartType0/SmoothingLength"][:]
ids = fileInput["/PartType0/ParticleIDs"][:]
boxSize = fileInput["/Header"].attrs["BoxSize"][0]
numPart = size(h)
fileInput.close()

# Now generate the rest
m = ones(numPart) * rho0 * boxSize ** 3 / numPart
u = zeros(numPart)
v = zeros((numPart, 3))
mat = zeros(numPart)

for i in range(numPart):

    x = coords[i, 0]
    y = coords[i, 1]

    r2 = (x - boxSize / 2) ** 2 + (y - boxSize / 2) ** 2
    r = sqrt(r2)

    v_phi = 0.0
    if r < 0.2:
        v_phi = 5.0 * r
    elif r < 0.4:
        v_phi = 2.0 - 5.0 * r
    else:
        v_phi = 0.0
    v[i, 0] = -v_phi * (y - boxSize / 2) / r
    v[i, 1] = v_phi * (x - boxSize / 2) / r
    v[i, 2] = 0.0

    P = P0
    if r < 0.2:
        P = P + 5.0 + 12.5 * r2
    elif r < 0.4:
        P = P + 9.0 + 12.5 * r2 - 20.0 * r + 4.0 * log(r / 0.2)
    else:
        P = P + 3.0 + 4.0 * log(2.0)
    u[i] = P / ((gamma - 1.0) * rho0)


# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = fileOutput.create_group("/PartType0")
ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
ds[()] = coords
ds = grp.create_dataset("Velocities", (numPart, 3), "f")
ds[()] = v
ds = grp.create_dataset("Masses", (numPart, 1), "f")
ds[()] = m.reshape((numPart, 1))
ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
ds[()] = h.reshape((numPart, 1))
ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
ds[()] = u.reshape((numPart, 1))
ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
ds[()] = ids.reshape((numPart, 1))
ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
ds[()] = mat.reshape((numPart, 1))

fileOutput.close()
