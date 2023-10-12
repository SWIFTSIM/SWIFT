###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
import numpy as np

# Generates a swift IC file for the 3D Sod Shock in a periodic box

# Parameters
gamma = 2.0  # Gas adiabatic index
x_min = -1.0
x_max = 1.0
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state

Bx_L = 0.75
Bx_R = 0.75
By_L = 1.0
By_R = -1.0

fileName = "BrioWu.hdf5"

# ---------------------------------------------------
boxSize = x_max - x_min

glass_L = h5py.File("FCCglassCube_48.hdf5", "r")
glass_R = h5py.File("FCCglassCube_24.hdf5", "r")

times = 11
scale = boxSize / (2 * times)

pos_L = glass_L["/PartType0/Coordinates"][:, :] * scale
pos_R = glass_R["/PartType0/Coordinates"][:, :] * scale
h_L = glass_L["/PartType0/SmoothingLength"][:] * scale
h_R = glass_R["/PartType0/SmoothingLength"][:] * scale

numPart_glass_L = np.size(h_L)
numPart_glass_R = np.size(h_R)
numPart_L = numPart_glass_L * times
numPart_R = numPart_glass_R * times

pos_LL = np.zeros((numPart_L, 3))
pos_RR = np.zeros((numPart_R, 3))
h_LL = np.zeros(numPart_L)
h_RR = np.zeros(numPart_R)

for ii in range(times):
    bound_L = int(ii * numPart_glass_L)
    bound_R = int(ii * numPart_glass_R)
    pos_LL[bound_L : bound_L + numPart_glass_L, :] = pos_L + np.array(
        [ii * scale, 0.0, 0.0]
    )
    pos_RR[bound_R : bound_R + numPart_glass_R, :] = pos_R + np.array(
        [ii * scale, 0.0, 0.0]
    )
    h_LL[bound_L : bound_L + numPart_glass_L] = h_L
    h_RR[bound_R : bound_R + numPart_glass_R] = h_R

pos = np.append(pos_LL + np.array([x_min, 0.0, 0.0]), pos_RR, axis=0)
h = np.append(h_LL, h_RR)

numPart_L = np.size(h_LL)
numPart_R = np.size(h_RR)
numPart = np.size(h)

print(numPart)

vol_L = (boxSize / 2) * scale * scale
vol_R = (boxSize / 2) * scale * scale

# Generate extra arrays
v = np.zeros((numPart, 3))
B = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

B[:, 0] = Bx_L

for i in range(numPart):
    x = pos[i, 0]

    if x < 0:  # left
        u[i] = P_L / (rho_L * (gamma - 1.0))
        m[i] = rho_L * vol_L / numPart_L
        v[i, 0] = v_L
        B[i, 1] = By_L
    else:  # right
        u[i] = P_R / (rho_R * (gamma - 1.0))
        m[i] = rho_R * vol_R / numPart_R
        v[i, 0] = v_R
        B[i, 1] = By_R

# Shift particles
pos[:, 0] -= x_min

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, scale, scale]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
