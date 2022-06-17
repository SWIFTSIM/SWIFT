###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a swift IC file for the 1D Toro (2009) test 2 in a periodic box
# Fun fact: because of the periodic boundaries, we get a second test for free

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
numPart = 1000  # Number of particles
x_min = -1.0
x_max = 1.0
rho_L = 1.0  # Density left state
rho_R = 1.0  # Density right state
v_L = -2.0  # Velocity left state
v_R = 2.0  # Velocity right state
P_L = 0.4  # Pressure left state
P_R = 0.4  # Pressure right state
fileName = "toroTest2.hdf5"

unit_l_in_cgs = 3.086e18
unit_m_in_cgs = 2.94e55
unit_t_in_cgs = 3.086e18
a_beg = 0.001

# ---------------------------------------------------

boxSize = x_max - x_min
delta = boxSize / numPart

# Build the arrays
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
h = zeros(numPart)
u = zeros(numPart)

# Set the particles on the left
for i in range(numPart // 2):
    coords[i, 0] = x_min + (i + 0.5) * delta
    u[i] = P_L / (rho_L * (gamma - 1.0))
    h[i] = 1.2348 * delta
    m[i] = delta * rho_L
    v[i, 0] = v_L

# Set the particles on the right
for j in range(numPart // 2):
    i = numPart // 2 + j
    coords[i, 0] = x_min + (i + 0.5) * delta
    u[i] = P_R / (rho_R * (gamma - 1.0))
    h[i] = 1.2348 * delta
    m[i] = delta * rho_R
    v[i, 0] = v_R

u /= a_beg ** (3.0 * (gamma - 1.0))
v /= a_beg

# Shift particles
coords[:, 0] -= x_min

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 1

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_l_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_m_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_t_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=coords, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
