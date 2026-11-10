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

import unyt
import swiftsimio as sw
from numpy import *

# Generates a swift IC file for the 3D Sod Shock in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
x_min = -1.0
x_max = 1.0
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state
fileName = "sodShock.hdf5"


# ---------------------------------------------------
boxSize = x_max - x_min

glass_L = sw.load("glassCube_64.hdf5")
glass_R = sw.load("glassCube_32.hdf5")

pos_L = glass_L.gas.coordinates.value
pos_R = glass_R.gas.coordinates.value
h_L = glass_L.gas.smoothing_lengths.value
h_R = glass_R.gas.smoothing_lengths.value

radius_L = sqrt(
    (pos_L[:, 0] - 0.5) ** 2 + (pos_L[:, 1] - 0.5) ** 2 + (pos_L[:, 2] - 0.5) ** 2
)
index_L = radius_L < 0.25
pos_LL = pos_L[index_L]
h_LL = h_L[index_L]

radius_R = sqrt(
    (pos_R[:, 0] - 0.5) ** 2 + (pos_R[:, 1] - 0.5) ** 2 + (pos_R[:, 2] - 0.5) ** 2
)
index_R = radius_R > 0.25
pos_RR = pos_R[index_R]
h_RR = h_R[index_R]

# Merge things
pos = append(pos_LL, pos_RR, axis=0)
h = append(h_LL, h_RR)

numPart_L = size(h_LL)
numPart_R = size(h_RR)
numPart = size(h)

vol_L = 4.0 * pi / 3.0 * 0.25**3
vol_R = 1.0 - 4.0 * pi / 3.0 * 0.25**3

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)

for i in range(numPart):
    x = sqrt((pos[i, 0] - 0.5) ** 2 + (pos[i, 1] - 0.5) ** 2 + (pos[i, 2] - 0.5) ** 2)

    if x < 0.25:  # left
        u[i] = P_L / (rho_L * (gamma - 1.0))
        m[i] = rho_L * vol_L / numPart_L
        v[i, 0] = v_L
    else:  # right
        u[i] = P_R / (rho_R * (gamma - 1.0))
        m[i] = rho_R * vol_R / numPart_R
        v[i, 0] = v_R

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [1.0, 1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(uint64)

w.write(fileName)
