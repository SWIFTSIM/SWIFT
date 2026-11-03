###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a swift IC file for the 1D Sod Shock in a periodic box

unit_l_in_cgs = 3.086e18
unit_m_in_cgs = 2.94e55
unit_t_in_cgs = 3.086e18

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
numPart_L = 800  # Number of particles in the left state
x_min = -1.0
x_max = 1.0
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state
a_beg = 0.001
fileName = "sodShock.hdf5"


# ---------------------------------------------------

# Find how many particles we actually have
boxSize = x_max - x_min
numPart_R = int(numPart_L * (rho_R / rho_L))
numPart = numPart_L + numPart_R

# Now get the distances
delta_L = (boxSize / 2) / numPart_L
delta_R = (boxSize / 2) / numPart_R
offset_L = delta_L / 2
offset_R = delta_R / 2

# Build the arrays
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
h = zeros(numPart)
u = zeros(numPart)

# Set the particles on the left
for i in range(numPart_L):
    coords[i, 0] = x_min + offset_L + i * delta_L
    u[i] = P_L / (rho_L * (gamma - 1.0))
    h[i] = 1.2348 * delta_L
    m[i] = boxSize * rho_L / (2.0 * numPart_L)
    v[i, 0] = v_L

# Set the particles on the right
for j in range(numPart_R):
    i = numPart_L + j
    coords[i, 0] = offset_R + j * delta_R
    u[i] = P_R / (rho_R * (gamma - 1.0))
    h[i] = 1.2348 * delta_R
    m[i] = boxSize * rho_R / (2.0 * numPart_R)
    v[i, 0] = v_R

# Shift particles
coords[:, 0] -= x_min

u /= a_beg ** (3.0 * (gamma - 1.0))

# File
UL = unyt.Unit(unit_l_in_cgs * unyt.cm)
UM = unyt.Unit(unit_m_in_cgs * unyt.g)
Ut = unyt.Unit(unit_t_in_cgs * unyt.s)
units = unyt.UnitSystem("comoving_sod_1d_units", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize], UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=1,
)

w.gas.coordinates = coords * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(uint64)

w.write(fileName)
