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

import unyt
import swiftsimio as sw
import numpy as np

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
coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
h = np.zeros(numPart)
u = np.zeros(numPart)

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
UL = unyt.Unit(unit_l_in_cgs * unyt.cm)
UM = unyt.Unit(unit_m_in_cgs * unyt.g)
Ut = unyt.Unit(unit_t_in_cgs * unyt.s)
units = unyt.UnitSystem("comoving_toro_units", UL, UM, Ut)

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
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)
