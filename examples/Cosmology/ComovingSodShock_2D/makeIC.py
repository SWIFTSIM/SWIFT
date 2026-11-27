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
import numpy as np

# Generates a swift IC file for the 2D Sod Shock in a periodic box

unit_l_in_cgs = 3.086e18
unit_m_in_cgs = 2.94e55
unit_t_in_cgs = 3.086e18

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
x_min = -1.0
x_max = 1.0
rho_L = 1.0  # Density left state
rho_R = 0.140625  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state
a_beg = 0.001
fileName = "sodShock.hdf5"


# ---------------------------------------------------
boxSize = x_max - x_min

glass_L = sw.load("glassPlane_128.hdf5")
glass_R = sw.load("glassPlane_48.hdf5")

pos_L = glass_L.gas.coordinates.value * 0.5
pos_R = glass_R.gas.coordinates.value * 0.5
h_L = glass_L.gas.smoothing_lengths.value * 0.5
h_R = glass_R.gas.smoothing_lengths.value * 0.5

# Merge things
aa = pos_L - np.array([0.5, 0.0, 0.0])
pos_LL = np.append(pos_L, pos_L + np.array([0.5, 0.0, 0.0]), axis=0)
pos_RR = np.append(pos_R, pos_R + np.array([0.5, 0.0, 0.0]), axis=0)
pos = np.append(pos_LL - np.array([1.0, 0.0, 0.0]), pos_RR, axis=0)
h_LL = np.append(h_L, h_L)
h_RR = np.append(h_R, h_R)
h = np.append(h_LL, h_RR)

numPart_L = np.size(h_LL)
numPart_R = np.size(h_RR)
numPart = np.size(h)

vol_L = 0.5
vol_R = 0.5

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

for i in range(numPart):
    x = pos[i, 0]

    if x < 0:  # left
        u[i] = P_L / (rho_L * (gamma - 1.0))
        m[i] = rho_L * vol_L / numPart_L
        v[i, 0] = v_L
    else:  # right
        u[i] = P_R / (rho_R * (gamma - 1.0))
        m[i] = rho_R * vol_R / numPart_R
        v[i, 0] = v_R

# Shift particles
pos[:, 0] -= x_min

u /= a_beg ** (3.0 * (gamma - 1.0))

# File
UL = unyt.Unit(unit_l_in_cgs * unyt.cm)
UM = unyt.Unit(unit_m_in_cgs * unyt.g)
Ut = unyt.Unit(unit_t_in_cgs * unyt.s)
units = unyt.UnitSystem("comoving_sod_2d_units", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize, 0.5], UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)
