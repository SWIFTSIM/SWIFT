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
import numpy as np

# Generates a swift IC file for the Square test in a periodic box

# Parameters
L = 64  # Number of particles on the side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 4.0  # Gas central density
rho1 = 1.0  # Gas outskirt density
P0 = 2.5  # Gas central pressure
P1 = 2.5  # Gas central pressure
vx = 0.0  # Random velocity for all particles
vy = 0.0
fileOutputName = "square.hdf5"
# ---------------------------------------------------

vol = 1.0

numPart_out = L * L
numPart_in = int(L * L * rho0 / rho1 / 4)

L_out = int(np.sqrt(numPart_out))
L_in = int(np.sqrt(numPart_in))

pos_out = np.zeros((numPart_out, 3))
for i in range(L_out):
    for j in range(L_out):
        index = i * L_out + j
        pos_out[index, 0] = i / (float(L_out)) + 1.0 / (2.0 * L_out)
        pos_out[index, 1] = j / (float(L_out)) + 1.0 / (2.0 * L_out)
h_out = np.ones(numPart_out) * (1.0 / L_out) * 1.2348
m_out = np.ones(numPart_out) * vol * rho1 / numPart_out
u_out = np.ones(numPart_out) * P1 / (rho1 * (gamma - 1.0))

pos_in = np.zeros((numPart_in, 3))
for i in range(L_in):
    for j in range(L_in):
        index = i * L_in + j
        pos_in[index, 0] = 0.25 + i / float(2.0 * L_in) + 1.0 / (2.0 * 2.0 * L_in)
        pos_in[index, 1] = 0.25 + j / float(2.0 * L_in) + 1.0 / (2.0 * 2.0 * L_in)
h_in = np.ones(numPart_in) * (1.0 / L_in) * 1.2348
m_in = np.ones(numPart_in) * 0.25 * vol * rho0 / numPart_in
u_in = np.ones(numPart_in) * P0 / (rho0 * (gamma - 1.0))

# Remove the central particles
select_out = np.logical_or(
    np.logical_or(pos_out[:, 0] < 0.25, pos_out[:, 0] > 0.75),
    np.logical_or(pos_out[:, 1] < 0.25, pos_out[:, 1] > 0.75),
)
pos_out = pos_out[select_out, :]
h_out = h_out[select_out]
u_out = u_out[select_out]
m_out = m_out[select_out]

# Add the central region
pos = np.append(pos_out, pos_in, axis=0)
h = np.append(h_out, h_in, axis=0)
u = np.append(u_out, u_in)
m = np.append(m_out, m_in)
numPart = np.size(h)
ids = np.linspace(1, numPart, numPart)
vel = np.zeros((numPart, 3))
vel[:, 0] = vx
vel[:, 1] = vy


# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [vol, vol], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = vel * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)
