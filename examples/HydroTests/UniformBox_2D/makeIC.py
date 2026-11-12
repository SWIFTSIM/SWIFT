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
import sys
import numpy as np

# Generates a swift IC file containing a cartesian distribution of particles
# at a constant density and pressure in a square (2D) box

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1.0
L = int(sys.argv[1])  # Number of particles along one axis
rho = 1.0  # Density
P = 1.0  # Pressure
gamma = 5.0 / 3.0  # Gas adiabatic index
fileName = "uniformPlane.hdf5"


# ---------------------------------------------------
numPart = L**2
mass = boxSize**2 * rho / numPart
internalEnergy = P / ((gamma - 1.0) * rho)

# Generate particles
coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart)
h = np.zeros(numPart)
u = np.zeros(numPart)
ids = np.zeros(numPart, dtype="L")

for i in range(L):
    for j in range(L):
        index = i * L + j
        x = i * boxSize / L + boxSize / (2.0 * L)
        y = j * boxSize / L + boxSize / (2.0 * L)
        z = 0
        coords[index, 0] = x
        coords[index, 1] = y
        coords[index, 2] = z
        v[index, 0] = 0.0
        v[index, 1] = 0.0
        v[index, 2] = 0.0
        m[index] = mass
        h[index] = 1.23485 * boxSize / L
        u[index] = internalEnergy
        ids[index] = index


# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [boxSize, boxSize], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = coords * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = (ids + 1).astype(np.uint64)

w.write(fileName)
