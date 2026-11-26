###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a regular 1D grid that contains the Woodward & Colella (1984)
# interacting blast waves test.
# Note that the original test requires reflective boundary conditions. Since we
# do not have those in SWIFT, we double the computational domain and mirror the
# setup in the second half.

import numpy as np
import h5py
import unyt
import swiftsimio as sw

fileName = "interactingBlastWaves.hdf5"
numPart = 800
boxSize = 2.0

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart) + boxSize / numPart
h = np.zeros(numPart) + 2.0 * boxSize / numPart
u = np.zeros(numPart)
ids = np.arange(numPart, dtype="L")
rho = np.ones(numPart)

for i in range(numPart):
    coords[i, 0] = (i + 0.5) * boxSize / numPart
    if coords[i, 0] < 0.1 or coords[i, 0] > 1.9:
        u[i] = 2500.0
    elif coords[i, 0] > 0.9 and coords[i, 0] < 1.1:
        u[i] = 250.0
    else:
        u[i] = 0.025

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [boxSize], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=1,
)

w.gas.coordinates = coords * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
# SWIFT does not accept a particle ID of 0
w.gas.particle_ids = ids + 1

w.write(fileName)

# Add the extra fields not handled by the Writer
with h5py.File(fileName, "r+") as f:
    f["PartType0"].create_dataset("Density", data=rho, dtype="f")
