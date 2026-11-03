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
###############################################################################

import numpy as np
import h5py
import unyt
import swiftsimio as sw

# Generates an overdensity within a vacuum to test the vacuum resolving
# capabilities of the code

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index

fileName = "vacuum.hdf5"

# ---------------------------------------------------
glass = h5py.File("glassPlane_128.hdf5", "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:] * 0.3

# Shrink the glass to half its size, move it to the centre of the box, and
# clone it in the y direction
pos = 0.5 * pos + np.array([0.25, 0.0, 0.0])
h *= 0.5
pos = np.append(pos, pos + np.array([0.0, 0.5, 0.0]), axis=0)
h = np.append(h, h)

numPart = len(h)
vol = 0.5

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

m[:] = 1.0 * vol / numPart
u[:] = 1.0 / (1.0 * (gamma - 1.0))

# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)
