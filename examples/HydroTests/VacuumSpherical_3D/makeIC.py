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
import unyt
import swiftsimio as sw

# Generates an overdensity within a vacuum to test the vacuum resolving
# capabilities of the code

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index

fileName = "vacuum.hdf5"

# ---------------------------------------------------
glass = sw.load("glassCube_64.hdf5")

# Read particle positions and h from the glass
pos = glass.gas.coordinates.value
h = glass.gas.smoothing_lengths.value * 0.3

# Make 8 copies of the glass to get more particles
pos *= 0.5
h *= 0.5
pos = np.append(pos, pos + np.array([0.5, 0.0, 0.0]), axis=0)
pos = np.append(pos, pos + np.array([0.0, 0.5, 0.0]), axis=0)
pos = np.append(pos, pos + np.array([0.0, 0.0, 0.5]), axis=0)
h = np.append(h, h)
h = np.append(h, h)
h = np.append(h, h)

radius = np.sqrt(
    (pos[:, 0] - 0.5) ** 2 + (pos[:, 1] - 0.5) ** 2 + (pos[:, 2] - 0.5) ** 2
)
index = radius < 0.25
pos = pos[index]
h = h[index]

numPart = len(h)
vol = 4.0 * np.pi / 3.0 * 0.25**3

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
        [1.0, 1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)
