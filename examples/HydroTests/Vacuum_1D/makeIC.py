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

# Generates an overdensity within a vacuum to test the vacuum handling
# capabilities of the code

import numpy as np
import h5py
import unyt
import swiftsimio as sw

fileName = "vacuum.hdf5"
numPart = 100
boxSize = 1.0
gamma = 5.0 / 3.0

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart)
h = np.zeros(numPart)
u = np.zeros(numPart)
ids = np.arange(1, numPart + 1, dtype="L")
rho = np.zeros(numPart)

# first set the positions, as we try to do a reasonable volume estimate to
# set the masses
for i in range(numPart):
    # we only generate particles in the range [0.25, 0.75]
    coords[i, 0] = 0.25 + 0.5 * (i + 0.5) / numPart
    rho[i] = 1.0
    P = 1.0
    u[i] = P / (gamma - 1.0) / rho[i]
    m[i] = rho[i] * 0.5 / numPart
    # reasonable smoothing length estimate
    h[i] = 1.0 / numPart

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
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)

# Extra fields not handled by the Writer
with h5py.File(fileName, "r+") as f:
    f["PartType0"].create_dataset("Density", data=rho, dtype="f")
