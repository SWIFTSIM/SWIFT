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

# Generates a swift IC file for the 1D Noh problem test in a periodic box

# Parameters
numPart = 1001
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1.0  # Background density
P0 = 1.0e-6  # Background pressure
fileName = "noh.hdf5"

# ---------------------------------------------------
coords = np.zeros((numPart, 3))
h = np.zeros(numPart)
vol = 2.0

for i in range(numPart):
    coords[i, 0] = i * vol / numPart + vol / (2.0 * numPart)
    h[i] = 1.2348 * vol / numPart

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

m[:] = rho0 * vol / numPart
u[:] = P0 / (rho0 * (gamma - 1))
v[coords[:, 0] < vol / 2.0, 0] = 1
v[coords[:, 0] > vol / 2.0, 0] = -1

# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [vol], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
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
