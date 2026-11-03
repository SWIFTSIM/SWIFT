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

import h5py
import unyt
import swiftsimio as sw
from numpy import *

# Generates a swift IC file for the 3D Noh problem in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1.0  # Background density
P0 = 1.0e-6  # Background pressure
fileName = "noh.hdf5"

# ---------------------------------------------------
glass = h5py.File("glassCube_64.hdf5", "r")

vol = 8.0

pos = glass["/PartType0/Coordinates"][:, :] * vol ** (1.0 / 3.0)
h = glass["/PartType0/SmoothingLength"][:] * vol ** (1.0 / 3.0)
numPart = size(h)

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)

m = zeros(numPart)
u = zeros(numPart)
m[:] = rho0 * vol / numPart
u[:] = P0 / (rho0 * (gamma - 1))

# Make radial velocities
# r = sqrt((pos[:,0]-1.)**2 + (pos[:,1]-1.)**2)
# theta = arctan2((pos[:,1]-1.), (pos[:,0]-1.))
v[:, 0] = -(pos[:, 0] - 1.0)
v[:, 1] = -(pos[:, 1] - 1.0)
v[:, 2] = -(pos[:, 2] - 1.0)

norm_v = sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
v[:, 0] /= norm_v
v[:, 1] /= norm_v
v[:, 2] /= norm_v


# Densities (each particle has the rho0 value as density, since each mass
# is rho0 * vol / numpart. The vol / numpart can be seen as the particle value,
# i.e. the whole volume is split into numpart)
densities = np.ones(numPart, dtype=np.float32)

# Material IDs
mat_ids = np.zeros(numPart)

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [vol ** (1.0 / 3.0), vol ** (1.0 / 3.0), vol ** (1.0 / 3.0)],
        unyt.cm,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(uint64)

w.write(fileName)

# Add the extra fields not handled by the Writer
with h5py.File(fileName, "r+") as f:
    f["PartType0"].create_dataset("MaterialIDs", data=mat_ids, dtype="L")
    f["PartType0"].create_dataset("Density", data=densities, dtype="f")
