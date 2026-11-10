################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

import h5py
import unyt
import swiftsimio as sw
import numpy as np

# Generates a swift IC file for the Gresho-Chan vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
rho0 = 1  # Gas density
P0 = 0.0  # Constant additional pressure (should have no impact on the dynamics)
fileOutputName = "greshoVortex.hdf5"
fileGlass = "glassCube_64.hdf5"
# ---------------------------------------------------

# Get position and smoothing lengths from the glass
fileInput = sw.load(fileGlass)
coords = fileInput.gas.coordinates.value
h = fileInput.gas.smoothing_lengths.value
ids = fileInput.gas.particle_ids.value
boxSize = fileInput.metadata.boxsize[0].value
numPart = np.size(h)

# Now generate the rest
m = np.ones(numPart) * rho0 * boxSize**3 / numPart
u = np.zeros(numPart)
v = np.zeros((numPart, 3))

# Add rho and material ids
rho = np.ones(numPart)  # given rho0=1, this is easily done
material_ids = np.zeros(numPart)  # assuming ideal gas particles

for i in range(numPart):

    x = coords[i, 0]
    y = coords[i, 1]

    r2 = (x - boxSize / 2) ** 2 + (y - boxSize / 2) ** 2
    r = np.sqrt(r2)

    v_phi = 0.0
    if r < 0.2:
        v_phi = 5.0 * r
    elif r < 0.4:
        v_phi = 2.0 - 5.0 * r
    else:
        v_phi = 0.0
    v[i, 0] = -v_phi * (y - boxSize / 2) / r
    v[i, 1] = v_phi * (x - boxSize / 2) / r
    v[i, 2] = 0.0

    P = P0
    if r < 0.2:
        P = P + 5.0 + 12.5 * r2
    elif r < 0.4:
        P = P + 9.0 + 12.5 * r2 - 20.0 * r + 4.0 * np.log(r / 0.2)
    else:
        P = P + 3.0 + 4.0 * np.log(2.0)
    u[i] = P / ((gamma - 1.0) * rho0)


# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [boxSize, boxSize, boxSize],
        unyt.cm,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

w.gas.coordinates = coords * unyt.cm * w.a.comoving
w.gas.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)

# Add the extra fields not handled by the Writer
with h5py.File(fileOutputName, "r+") as f:
    f["PartType0"].create_dataset("MaterialIDs", data=material_ids, dtype="i")
    f["PartType0"].create_dataset("Density", data=rho, dtype="f")
