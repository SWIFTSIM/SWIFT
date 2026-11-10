###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Generates a swift IC file containing a perturbed cartesian distribution of particles
# at a constant density and pressure in a cubic box

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1.0
L = int(sys.argv[1])  # Number of particles along one axis
rho = 1.0  # Density
P = 1.0  # Pressure
gamma = 5.0 / 3.0  # Gas adiabatic index
pert = 0.01  # Perturbation scale (in units of the interparticle separation)
fileName = "perturbedBox.hdf5"


# ---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.0) * rho)

# Generate particles
coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros((numPart, 1))
h = np.zeros((numPart, 1))
u = np.zeros((numPart, 1))
ids = np.zeros((numPart, 1), dtype="L")

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i * L * L + j * L + k
            x = (
                i * boxSize / L
                + boxSize / (2 * L)
                + np.random.random() * pert * boxSize / (2.0 * L)
            )
            y = (
                j * boxSize / L
                + boxSize / (2 * L)
                + np.random.random() * pert * boxSize / (2.0 * L)
            )
            z = (
                k * boxSize / L
                + boxSize / (2 * L)
                + np.random.random() * pert * boxSize / (2.0 * L)
            )
            coords[index, 0] = x
            coords[index, 1] = y
            coords[index, 2] = z
            v[index, 0] = 0.0
            v[index, 1] = 0.0
            v[index, 2] = 0.0
            m[index] = mass
            h[index] = 1.1255 * boxSize / L
            u[index] = internalEnergy
            ids[index] = index


# --------------------------------------------------

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
w.gas.masses = m[:, 0] * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h[:, 0] * unyt.cm * w.a.comoving
w.gas.internal_energy = u[:, 0] * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = (ids[:, 0] + 1).astype(np.uint64)

w.write(fileName)
