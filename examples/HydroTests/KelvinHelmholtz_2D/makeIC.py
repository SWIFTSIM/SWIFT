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

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
L2 = 256  # Particles along one edge in the low-density region
gamma = 5.0 / 3.0  # Gas adiabatic index
P1 = 2.5  # Central region pressure
P2 = 2.5  # Outskirts pressure
v1 = 0.5  # Central region velocity
v2 = -0.5  # Outskirts vlocity
rho1 = 2  # Central density
rho2 = 1  # Outskirts density
omega0 = 0.1
sigma = 0.05 / np.sqrt(2)
fileOutputName = "kelvinHelmholtz.hdf5"
# ---------------------------------------------------

# Start by generating grids of particles at the two densities
numPart2 = L2 * L2
L1 = int(np.sqrt(numPart2 / rho2 * rho1))
numPart1 = L1 * L1

print("N2 =", numPart2, "N1 =", numPart1)
print("L2 =", L2, "L1 = ", L1)
print("rho2 =", rho2, "rho1 =", (float(L1 * L1)) / (float(L2 * L2)))

coords1 = np.zeros((numPart1, 3))
coords2 = np.zeros((numPart2, 3))
h1 = np.ones(numPart1) * 1.2348 / L1
h2 = np.ones(numPart2) * 1.2348 / L2
m1 = np.zeros(numPart1)
m2 = np.zeros(numPart2)
u1 = np.zeros(numPart1)
u2 = np.zeros(numPart2)
vel1 = np.zeros((numPart1, 3))
vel2 = np.zeros((numPart2, 3))

# Particles in the central region
for i in range(L1):
    for j in range(L1):

        index = i * L1 + j

        x = i / float(L1) + 1.0 / (2.0 * L1)
        y = j / float(L1) + 1.0 / (2.0 * L1)

        coords1[index, 0] = x
        coords1[index, 1] = y
        u1[index] = P1 / (rho1 * (gamma - 1.0))
        vel1[index, 0] = v1

# Particles in the outskirts
for i in range(L2):
    for j in range(L2):

        index = i * L2 + j

        x = i / float(L2) + 1.0 / (2.0 * L2)
        y = j / float(L2) + 1.0 / (2.0 * L2)

        coords2[index, 0] = x
        coords2[index, 1] = y
        u2[index] = P2 / (rho2 * (gamma - 1.0))
        vel2[index, 0] = v2


# Now concatenate arrays
where1 = np.abs(coords1[:, 1] - 0.5) < 0.25
where2 = np.abs(coords2[:, 1] - 0.5) > 0.25

coords = np.append(coords1[where1, :], coords2[where2, :], axis=0)

# print L2*(L2/2), L1*(L1/2)
# print shape(coords), shape(coords1[where1,:]), shape(coords2[where2,:])
# print shape(coords), shape(logical_not(coords1[where1,:])), shape(logical_not(coords2[where2,:]))

vel = np.append(vel1[where1, :], vel2[where2, :], axis=0)
h = np.append(h1[where1], h2[where2], axis=0)
m = np.append(m1[where1], m2[where2], axis=0)
u = np.append(u1[where1], u2[where2], axis=0)
numPart = np.size(h)
ids = np.linspace(1, numPart, numPart)
m[:] = (0.5 * rho1 + 0.5 * rho2) / float(numPart)

# Velocity perturbation
vel[:, 1] = (
    omega0
    * np.sin(4 * np.pi * coords[:, 0])
    * (
        np.exp(-((coords[:, 1] - 0.25) ** 2) / (2 * sigma**2))
        + np.exp(-((coords[:, 1] - 0.75) ** 2) / (2 * sigma**2))
    )
)

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [1.0, 1.0], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=2,
)

w.gas.coordinates = coords * unyt.cm * w.a.comoving
w.gas.velocities = vel * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)
