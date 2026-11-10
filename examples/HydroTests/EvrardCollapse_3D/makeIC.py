################################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import unyt
import swiftsimio as sw
import argparse as ap
import numpy as np

parser = ap.ArgumentParser(
    description="Generates a swift IC file for the Evrard collapse"
)

parser.add_argument(
    "-n",
    "--nparts",
    help="""
         Number of particles to be used in the Evrard collapse.
         """,
    required=False,
    default=100000,
)

args = vars(parser.parse_args())

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
M = 1.0  # total mass of the sphere
R = 1.0  # radius of the sphere
u0 = 0.05 / M  # initial thermal energy
fileName = "evrard.hdf5"
numPart = int(args["nparts"])

r = R * np.sqrt(np.random.random(numPart))
phi = 2.0 * np.pi * np.random.random(numPart)
cos_theta = 2.0 * np.random.random(numPart) - 1.0

sin_theta = np.sqrt(1.0 - cos_theta**2)
cos_phi = np.cos(phi)
sin_phi = np.sin(phi)

pos = np.zeros((numPart, 3))
pos[:, 0] = r * sin_theta * cos_phi
pos[:, 1] = r * sin_theta * sin_phi
pos[:, 2] = r * cos_theta

# shift particles to put the sphere in the centre of the box
pos += np.array([50.0 * R, 50.0 * R, 50.0 * R])

h = np.ones(numPart) * 2.0 * R / numPart ** (1.0 / 3.0)

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.ones(numPart) * M / numPart
u = np.ones(numPart) * u0

# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [100.0 * R, 100.0 * R, 100.0 * R],
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
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)
