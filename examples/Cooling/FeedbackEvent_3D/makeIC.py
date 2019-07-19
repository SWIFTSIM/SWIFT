###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

from unyt import cm, kpc, mh, msun, K, s, kb

import h5py
import numpy as np

# Generates a SWIFT IC file for the Feedback blast test.

# Parameters
gamma = 5.0 / 3.0
initial_density = 0.1 * mh / (cm ** 3)
initial_temperature = 1e4 * K
inject_temperature = 10 ** (7.5) * K
mu = 0.5
particle_mass = 1e6 * msun

unit_system = cosmo_units
file_name = "feedback.hdf5"

# Read in glass file
with h5py.File("glassCube_32.hdf5", "r") as handle:
    positions = handle["/PartType0/Coordinates"][:]
    h = handle["PartType0/SmoothingLength"][:] * 0.3

number_of_particles = len(h)
side_length = (number_of_particles * particle_mass / initial_density) ** (1 / 3)
side_length.convert_to_base(unit_system)

print(f"Your box has a side length of {side_length}")

# Find the central particle
central_particle = np.sum((positions - 0.5) ** 2, axis=1).argmin()

# Inject the feedback into that central particle
background_internal_energy = (
    (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * initial_temperature
)
heated_internal_energy = (1.0 / (mu * mh)) * (kb / (gamma - 1.0)) * inject_temperature
internal_energy = np.ones_like(h) * background_internal_energy
internal_energy[central_particle] = heated_internal_energy

# Now we have all the information we need to set up the initial conditions!
output = Writer(unit_system=unit_system, box_size=side_length)

output.gas.coordinates = positions * side_length
output.gas.velocities = np.zeros_like(positions) * cm / s
output.gas.smoothing_length = h * side_length
output.gas.internal_energy = internal_energy
output.gas.masses = np.ones_like(h) * particle_mass

output.write(file_name)
