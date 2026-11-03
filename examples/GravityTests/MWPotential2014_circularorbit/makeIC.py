###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
import numpy as np
import unyt
import swiftsimio as sw

box_size = 1000.0
N_PARTICLES = 3
print("Initial conditions written to 'circular_orbits_MW.hdf5'")

pos = np.zeros((3, 3))
pos[0, 2] = 5.0
pos[1, 1] = 5.0
pos[2, 0] = 30.0
pos += np.array(
    [box_size / 2, box_size / 2, box_size / 2]
)  # shifts the particles to the center of the box
vel = np.zeros((3, 3))
vel[0, 0] = 198.5424557586175
vel[1, 0] = 225.55900735974072
vel[2, 1] = 188.5272441374569

ids = np.array([1.0, 2.0, 3.0])
mass = np.array([1.0, 1.0, 1.0]) * 1e-10

# File
UL = unyt.Unit(3.086e21 * unyt.cm)
UM = unyt.Unit(1.98848e43 * unyt.g)
Ut = unyt.Unit(3.086e16 * unyt.s)
units = unyt.UnitSystem("mw_potential_2014_circular_orbit", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [box_size, box_size, box_size],
        UL,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

# Particle group
w.dark_matter.coordinates = pos * UL * w.a.comoving
w.dark_matter.velocities = vel * UL / Ut * w.a.comoving**0
w.dark_matter.masses = mass * UM * w.a.comoving**0
w.dark_matter.particle_ids = ids.astype(np.uint64)

w.write("circular_orbits_MW.hdf5")
