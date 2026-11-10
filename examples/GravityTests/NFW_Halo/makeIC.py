################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Ashley Kelly ()
#                    Folkert Nobels (nobels@strw.leidenuniv.nl)
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

C = 8.0
M_200 = 2.0


print("\nInitial conditions written to 'test_nfw.hdf5'")

pos = np.array([[8.0, 0.0, 0.0]]) + 500.0
vel = np.array([[0.0, 240.0, 5.0]])
ids = np.array([1.0])
mass = np.array([1.0])

# File
UL = unyt.Unit(3.086e21 * unyt.cm)
UM = unyt.Unit(1.988e33 * unyt.g)
Ut = unyt.Unit(3.086e16 * unyt.s)
unit_system = unyt.UnitSystem("nfw_halo", UL, UM, Ut)

w = sw.Writer(
    unit_system=unit_system,
    boxsize=sw.cosmo_array(
        [1000.0, 1000.0, 1000.0],
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

w.write("test_nfw.hdf5")
