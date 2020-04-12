###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import h5py as h5
import numpy as np
import math

num_snapshots = 61

# physical constants in cgs
NEWTON_GRAVITY_CGS = 6.67430e-8
SOLAR_MASS_IN_CGS = 1.98841e33
PARSEC_IN_CGS = 3.08567758149e18
YEAR_IN_CGS = 3.15569251e7

# choice of units
const_unit_length_in_cgs = PARSEC_IN_CGS
const_unit_mass_in_cgs = SOLAR_MASS_IN_CGS
const_unit_velocity_in_cgs = 1e5

# parameters of potential
surface_density = 10.0
scale_height = 100.0
centre = np.array([300.0, 300.0, 300.0])

# constants
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_G = (
    NEWTON_GRAVITY_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)


t = np.zeros(num_snapshots)
E_k = np.zeros(num_snapshots)
E_p = np.zeros(num_snapshots)
E_tot = np.zeros(num_snapshots)

# loop over the snapshots
for i in range(num_snapshots):

    filename = "Disc-Patch_%04d.hdf5" % i
    f = h5.File(filename, "r")

    # time
    t[i] = f["/Header"].attrs.get("Time")[0]

    # positions and velocities of the particles
    p = f["/PartType1/Coordinates"][:]
    v = f["/PartType1/Velocities"][:]

    # Kinetic energy
    E_k[i] = np.sum(0.5 * (v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2))

    # Potential energy
    d = p[:, 0] - centre[0]
    E_p[i] = np.sum(
        2.0
        * math.pi
        * const_G
        * surface_density
        * scale_height
        * np.log(np.cosh(np.abs(d) / scale_height))
    )

    # Total energy
    E_tot[i] = E_k[i] + E_p[i]


# Maximal change
delta_E = np.max(np.abs(E_tot - E_tot[0])) / E_tot[0]

print("Maximal relative energy change   :", delta_E * 100, "%")
