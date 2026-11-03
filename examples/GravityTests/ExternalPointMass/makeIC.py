###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 John A. Regan (john.a.regan@durham.ac.uk)
#                    Tom Theuns (tom.theuns@durham.ac.uk)
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

import sys
import numpy
import unyt
import swiftsimio as sw
import math
import random

# Generates a random distriution of particles, for motion in an external potential centred at (0,0,0)

# physical constants in cgs
NEWTON_GRAVITY_CGS = 6.67408e-8
SOLAR_MASS_IN_CGS = 1.98848e33
PARSEC_IN_CGS = 3.08567758e18

# choice of units
const_unit_length_in_cgs = 1000 * PARSEC_IN_CGS
const_unit_mass_in_cgs = SOLAR_MASS_IN_CGS
const_unit_velocity_in_cgs = 1e5

print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs)
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)
print("UnitTime_in_cgs:     ", const_unit_length_in_cgs / const_unit_velocity_in_cgs)

# derived units
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_G = (
    NEWTON_GRAVITY_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)
print("---------------------")
print("G in internal units: ", const_G)


# Parameters
periodic = 1  # 1 For periodic box
boxSize = 100.0  #
max_radius = boxSize / 4.0  # maximum radius of particles
Mass = 1e10
print("Mass at the centre:  ", Mass)

numPart = int(sys.argv[1])  # Number of particles
mass = 1.0

fileName = "PointMass.hdf5"


# --------------------------------------------------

# generate particle positions
radius = max_radius * (numpy.random.rand(numPart)) ** (1.0 / 3.0)
print("---------------------")
print("Radius: minimum = ", min(radius))
print("Radius: maximum = ", max(radius))
radius = numpy.sort(radius)
r = numpy.zeros((numPart, 3))
r[:, 0] = radius

# generate particle velocities
speed = numpy.sqrt(const_G * Mass / radius)
omega = speed / radius
period = 2.0 * math.pi / omega
print("---------------------")
print("Period: minimum = ", min(period))
print("Period: maximum = ", max(period))

v = numpy.zeros((numPart, 3))
v[:, 0] = -omega * r[:, 1]
v[:, 1] = omega * r[:, 0]

m = numpy.full((numPart,), mass)

ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)

# File
UL = unyt.Unit(3.0856776e21 * unyt.cm)
UM = unyt.Unit(1.98848e33 * unyt.g)
Ut = unyt.Unit(3.0856776e16 * unyt.s)
units = unyt.UnitSystem("external_point_mass", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize, boxSize, boxSize],
        UL,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

w.dark_matter.coordinates = r * UL * w.a.comoving
w.dark_matter.velocities = v * UL / Ut * w.a.comoving**0
w.dark_matter.masses = m * UM * w.a.comoving**0
w.dark_matter.particle_ids = ids.astype(numpy.uint64)

w.write(fileName)
