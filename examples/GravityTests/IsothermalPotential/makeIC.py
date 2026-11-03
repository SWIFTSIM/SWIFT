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

# Generates N particles in a spherical distribution centred on [0,0,0], to be moved in an isothermal potential
# usage: python3 makeIC.py 1000 0 : generate 1000 particles on circular orbits
#        python3 makeIC.py 1000 1 : generate 1000 particles with Lz/L uniform in [0,1]
# all particles move in the xy plane, and start at y=0

# physical constants in cgs
NEWTON_GRAVITY_CGS = 6.67408e-8
SOLAR_MASS_IN_CGS = 1.98848e33
PARSEC_IN_CGS = 3.08567758e18
YEAR_IN_CGS = 3.15569252e7

# choice of units
const_unit_length_in_cgs = 1000 * PARSEC_IN_CGS
const_unit_mass_in_cgs = SOLAR_MASS_IN_CGS
const_unit_velocity_in_cgs = 1e5

print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs)
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)


# rotation speed of isothermal potential [km/s]
vrot_kms = 200.0

# derived units
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_G = (
    NEWTON_GRAVITY_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)
print("G=", const_G)
vrot = vrot_kms * 1e5 / const_unit_velocity_in_cgs


# Parameters
periodic = 1  # 1 For periodic box
boxSize = 400.0  #  [kpc]
Radius = 100.0  # maximum radius of particles [kpc]
G = const_G

N = int(sys.argv[1])  # Number of particles
icirc = int(
    sys.argv[2]
)  # if = 0, all particles are on circular orbits, if = 1, Lz/Lcirc uniform in ]0,1[
L = N ** (1.0 / 3.0)

fileName = "Isothermal.hdf5"


# ---------------------------------------------------
numPart = N
mass = 1

# --------------------------------------------------

# set seed for random number
numpy.random.seed(1234)

# generate particle positions
radius = Radius * (numpy.random.rand(N)) ** (1.0 / 3.0)
ctheta = -1.0 + 2 * numpy.random.rand(N)
stheta = numpy.sqrt(1.0 - ctheta**2)
phi = 2 * math.pi * numpy.random.rand(N)
r = numpy.zeros((numPart, 3))
r[:, 0] = radius

#
speed = vrot
v = numpy.zeros((numPart, 3))
omega = speed / radius
period = 2.0 * math.pi / omega
print("period = minimum = ", min(period), " maximum = ", max(period))

omegav = omega
if icirc != 0:
    omegav = omega * numpy.random.rand(N)

v[:, 0] = -omegav * r[:, 1]
v[:, 1] = omegav * r[:, 0]

m = numpy.full((numPart,), mass, dtype="f")

ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)

# File
UL = unyt.Unit(const_unit_length_in_cgs * unyt.cm)
UM = unyt.Unit(const_unit_mass_in_cgs * unyt.g)
Ut = unyt.Unit(const_unit_time_in_cgs * unyt.s)
units = unyt.UnitSystem("isothermal_potential", UL, UM, Ut)

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

# Particle group
w.dark_matter.coordinates = r * UL * w.a.comoving
w.dark_matter.velocities = v * UL / Ut * w.a.comoving**0
w.dark_matter.masses = m * UM * w.a.comoving**0
w.dark_matter.particle_ids = ids.astype(numpy.uint64)

w.write(fileName)
