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

import numpy
import unyt
import swiftsimio as sw
import sys
import math

# Generates N particles in a box of [0:BoxSize, 0:BoxSize, 0:BoxSize], with the
# disc mid-plane at x = BoxSize/2 (i.e. up to 3 scale heights on either side).
#
# The potential is along the x-axis.
#
# see Creasey, Theuns & Bower, 2013, for the equations:
# disc parameters are: surface density sigma
#                      scale height b
# density: rho(x) = (sigma/2b) sech^2(x/b)
# isothermal velocity dispersion = <v_x^2> = b pi G sigma
# grad potential  = 2 pi G sigma tanh(x/b)
# potential       = 2 pi G sigma b ln(cosh(x/b)) + const
# Dynamical time  = sqrt(b / (G sigma))
# to obtain the 1/ch^2(x/b) profile from a uniform profile (a glass, say, or a uniform random variable), note that, when integrating in x
# \int 0^x dx/ch^2(x) = tanh(x)-tanh(0) = \int_0^u du = u (where the last integral refers to a uniform density distribution), so that x = atanh(u)
#
# usage: python3 makeIC.py 1000

# physical constants in cgs
NEWTON_GRAVITY_CGS = 6.67430e-8
SOLAR_MASS_IN_CGS = 1.98841e33
PARSEC_IN_CGS = 3.08567758149e18
YEAR_IN_CGS = 3.15569251e7

# choice of units
const_unit_length_in_cgs = PARSEC_IN_CGS
const_unit_mass_in_cgs = SOLAR_MASS_IN_CGS
const_unit_velocity_in_cgs = 1e5

print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs)
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)


# parameters of potential
surface_density = 10.0
scale_height = 100.0

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
v_disp = numpy.sqrt(scale_height * math.pi * const_G * surface_density)
t_dyn = numpy.sqrt(scale_height / (const_G * surface_density))
print("dynamical time = ", t_dyn)
print(" velocity dispersion = ", v_disp)

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 600.0  #
Radius = 100.0  # maximum radius of particles [pc] (not used)
G = const_G
mass = 1
numPart = int(sys.argv[1])  # Number of particles
fileName = "Disc-Patch.hdf5"

# --------------------------------------------------

# set seed for random number
numpy.random.seed(1234)

# generate particle positions
r = numpy.zeros((numPart, 3))
x = scale_height * numpy.arctanh(numpy.random.rand(2 * numPart))
gd = x < boxSize / 2
r[:, 0] = x[gd][0:numPart]
random = numpy.random.rand(numPart) > 0.5
r[random, 0] *= -1
r[:, 0] += 0.5 * boxSize

r[:, 1] = numpy.random.rand(numPart) * boxSize
r[:, 2] = numpy.random.rand(numPart) * boxSize

# generate particle velocities
v = numpy.zeros((numPart, 3))

m = numpy.ones((numPart,), dtype=numpy.float32) * mass

ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)

# File
UL = unyt.Unit(const_unit_length_in_cgs * unyt.cm)
UM = unyt.Unit(const_unit_mass_in_cgs * unyt.g)
Ut = unyt.Unit(const_unit_time_in_cgs * unyt.s)
units = unyt.UnitSystem("disc_patch_gravity_only", UL, UM, Ut)

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
