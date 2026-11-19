################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
from optparse import OptionParser
from astropy import units
from astropy import constants

# Fixed seed so that the initial conditions are reproducible
np.random.seed(1)


def parse_options():

    usage = "usage: %prog [options] file"
    parser = OptionParser(usage=usage)

    parser.add_option(
        "--lJ",
        action="store",
        dest="lJ",
        type="float",
        default=0.250,
        help="Jeans wavelength in box size unit",
    )

    parser.add_option(
        "--rho",
        action="store",
        dest="rho",
        type="float",
        default=0.1,
        help="Mean gas density in atom/cm3",
    )

    parser.add_option(
        "--mass",
        action="store",
        dest="mass",
        type="float",
        default=50,
        help="Gas particle mass in solar mass",
    )

    parser.add_option(
        "--level",
        action="store",
        dest="level",
        type="int",
        default=5,
        help="Resolution level: N = (2**l)**3",
    )

    parser.add_option(
        "-o",
        action="store",
        dest="outputfilename",
        type="string",
        default="box.dat",
        help="output filename",
    )

    options, args = parser.parse_args()

    files = args

    return files, options


########################################
# main
########################################

files, opt = parse_options()

# define standard units
UnitMass_in_cgs = 1.989e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.085678e21  # kpc in centimeters
UnitVelocity_in_cgs = 1e5  # km/s in centimeters per second
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm
UnitTime = UnitTime_in_cgs * units.s
UnitVelocity = UnitVelocity_in_cgs * units.cm / units.s


# Number of particles
N = (2**opt.level) ** 3  # number of particles

# Mean density
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm**3

# Gas particle mass
m = opt.mass  # in solar mass
m = m * units.Msun

# Gas mass in the box
M = N * m

# Size of the box
L = (M / rho) ** (1 / 3.0)

# Jeans wavelength in box size unit
lJ = opt.lJ
lJ = lJ * L

# Gravitational constant
G = constants.G

# Jeans wave number
kJ = 2 * np.pi / lJ
# Velocity dispersion
sigma = np.sqrt(4 * np.pi * G * rho) / kJ


print("Number of particles                   : {}".format(N))
print("Equivalent velocity dispertion        : {}".format(sigma.to(units.m / units.s)))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
sigma = sigma.to(UnitVelocity).value

# Generate the particles
pos = np.random.random([N, 3]) * np.array([L, L, L])
vel = np.zeros((N, 3))
mass = np.ones(N) * m
u = np.ones(N) * sigma**2
ids = np.arange(1, N + 1)
h = np.ones(N) * 3 * L / N ** (1 / 3.0)

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))


# File
UL = unyt.Unit(UnitLength_in_cgs * unyt.cm)
UM = unyt.Unit(UnitMass_in_cgs * unyt.g)
Ut = unyt.Unit(UnitTime_in_cgs * unyt.s)
unit_system = unyt.UnitSystem("jeans_fragmentation", UL, UM, Ut)

w = sw.Writer(
    unit_system=unit_system,
    boxsize=sw.cosmo_array(
        [L, L, L], UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Particle group
w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = vel * UL / Ut * w.a.comoving**0
w.gas.masses = mass * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(opt.outputfilename)

print("{} saved.".format(opt.outputfilename))
