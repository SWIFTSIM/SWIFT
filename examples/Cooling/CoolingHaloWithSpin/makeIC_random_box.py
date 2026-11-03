###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
import unyt
import swiftsimio as sw
import numpy as np
import math
import random

# Generates N particles in a spherically symmetric distribution with density profile ~r^(-2)
# usage: python3 makeIC.py 1000: generate 1000 particles

# Some constants

OMEGA = 0.3  # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
h = 0.67777  # hubble parameter
gamma = 5.0 / 3.0
eta = 1.2349

# First set unit velocity and then the circular velocity parameter for the isothermal potential
const_unit_velocity_in_cgs = 1.0e5  # kms^-1

v_c = 200.0
v_c_cgs = v_c * const_unit_velocity_in_cgs

# Now we use this to get the virial mass and virial radius, which we will set to be the unit mass and radius

# Find H_0, the inverse Hubble time, in cgs

H_0_cgs = 100.0 * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# From this we can find the virial radius, the radius within which the average density of the halo is
# 200. * the mean matter density

r_vir_cgs = v_c_cgs / (10.0 * H_0_cgs * np.sqrt(OMEGA))

# Now get the virial mass

M_vir_cgs = r_vir_cgs * v_c_cgs**2 / CONST_G_CGS

# Now set the unit length and mass

const_unit_mass_in_cgs = M_vir_cgs
const_unit_length_in_cgs = r_vir_cgs

print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs)
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)

# derived quantities
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
print("UnitTime_in_cgs:     ", const_unit_time_in_cgs)
const_G = (
    CONST_G_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)
print("G=", const_G)

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 4.0
G = const_G
N = int(sys.argv[1])  # Number of particles

# Create the file
filename = "random_box.hdf5"

# Units
UL = unyt.Unit(const_unit_length_in_cgs * unyt.cm)
UM = unyt.Unit(const_unit_mass_in_cgs * unyt.g)
Ut = unyt.Unit(const_unit_time_in_cgs * unyt.s)
units = unyt.UnitSystem("cooling_halo_spin_random_box_units", UL, UM, Ut)

# Header
w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize] * 3, UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# set seed for random number
np.random.seed(1234)

# Particle group
# Positions

# distribute particles randomly in the box
coords = np.zeros((N, 3))
coords[:, 0] = boxSize * np.random.rand(N)
coords[:, 1] = boxSize * np.random.rand(N)
coords[:, 2] = boxSize * np.random.rand(N)
# shift to centre of box
# coords += np.full((N,3),boxSize/2.)
print("x range = (%f,%f)" % (np.min(coords[:, 0]), np.max(coords[:, 0])))
print("y range = (%f,%f)" % (np.min(coords[:, 1]), np.max(coords[:, 1])))
print("z range = (%f,%f)" % (np.min(coords[:, 2]), np.max(coords[:, 2])))

print(np.mean(coords[:, 0]))
print(np.mean(coords[:, 1]))
print(np.mean(coords[:, 2]))

w.gas.coordinates = coords * UL * w.a.comoving

# All velocities set to zero
v = np.zeros((N, 3))
w.gas.velocities = v * UL / Ut * w.a.comoving**0

# All particles of equal mass
mass = 1.0 / N
m = np.full((N,), mass)
w.gas.masses = m * UM * w.a.comoving**0

# Smoothing lengths
l = (boxSize**3 / N) ** (1.0 / 3.0)  # local mean inter-particle separation
h = np.full((N,), eta * l)
w.gas.smoothing_lengths = h * UL * w.a.comoving

# Internal energies
u = v_c**2 / (2.0 * (gamma - 1.0))
u = np.full((N,), u)
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0

# Particle IDs
ids = 1 + np.linspace(0, N, N, endpoint=False, dtype="L")
w.gas.particle_ids = ids.astype(np.uint64)

w.write(filename)
