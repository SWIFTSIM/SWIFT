###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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

import h5py
import sys
import numpy
import math
import random
import numpy as np

# Generates N particles in a spherical distribution centred on [0,0,0], to be moved in an isothermal potential
# usage: python makeIC.py 1000 0 : generate 1000 particles on circular orbits
#        python makeIC.py 1000 1 : generate 1000 particles with Lz/L uniform in [0,1]
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


# Properties of the Hernquist potential
Mass = 1e15
scaleLength = 30.0  # kpc


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


def hernquistcircvel(r, M, a):
    """ Function that calculates the circular velocity in a 
    Hernquist potential.
    @param r: radius from centre of potential
    @param M: mass of the Hernquist potential
    @param a: Scale length of the potential
    @return: circular velocity
    """
    return (const_G * M * r) ** 0.5 / (r + a)


# Parameters
periodic = 1  # 1 For periodic box
boxSize = 400.0  #  [kpc]
Radius = 100.0  # maximum radius of particles [kpc]
G = const_G

N = 5
L = N ** (1.0 / 3.0)

fileName = "Hernquist.hdf5"


# ---------------------------------------------------
numPart = N
mass = 1

# --------------------------------------------------

# File
file = h5py.File(fileName, "w")

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = (
    const_unit_length_in_cgs / const_unit_velocity_in_cgs
)
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [0, numPart, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, numPart, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# set seed for random number
numpy.random.seed(1234)

# Particle group
grp1 = file.create_group("/PartType1")
# generate particle positions
# radius = Radius * (numpy.random.rand(N))**(1./3.) + 10.
radius = np.zeros(N)
radius[0] = 10
radius[1] = 20
radius[2] = 30
radius[3] = 40
radius[4] = 50
# this part is not even used:
# ctheta = -1. + 2 * numpy.random.rand(N)
# stheta = numpy.sqrt(1.-ctheta**2)
# phi    =  2 * math.pi * numpy.random.rand(N)
# end
r = numpy.zeros((numPart, 3))
r[:, 0] = radius

# import matplotlib.pyplot as plt
# plt.plot(r[:,0],'.')
# plt.show()

# print('Mass = ', Mass)
# print('radius = ', radius)
# print('scaleLength = ',scaleLength)
#
v = numpy.zeros((numPart, 3))
# v[:,0] = hernquistcircvel(radius,Mass,scaleLength)
omega = v[:, 0] / radius
period = 2.0 * math.pi / omega
print("period = minimum = ", min(period), " maximum = ", max(period))
print("Circular velocity = minimum =", min(v[:, 0]), " maximum = ", max(v[:, 0]))

omegav = omega

v[:, 0] = -omegav * r[:, 1]
v[:, 1] = omegav * r[:, 0]

ds = grp1.create_dataset("Velocities", (numPart, 3), "f", data=v)

m = numpy.full((numPart,), mass, dtype="f")
ds = grp1.create_dataset("Masses", (numPart,), "f", data=m)

ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)
ds = grp1.create_dataset("ParticleIDs", (numPart,), "L", data=ids)

ds = grp1.create_dataset("Coordinates", (numPart, 3), "d", data=r)


file.close()
