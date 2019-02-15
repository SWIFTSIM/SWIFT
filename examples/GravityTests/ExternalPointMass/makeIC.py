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

import h5py
import sys
import numpy
import math
import random

# Generates a random distriution of particles, for motion in an external potential centred at (0,0,0)

# physical constants in cgs
NEWTON_GRAVITY_CGS  = 6.67408e-8
SOLAR_MASS_IN_CGS   = 1.98848e33
PARSEC_IN_CGS       = 3.08567758e18

# choice of units
const_unit_length_in_cgs   =   (1000*PARSEC_IN_CGS)
const_unit_mass_in_cgs     =   (SOLAR_MASS_IN_CGS)
const_unit_velocity_in_cgs =   (1e5)

print("UnitMass_in_cgs:     ", const_unit_mass_in_cgs) 
print("UnitLength_in_cgs:   ", const_unit_length_in_cgs)
print("UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs)
print("UnitTime_in_cgs:     ", const_unit_length_in_cgs / const_unit_velocity_in_cgs)

# derived units
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
const_G                = ((NEWTON_GRAVITY_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print('---------------------')
print('G in internal units: ', const_G)


# Parameters
periodic   = 1            # 1 For periodic box
boxSize    = 100.         # 
max_radius = boxSize / 4. # maximum radius of particles
Mass       = 1e10         
print("Mass at the centre:  ", Mass)

numPart = int(sys.argv[1])  # Number of particles
mass    = 1.

fileName = "PointMass.hdf5" 



#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [0, numPart, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, numPart, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3


#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 3.0856776e21
grp.attrs["Unit mass in cgs (U_M)"] = 1.98848e33
grp.attrs["Unit time in cgs (U_t)"] = 3.0856776e16
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp1 = file.create_group("/PartType1")

#generate particle positions
radius = max_radius * (numpy.random.rand(numPart))**(1./3.)
print('---------------------')
print('Radius: minimum = ',min(radius))
print('Radius: maximum = ',max(radius))
radius = numpy.sort(radius)
r      = numpy.zeros((numPart, 3))
r[:,0] = radius

#generate particle velocities
speed  = numpy.sqrt(const_G * Mass / radius)
omega  = speed / radius
period = 2.*math.pi/omega
print('---------------------')
print('Period: minimum = ',min(period))
print('Period: maximum = ',max(period))

v      = numpy.zeros((numPart, 3))
v[:,0] = -omega * r[:,1]
v[:,1] =  omega * r[:,0]

ds = grp1.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
v = numpy.zeros(1)

m = numpy.full((numPart, ), mass)
ds = grp1.create_dataset('Masses', (numPart,), 'f')
ds[()] = m
m = numpy.zeros(1)

ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)
ds = grp1.create_dataset('ParticleIDs', (numPart, ), 'L')
ds[()] = ids

ds = grp1.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = r

file.close()
