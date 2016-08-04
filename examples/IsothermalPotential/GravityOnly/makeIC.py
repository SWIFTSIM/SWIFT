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

# Generates N particles in a spherical distribution centred on [0,0,0], to be moved in an isothermal potential
# usage: python makeIC.py 1000 0 1: generate 1000 particles on circular orbits
#        python makeIC.py 1000 1 1: generate 1000 particles with Lz/L uniform in [0,1]
# all particles move in the xy plane, and start at y=0

# physical constants in cgs
NEWTON_GRAVITY_CGS  = 6.672e-8
SOLAR_MASS_IN_CGS   = 1.9885e33
PARSEC_IN_CGS       = 3.0856776e18
PROTON_MASS_IN_CGS  = 1.6726231e24
YEAR_IN_CGS         = 3.154e+7

# choice of units
const_unit_length_in_cgs   =   (1000*PARSEC_IN_CGS)
const_unit_mass_in_cgs     =   (SOLAR_MASS_IN_CGS)
const_unit_velocity_in_cgs =   (1e5)

print "UnitMass_in_cgs:     ", const_unit_mass_in_cgs 
print "UnitLength_in_cgs:   ", const_unit_length_in_cgs
print "UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs


# rotation speed of isothermal potential [km/s]
vrot_kms = 200.

# derived units
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
const_G                = ((NEWTON_GRAVITY_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print 'G=', const_G
vrot   = vrot_kms * 1e5 / const_unit_velocity_in_cgs


# Parameters
periodic= 1            # 1 For periodic box
boxSize = 400.         #  [kpc]
Radius  = 100.          # maximum radius of particles [kpc]
G       = const_G 

N       = int(sys.argv[1])  # Number of particles
icirc   = int(sys.argv[2])  # if = 0, all particles are on circular orbits, if = 1, Lz/Lcirc uniform in ]0,1[
three_d = int(sys.argv[3])  # if = 0, particles distributed in 3 dimensions, if 1, all along one straight line
L       = N**(1./3.)

# these are not used but necessary for I/O
rho = 2.              # Density
P = 1.                # Pressure
gamma = 5./3.         # Gas adiabatic index
fileName = "Isothermal.hdf5" 


#---------------------------------------------------
numPart        = N
mass           = 1
internalEnergy = P / ((gamma - 1.)*rho)

#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs 
grp.attrs["Unit time in cgs (U_t)"] = const_unit_length_in_cgs / const_unit_velocity_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

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


#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# set seed for random number
numpy.random.seed(1234)

#Particle group
#grp0 = file.create_group("/PartType0")
grp1 = file.create_group("/PartType1")
#generate particle positions
radius = Radius * (numpy.random.rand(N))**(1./3.) 
ctheta = -1. + 2 * numpy.random.rand(N)
stheta = numpy.sqrt(1.-ctheta**2)
phi    =  2 * math.pi * numpy.random.rand(N)
r      = numpy.zeros((numPart, 3))

if (three_d == 0):
    r[:,0] = radius * stheta * numpy.cos(phi)
    r[:,1] = radius * stheta * numpy.sin(phi)
    r[:,2] = radius * ctheta
else:
    r[:,0] = radius
#

speed  = vrot
v      = numpy.zeros((numPart, 3))

if(three_d == 0):
    #make N random angular momentum vectors
    omega = numpy.random.rand(N,3)
    #cross with r to get velocities perpendicular to displacement vector
    v = numpy.cross(omega,r)
    #rescale v
    mag_v = numpy.zeros(N)
    for i in range(N):
        mag_v[i] = numpy.sqrt(v[i,0]**2 + v[i,1]**2 + v[i,2]**2)
        v[i,:] /= mag_v[i]
    v *= speed
    if (icirc != 0):
        for i in range(N):
            v[i,:] *= numpy.random.rand()
    
else:
    omega  = speed / radius
    period = 2.*math.pi/omega
    print 'period = minimum = ',min(period), ' maximum = ',max(period)

    omegav = omega
    if (icirc != 0):
        omegav = omega * numpy.random.rand(N)
        
    v[:,0] = -omegav * r[:,1]
    v[:,1] =  omegav * r[:,0]

ds = grp1.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
v = numpy.zeros(1)

m = numpy.full((numPart, ), mass)
ds = grp1.create_dataset('Masses', (numPart,), 'f')
ds[()] = m
m = numpy.zeros(1)

h = numpy.full((numPart, ), 1.1255 * boxSize / L)
ds = grp1.create_dataset('SmoothingLength', (numPart,), 'f')
ds[()] = h
h = numpy.zeros(1)

u = numpy.full((numPart, ), internalEnergy)
ds = grp1.create_dataset('InternalEnergy', (numPart,), 'f')
ds[()] = u
u = numpy.zeros(1)


ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False, dtype='L')
ds = grp1.create_dataset('ParticleIDs', (numPart, ), 'L')
ds[()] = ids

ds = grp1.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = r


file.close()
