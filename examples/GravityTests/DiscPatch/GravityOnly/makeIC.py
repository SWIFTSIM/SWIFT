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

# Generates N particles in a box of [0:BoxSize,0:BoxSize,-2scale_height:2scale_height]
# see Creasey, Theuns & Bower, 2013, for the equations:
# disc parameters are: surface density sigma
#                      scale height b
# density: rho(z) = (sigma/2b) sech^2(z/b)
# isothermal velocity dispersion = <v_z^2? = b pi G sigma
# grad potential  = 2 pi G sigma tanh(z/b)
# potential       = ln(cosh(z/b)) + const
# Dynamical time  = sqrt(b / (G sigma))
# to obtain the 1/ch^2(z/b) profile from a uniform profile (a glass, say, or a uniform random variable), note that, when integrating in z
# \int 0^z dz/ch^2(z) = tanh(z)-tanh(0) = \int_0^x dx = x (where the last integral refers to a uniform density distribution), so that z = atanh(x)
# usage: python makeIC.py 1000 

# physical constants in cgs
NEWTON_GRAVITY_CGS  = 6.67408e-8
SOLAR_MASS_IN_CGS   = 1.98848e33
PARSEC_IN_CGS       = 3.08567758e18
YEAR_IN_CGS         = 3.15569252e7

# choice of units
const_unit_length_in_cgs   =   (PARSEC_IN_CGS)
const_unit_mass_in_cgs     =   (SOLAR_MASS_IN_CGS)
const_unit_velocity_in_cgs =   (1e5)

print "UnitMass_in_cgs:     ", const_unit_mass_in_cgs 
print "UnitLength_in_cgs:   ", const_unit_length_in_cgs
print "UnitVelocity_in_cgs: ", const_unit_velocity_in_cgs


# parameters of potential
surface_density = 10.
scale_height    = 100.

# derived units
const_unit_time_in_cgs = (const_unit_length_in_cgs / const_unit_velocity_in_cgs)
const_G                = ((NEWTON_GRAVITY_CGS*const_unit_mass_in_cgs*const_unit_time_in_cgs*const_unit_time_in_cgs/(const_unit_length_in_cgs*const_unit_length_in_cgs*const_unit_length_in_cgs)))
print 'G=', const_G
v_disp                 = numpy.sqrt(scale_height * math.pi * const_G * surface_density)
t_dyn                  = numpy.sqrt(scale_height / (const_G * surface_density))
print 'dynamical time = ',t_dyn
print ' velocity dispersion = ',v_disp

# Parameters
periodic= 1             # 1 For periodic box
boxSize = 600.          #  
Radius  = 100.          # maximum radius of particles [kpc]
G       = const_G 

N       = int(sys.argv[1])  # Number of particles

# these are not used but necessary for I/O
rho = 2.              # Density
P = 1.                # Pressure
gamma = 5./3.         # Gas adiabatic index
fileName = "Disc-Patch.hdf5" 


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
grp.attrs["Dimension"] = 3

# set seed for random number
numpy.random.seed(1234)

#Particle group
#grp0 = file.create_group("/PartType0")
grp1 = file.create_group("/PartType1")

#generate particle positions
r      = numpy.zeros((numPart, 3))
r[:,0] = numpy.random.rand(N) * boxSize
r[:,1] = numpy.random.rand(N) * boxSize
z      = scale_height * numpy.arctanh(numpy.random.rand(2*N))
gd     = z < boxSize / 2
r[:,2] = z[gd][0:N]
random = numpy.random.rand(N) > 0.5
r[random,2] *= -1
r[:,2] += 0.5 * boxSize

#generate particle velocities
v      = numpy.zeros((numPart, 3))
v      = numpy.zeros(1)
#v[:,2] = 


ds = grp1.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v


m = numpy.ones((numPart, ), dtype=numpy.float32) * mass
ds = grp1.create_dataset('Masses', (numPart,), 'f')
ds[()] = m
m = numpy.zeros(1)


ids = 1 + numpy.linspace(0, numPart, numPart, endpoint=False)
ds = grp1.create_dataset('ParticleIDs', (numPart, ), 'L')
ds[()] = ids

ds = grp1.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = r


file.close()
