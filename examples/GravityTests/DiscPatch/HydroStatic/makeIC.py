###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 John A. Regan (john.a.regan@durham.ac.uk)
#                    Tom Theuns (tom.theuns@durham.ac.uk)
#               2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#                    Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
import numpy as np
import math
import random

# Generates a disc-patch in hydrostatic equilibrium
#
# See Creasey, Theuns & Bower, 2013, MNRAS, Volume 429, Issue 3, p.1922-1948
#
#
# Disc parameters are: surface density  -- sigma
#                      scale height -- b
#                      gas adiabatic index -- gamma
#
# Problem parameters are: Ratio height/width of the box -- z_factor
#                         Size of the patch -- side_length

# Parameters of the gas disc
surface_density = 10.
scale_height    = 100.
gas_gamma       = 5./3.

# Parameters of the problem
x_factor        = 2
side_length     = 400.

# File
fileName = "Disc-Patch.hdf5"

####################################################################

# physical constants in cgs
NEWTON_GRAVITY_CGS  = 6.67408e-8
SOLAR_MASS_IN_CGS   = 1.98848e33
PARSEC_IN_CGS       = 3.08567758e18
PROTON_MASS_IN_CGS  = 1.672621898e-24
BOLTZMANN_IN_CGS    = 1.38064852e-16
YEAR_IN_CGS         = 3.15569252e7

# choice of units
unit_length_in_cgs   =   (PARSEC_IN_CGS)
unit_mass_in_cgs     =   (SOLAR_MASS_IN_CGS)
unit_velocity_in_cgs =   (1e5)
unit_time_in_cgs     =   unit_length_in_cgs / unit_velocity_in_cgs

print "UnitMass_in_cgs:     %.5e"%unit_mass_in_cgs
print "UnitLength_in_cgs:   %.5e"%unit_length_in_cgs
print "UnitVelocity_in_cgs: %.5e"%unit_velocity_in_cgs
print "UnitTime_in_cgs:     %.5e"%unit_time_in_cgs
print ""

# Derived units
const_G  = NEWTON_GRAVITY_CGS * unit_mass_in_cgs * unit_time_in_cgs**2 * \
           unit_length_in_cgs**-3
const_mp = PROTON_MASS_IN_CGS * unit_mass_in_cgs**-1
const_kb = BOLTZMANN_IN_CGS * unit_mass_in_cgs**-1 * unit_length_in_cgs**-2 * \
           unit_time_in_cgs**2

print "--- Some constants [internal units] ---"
print "G_Newton:    %.5e"%const_G
print "m_proton:    %.5e"%const_mp
print "k_boltzmann: %.5e"%const_kb
print ""

# derived quantities
temp       = math.pi * const_G * surface_density * scale_height * const_mp / \
             const_kb
u_therm    = const_kb * temp / ((gas_gamma-1) * const_mp)
v_disp     = math.sqrt(2 * u_therm)
soundspeed = math.sqrt(u_therm / (gas_gamma * (gas_gamma-1.)))
t_dyn      = math.sqrt(scale_height / (const_G * surface_density))
t_cross    = scale_height / soundspeed

print "--- Properties of the gas [internal units] ---"
print "Gas temperature:     %.5e"%temp
print "Gas thermal_energy:  %.5e"%u_therm
print "Dynamical time:      %.5e"%t_dyn
print "Sound crossing time: %.5e"%t_cross
print "Gas sound speed:     %.5e"%soundspeed
print "Gas 3D vel_disp:     %.5e"%v_disp
print ""

# Problem properties
boxSize_x = side_length
boxSize_y = boxSize_x
boxSize_z = boxSize_x
boxSize_x *= x_factor
volume = boxSize_x * boxSize_y * boxSize_z
M_tot = boxSize_y * boxSize_z * surface_density * \
        math.tanh(boxSize_x / (2. * scale_height))
density = M_tot / volume
entropy = (gas_gamma - 1.) * u_therm / density**(gas_gamma - 1.)

print "--- Problem properties [internal units] ---"
print "Box:        [%.1f, %.1f, %.1f]"%(boxSize_x, boxSize_y, boxSize_z)
print "Volume:     %.5e"%volume
print "Total mass: %.5e"%M_tot
print "Density:    %.5e"%density
print "Entropy:    %.5e"%entropy
print ""

####################################################################

# Read glass pre-ICs
infile  = h5py.File('glassCube_32.hdf5', "r")
one_glass_pos = infile["/PartType0/Coordinates"][:,:]
one_glass_h   = infile["/PartType0/SmoothingLength"][:]

# Rescale to the problem size
one_glass_pos *= side_length
one_glass_h   *= side_length

# Now create enough copies to fill the volume in x
pos = np.copy(one_glass_pos)
h = np.copy(one_glass_h)
for i in range(1, x_factor):
    one_glass_pos[:,0] += side_length
    pos = np.append(pos, one_glass_pos, axis=0)
    h   = np.append(h, one_glass_h, axis=0)

# Compute further properties of ICs
numPart = np.size(h)
mass = M_tot / numPart

print "--- Particle properties [internal units] ---"
print "Number part.: ", numPart
print "Part. mass:   %.5e"%mass
print ""

# Create additional arrays
u    = np.ones(numPart) * u_therm
mass = np.ones(numPart) * mass
vel  = np.zeros((numPart, 3))
ids  = 1 + np.linspace(0, numPart, numPart, endpoint=False)

####################################################################
# Create and write output file

#File
file = h5py.File(fileName, 'w')

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_time_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize_x, boxSize_y, boxSize_z]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# write gas particles
grp0   = file.create_group("/PartType0")

ds = grp0.create_dataset('Coordinates', (numPart, 3), 'f', data=pos)
ds = grp0.create_dataset('Velocities', (numPart, 3), 'f')
ds = grp0.create_dataset('Masses', (numPart,), 'f', data=mass)
ds = grp0.create_dataset('SmoothingLength', (numPart,), 'f', data=h)
ds = grp0.create_dataset('InternalEnergy', (numPart,), 'f', data=u)
ds = grp0.create_dataset('ParticleIDs', (numPart, ), 'L', data=ids)

####################################################################

print "--- Runtime parameters (YAML file): ---"
print "DiscPatchPotential:surface_density:    ", surface_density
print "DiscPatchPotential:scale_height:       ", scale_height
print "DiscPatchPotential:x_disc:             ", 0.5 * boxSize_x
print "EoS:isothermal_internal_energy: %ef"%u_therm
print ""
