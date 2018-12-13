###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 #               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
from numpy import *

# Generates a swift IC file for the 1D Sod Shock in a periodic box

unit_l_in_cgs = 3.086e18
unit_m_in_cgs = 2.94e55
unit_t_in_cgs = 3.086e18

# Parameters
gamma = 5./3.          # Gas adiabatic index
numPart_L = 800        # Number of particles in the left state
x_min = -1.
x_max = 1.
rho_L = 1.             # Density left state
rho_R = 0.125          # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state
a_beg = 0.001
fileName = "sodShock.hdf5" 


#---------------------------------------------------

# Find how many particles we actually have
boxSize = x_max - x_min
numPart_R = int(numPart_L * (rho_R / rho_L))
numPart = numPart_L + numPart_R

# Now get the distances
delta_L = (boxSize/2)  / numPart_L
delta_R = (boxSize/2)  / numPart_R
offset_L = delta_L / 2
offset_R = delta_R / 2

# Build the arrays
coords = zeros((numPart, 3))
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
h = zeros(numPart)
u = zeros(numPart)

# Set the particles on the left
for i in range(numPart_L):
    coords[i,0] = x_min + offset_L + i * delta_L
    u[i] = P_L / (rho_L * (gamma - 1.))
    h[i] = 1.2348 * delta_L
    m[i] = boxSize * rho_L / (2. * numPart_L)
    v[i,0] = v_L
    
# Set the particles on the right
for j in range(numPart_R):
    i = numPart_L + j
    coords[i,0] = offset_R + j * delta_R
    u[i] = P_R / (rho_R * (gamma - 1.))
    h[i] = 1.2348 * delta_R
    m[i] = boxSize * rho_R / (2. * numPart_R)
    v[i,0] = v_R

# Shift particles
coords[:,0] -= x_min

u /= (a_beg**(3. * (gamma - 1.)))
    
#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 1

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_l_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_m_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_t_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=coords, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')


file.close()


