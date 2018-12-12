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

# Generates a swift IC file for the 2D Sod Shock in a periodic box

unit_l_in_cgs = 3.086e18
unit_m_in_cgs = 2.94e55
unit_t_in_cgs = 3.086e18

# Parameters
gamma = 5./3.          # Gas adiabatic index
x_min = -1.
x_max = 1.
rho_L = 1.             # Density left state
rho_R = 0.140625       # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state
a_beg = 0.001
fileName = "sodShock.hdf5" 


#---------------------------------------------------
boxSize = (x_max - x_min)

glass_L = h5py.File("glassPlane_128.hdf5", "r")
glass_R = h5py.File("glassPlane_48.hdf5", "r")

pos_L = glass_L["/PartType0/Coordinates"][:,:] * 0.5
pos_R = glass_R["/PartType0/Coordinates"][:,:] * 0.5
h_L = glass_L["/PartType0/SmoothingLength"][:] * 0.5
h_R = glass_R["/PartType0/SmoothingLength"][:] * 0.5

# Merge things
aa = pos_L - array([0.5, 0., 0.])
pos_LL = append(pos_L, pos_L + array([0.5, 0., 0.]), axis=0)
pos_RR = append(pos_R, pos_R + array([0.5, 0., 0.]), axis=0)
pos = append(pos_LL - array([1.0, 0., 0.]), pos_RR, axis=0)
h_LL = append(h_L, h_L)
h_RR = append(h_R, h_R)
h = append(h_LL, h_RR)

numPart_L = size(h_LL)
numPart_R = size(h_RR)
numPart = size(h)

vol_L = 0.5
vol_R = 0.5

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)

for i in range(numPart):
    x = pos[i,0]

    if x < 0: #left
        u[i] = P_L / (rho_L * (gamma - 1.))
        m[i] = rho_L * vol_L / numPart_L
        v[i,0] = v_L
    else:     #right
        u[i] = P_R / (rho_R * (gamma - 1.))
        m[i] = rho_R * vol_R / numPart_R
        v[i,0] = v_R
        
# Shift particles
pos[:,0] -= x_min

u /= (a_beg**(3. * (gamma - 1.)))

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, 0.5, 1.0]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 2

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_l_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_m_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_t_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')


file.close()
