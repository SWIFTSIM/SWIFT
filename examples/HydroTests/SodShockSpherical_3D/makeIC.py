###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a swift IC file for the 3D Sod Shock in a periodic box

# Parameters
gamma = 5./3.          # Gas adiabatic index
x_min = -1.
x_max = 1.
rho_L = 1.             # Density left state
rho_R = 0.125          # Density right state
v_L = 0.               # Velocity left state
v_R = 0.               # Velocity right state
P_L = 1.               # Pressure left state
P_R = 0.1              # Pressure right state
fileName = "sodShock.hdf5" 


#---------------------------------------------------
boxSize = (x_max - x_min)

glass_L = h5py.File("glassCube_64.hdf5", "r")
glass_R = h5py.File("glassCube_32.hdf5", "r")

pos_L = glass_L["/PartType0/Coordinates"][:,:]
pos_R = glass_R["/PartType0/Coordinates"][:,:]
h_L = glass_L["/PartType0/SmoothingLength"][:]
h_R = glass_R["/PartType0/SmoothingLength"][:]

radius_L = sqrt((pos_L[:,0] - 0.5)**2 + (pos_L[:,1] - 0.5)**2 + \
                (pos_L[:,2] - 0.5)**2)
index_L = radius_L < 0.25
pos_LL = pos_L[index_L]
h_LL = h_L[index_L]

radius_R = sqrt((pos_R[:,0] - 0.5)**2 + (pos_R[:,1] - 0.5)**2 + \
                (pos_R[:,2] - 0.5)**2)
index_R = radius_R > 0.25
pos_RR = pos_R[index_R]
h_RR = h_R[index_R]

# Merge things
pos = append(pos_LL, pos_RR, axis=0)
h = append(h_LL, h_RR)

numPart_L = size(h_LL)
numPart_R = size(h_RR)
numPart = size(h)

vol_L = 4. * pi / 3. * 0.25**3
vol_R = 1. - 4. * pi / 3. * 0.25**3

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = zeros(numPart)
u = zeros(numPart)

for i in range(numPart):
    x = sqrt((pos[i,0]-0.5)**2+(pos[i,1]-0.5)**2+(pos[i,2]-0.5)**2)

    if x < 0.25: #left
        u[i] = P_L / (rho_L * (gamma - 1.))
        m[i] = rho_L * vol_L / numPart_L
        v[i,0] = v_L
    else:     #right
        u[i] = P_R / (rho_R * (gamma - 1.))
        m[i] = rho_R * vol_R / numPart_R
        v[i,0] = v_R
        
#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 1.]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
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
