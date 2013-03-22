###############################################################################
 # This file is part of SWIFT.
 # Coypright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
import random
from numpy import *

# Generates a swift IC file for the Sod Shock in a periodic box

# Parameters
periodic= 1      # 1 For periodic box
boxSize = 1.
L = 100           # Number of particles along one axis
P1 = 1.           # Pressure left state
P2 = 0.1795       # Pressure right state
gamma = 5./3.     # Gas adiabatic index
fileName = "sodShock.hdf5" 


#---------------------------------------------------

#Read in high density glass
glass1 = h5py.File("../Glass/glass_200000.hdf5")
pos1 = glass1["/PartType1/Coordinates"][:,:]
pos1 = pos1 / 20. # Particles are in [0:0.5, 0:0.5, 0:0.5]

toRemove = array(0)
for i in range(size(pos1)/3-1):
    if pos1[i,0] == pos1[i+1,0]:
        if pos1[i,1] == pos1[i+1,1]:
            if pos1[i,2] == pos1[i+1,2]:
                toRemove = append(toRemove, i)

pos1 = delete(pos1, toRemove, 0)


#Read in high density glass
glass2 = h5py.File("../Glass/glass_50000.hdf5")
pos2 = glass2["/PartType1/Coordinates"][:,:]
pos2 = pos2 / 20. # Particles are in [0:0.5, 0:0.5, 0:0.5]

toRemove = array(0)
for i in range(size(pos2)/3-1):
    if pos2[i,0] == pos2[i+1,0]:
        if pos2[i,1] == pos2[i+1,1]:
            if pos2[i,2] == pos2[i+1,2]:
                toRemove = append(toRemove, i)
pos2 = delete(pos2, toRemove, 0)


#Generate high density region
rho1 = size(pos1) / size(pos2) 
coord1 = append(pos1, pos1 + [0, 0.5, 0], 0)
coord1 = append(coord1, pos1 + [0, 0, 0.5], 0)
coord1 = append(coord1, pos1 + [0, 0.5, 0.5], 0)
N1 = size(coord1)/3
v1 = zeros((N1, 3))
h1 = ones((N1, 1)) * 2.251 * 0.5 / (size(pos1)/3)**(1./3.)
u1 = ones((N1, 1)) * P1 / ((gamma - 1.) * rho1)
m1 = ones((N1, 1)) * 0.5 * rho1 / N1

#Generate low density region
rho2 = 1.
coord2 = append(pos2, pos2 + [0, 0.5, 0], 0)
coord2 = append(coord2, pos2 + [0, 0, 0.5], 0)
coord2 = append(coord2, pos2 + [0, 0.5, 0.5], 0)
N2 = size(coord2)/3
v2 = zeros((N2, 3))
h2 = ones((N2, 1)) * 2.251 * 0.5 / (size(pos2)/3)**(1./3.)
u2 = ones((N2, 1)) * P2 / ((gamma - 1.) * rho2)
m2 = ones((N2, 1)) * 0.5 * rho2 / N2

#Merge arrays
numPart = N1 + N2
coords = append(coord1, coord2+[0.5, 0., 0.], 0)
v = append(v1, v2,0)
h = append(h1, h2,0)
u = append(u1, u2,0)
m = append(m1, m2,0)
ids = zeros((numPart, 1), dtype='L')
for i in range(numPart):
    ids[i] = i

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
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

#Particle group
grp = file.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = coords
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
ds = grp.create_dataset('Masses', (numPart,1), 'f')
ds[()] = m
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u
ds = grp.create_dataset('ParticleIDs', (numPart,1), 'L')
ds[()] = ids[:]

file.close()


