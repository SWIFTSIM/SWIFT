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

# Generates a swift IC file for the Sedov blast test in a periodic cubic box

# Parameters
periodic= 1      # 1 For periodic box
boxSize = 1.
L = 101           # Number of particles along one axis
rho = 1.          # Density
P = 1.e-5         # Pressure
E0= 1.e5          # Energy of the explosion
pert = 0.1
gamma = 5./3.     # Gas adiabatic index
fileName = "sedov.hdf5" 


#---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.)*rho)

if L%2 == 0:
    print "Number of particles along each dimension must be odd."
    exit()

#Generate particles
coords = zeros((numPart, 3))
v      = zeros((numPart, 3))
m      = zeros((numPart, 1))
h      = zeros((numPart, 1))
u      = zeros((numPart, 1))
ids    = zeros((numPart, 1), dtype='L')

for i in range(L):
    for j in range(L):
        for k in range(L):
            index = i*L*L + j*L + k
            x = i * boxSize / L + boxSize / (2*L)
            y = j * boxSize / L + boxSize / (2*L)
            z = k * boxSize / L + boxSize / (2*L)
            coords[index,0] = x
            coords[index,1] = y
            coords[index,2] = z
            v[index,0] = 0.
            v[index,1] = 0.
            v[index,2] = 0.
            m[index] = mass
            h[index] = 2.251 * boxSize / L
            u[index] = internalEnergy
            ids[index] = index
            if sqrt((x - boxSize/2.)**2 + (y - boxSize/2.)**2 + (z - boxSize/2.)**2) < 2.01 * boxSize/L:
                u[index] = u[index] + E0 / 33.
            coords[index,0] = x + random.random() * pert * boxSize/(2.*L)
            coords[index,1] = y + random.random() * pert * boxSize/(2.*L)
            coords[index,2] = z + random.random() * pert * boxSize/(2.*L)


#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = numPart

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

#Particle group
grp = file.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = coords
ds = grp.create_dataset('Velocity', (numPart, 3), 'f')
ds[()] = v
ds = grp.create_dataset('Mass', (numPart,1), 'f')
ds[()] = m
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids

file.close()
