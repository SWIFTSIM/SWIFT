###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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
boxSize = 10.
L = 128           # Number of particles boxes along one axis
rho = 1.          # Density
P = 1.e-5         # Pressure
E0= 1.e2          # Energy of the explosion
pert = 0.025
gamma = 5./3.     # Gas adiabatic index
fileName = "perturbedBox.hdf5" 


#---------------------------------------------------
numPart = 4*(L**3) 
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.)*rho)
off = array( [ [ 0.0 , 0.0 , 0.0 ] , [ 0.0 , 0.5 , 0.5 ] , [ 0.5 , 0.0 , 0.5 ] , [ 0.5 , 0.5 , 0.0 ] ] );
hbox = boxSize / L

# if L%2 == 0:
#     print "Number of particles along each dimension must be odd."
#     exit()

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
            x = (i + 0.25) * hbox
            y = (j + 0.25) * hbox
            z = (k + 0.25) * hbox
            for ell in range(4): 
                index = 4*(i*L*L + j*L + k) + ell
                coords[index,0] = x + off[ell,0] * hbox
                coords[index,1] = y + off[ell,1] * hbox
                coords[index,2] = z + off[ell,2] * hbox
                v[index,0] = 0.
                v[index,1] = 0.
                v[index,2] = 0.
                m[index] = mass
                h[index] = 2.251 / 4 * hbox
                u[index] = internalEnergy
                ids[index] = index
                coords[index,0] += random.random() * pert * hbox
                coords[index,1] += random.random() * pert * hbox
                coords[index,2] += random.random() * pert * hbox


#--------------------------------------------------

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
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids

file.close()
