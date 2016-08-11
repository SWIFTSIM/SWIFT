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
L = 101           # Number of particles along one axis
rho = 1.          # Density
P = 1.e-5         # Pressure
E0= 1.e2          # Energy of the explosion
pert = 0.1
gamma = 5./3.     # Gas adiabatic index
eta = 1.2349      # 48 ngbs with cubic spline kernel
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
m      = zeros(numPart)
h      = zeros(numPart)
u      = zeros(numPart)
ids    = zeros(numPart, dtype='L')

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
            h[index] = eta * boxSize / L
            u[index] = internalEnergy
            ids[index] = index + 1
            if sqrt((x - boxSize/2.)**2 + (y - boxSize/2.)**2 + (z - boxSize/2.)**2) < 2.01 * boxSize/L:
                u[index] = u[index] + E0 / (33. * mass)
                print "Particle " , index , " set to detonate."
            coords[index,0] = x + random.random() * pert * boxSize/(2.*L)
            coords[index,1] = y + random.random() * pert * boxSize/(2.*L)
            coords[index,2] = z + random.random() * pert * boxSize/(2.*L)


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

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
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
