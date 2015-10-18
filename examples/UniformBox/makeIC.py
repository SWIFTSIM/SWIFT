###############################################################################
 # This file is part of SWIFT.
 # Coypright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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
import sys
from numpy import *

# Generates a swift IC file containing a cartesian distribution of particles
# at a constant density and pressure in a cubic box

# Parameters
periodic= 1           # 1 For periodic box
boxSize = 1.
L = int(sys.argv[1])  # Number of particles along one axis
rho = 2.              # Density
P = 1.                # Pressure
gamma = 5./3.         # Gas adiabatic index
fileName = "uniformBox.hdf5" 


#---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.)*rho)

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
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]


#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

#Particle group
grp = file.create_group("/PartType0")

v  = zeros((numPart, 3))
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
v = zeros(1)

m = full((numPart, 1), mass)
ds = grp.create_dataset('Masses', (numPart,1), 'f')
ds[()] = m
m = zeros(1)

h = full((numPart, 1), 2.251 * boxSize / L)
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h
h = zeros(1)

u = full((numPart, 1), internalEnergy)
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u
u = zeros(1)


ids = linspace(0, numPart, numPart, endpoint=False, dtype='L').reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids
x      = ids % L;
y      = ((ids - x) / L) % L;
z      = (ids - x - L * y) / L**2;
coords = zeros((numPart, 3))
coords[:,0] = z[:,0] * boxSize / L + boxSize / (2*L)
coords[:,1] = y[:,0] * boxSize / L + boxSize / (2*L)
coords[:,2] = x[:,0] * boxSize / L + boxSize / (2*L)
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = coords

file.close()
