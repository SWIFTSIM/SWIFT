###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

# Generates a swift IC file containing a cartesian distribution of DM particles
# with a density of 1

# Parameters
periodic= 1           # 1 For periodic box
boxSize = 1.
rho = 1.
L = int(sys.argv[1])  # Number of particles along one axis
fileName = "uniformDMBox_%d.hdf5"%L 

#---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart

#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [0, numPart, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, numPart, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, mass, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

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
grp = file.create_group("/PartType1")

v  = zeros((numPart, 3))
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v
v = zeros(1)

m = full((numPart, 1), mass)
ds = grp.create_dataset('Masses', (numPart,1), 'f')
ds[()] = m
m = zeros(1)

ids = linspace(0, numPart, numPart, endpoint=False).reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids + 1
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
