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
import numpy as np

# Generates a swift IC file containing a cartesian distribution of DM particles
# with a density of 1

# Parameters
periodic= 1           # 1 For periodic box
boxSize = 1.
rho = 1.
L = int(sys.argv[1])  # Number of particles along one axis
fileName = "uniform_DM_box.hdf5"

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

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.


#Particle group
grp = file.create_group("/PartType1")

v  = np.zeros((numPart, 3))
ds = grp.create_dataset('Velocities', (numPart, 3), 'f', data=v)

m = np.full((numPart, 1), mass)
ds = grp.create_dataset('Masses', (numPart,1), 'f', data=m)

ids = np.linspace(0, numPart, numPart, endpoint=False).reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids + 1
coords = np.random.rand(numPart, 3) * boxSize
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd', data=coords)

file.close()
