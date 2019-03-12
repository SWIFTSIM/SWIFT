###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a regular 1D grid that contains the Woodward & Colella (1984)
# interacting blast waves test.
# Note that the original test requires reflective boundary conditions. Since we
# do not have those in SWIFT, we double the computational domain and mirror the
# setup in the second half.

import numpy as np
import h5py

fileName = "interactingBlastWaves.hdf5"
numPart = 800
boxSize = 2.

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart) + boxSize / numPart
h = np.zeros(numPart) + 2. * boxSize / numPart
u = np.zeros(numPart)
ids = np.arange(numPart, dtype = 'L')
rho = np.ones(numPart)

for i in range(numPart):
  coords[i,0] = (i + 0.5) * boxSize / numPart
  if coords[i,0] < 0.1 or coords[i,0] > 1.9:
    u[i] = 2500.
  elif coords[i,0] > 0.9 and coords[i,0] < 1.1:
    u[i] = 250.
  else:
    u[i] = 0.025

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
grp.create_dataset('Density', data=rho, dtype='f')

file.close()
