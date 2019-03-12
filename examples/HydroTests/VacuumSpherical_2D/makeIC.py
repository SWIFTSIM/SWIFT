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
###############################################################################

import numpy as np
import h5py

# Generates an overdensity within a vacuum to test the vacuum resolving
# capabilities of the code

# Parameters
gamma = 5. / 3.    # Gas adiabatic index

fileName = "vacuum.hdf5" 

#---------------------------------------------------
glass = h5py.File("glassPlane_128.hdf5", "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:,:]
h = glass["/PartType0/SmoothingLength"][:] * 0.3

# Make 4 copies of the glass to have more particles
pos *= 0.5
h *= 0.5
pos = np.append(pos, pos + np.array([0.5, 0., 0.]), axis = 0)
h = np.append(h, h)
pos = np.append(pos, pos + np.array([0., 0.5, 0.]), axis = 0)
h = np.append(h, h)

radius = np.sqrt((pos[:,0] - 0.5)**2 + (pos[:,1] - 0.5)**2)
index = radius < 0.25
pos = pos[index]
h = h[index]

numPart = len(h)
vol = np.pi * 0.25**2

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)

m[:] = 1. * vol / numPart
u[:] = 1. / (1. * (gamma - 1.))

#--------------------------------------------------

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
grp.attrs["Dimension"] = 2

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
