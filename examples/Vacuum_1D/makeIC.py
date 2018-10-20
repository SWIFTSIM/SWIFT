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

# Generates an overdensity within a vacuum to test the vacuum handling
# capabilities of the code

import numpy as np
import h5py

fileName = "vacuum.hdf5"
numPart = 100
boxSize = 1.
gamma = 5. / 3.

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart)
h = np.zeros(numPart)
u = np.zeros(numPart)
ids = np.arange(numPart, dtype = 'L')
rho = np.zeros(numPart)

# first set the positions, as we try to do a reasonable volume estimate to
# set the masses
for i in range(numPart):
  # we only generate particles in the range [0.25, 0.75]
  coords[i,0] = 0.25 + 0.5 * (i + 0.5) / numPart
  rho[i] = 1.
  P = 1.
  u[i] = P / (gamma - 1.) / rho[i]
  m[i] = rho[i] * 0.5 / numPart
  # reasonable smoothing length estimate
  h[i] = 1. / numPart

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
