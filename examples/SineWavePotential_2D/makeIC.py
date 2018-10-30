###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Generates a distorted 1D grid with a density profile that balances out the
# external sine wave potential if run with an isothermal equation of state.

import numpy as np
import h5py

# constant thermal energy
# the definition below assumes the same thermal energy is defined in const.h,
# and also that the code was configured with an adiabatic index of 5./3.
uconst = 20.2615290634
cs2 = 2.*uconst/3.
A = 10.

fileName = "sineWavePotential.hdf5"
numPart_1D = 50
boxSize = [1., 1.]
numPart = numPart_1D**2

coords = np.zeros((numPart, 3))
v = np.zeros((numPart, 3))
m = np.zeros(numPart) + 1.
h = np.zeros(numPart) + 2./numPart
u = np.zeros(numPart) + uconst
ids = np.arange(numPart, dtype = 'L')
rho = np.zeros(numPart)

# first set the positions, as we try to do a reasonable volume estimate to
# set the masses
for i in range(numPart_1D):
  for j in range(numPart_1D):
    coords[numPart_1D*i+j,0] = (i+0.5)/numPart_1D
    coords[numPart_1D*i+j,1] = (j+0.5)/numPart_1D

V = 1./numPart
for i in range(numPart):
  rho[i] = 1000.*np.exp(-0.5*A/np.pi/cs2*np.cos(2.*np.pi*coords[i,0]))
  m[i] = rho[i]*V

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
grp.create_dataset('Coordinates', data=coords, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')
grp.create_dataset('Density', data=rho, dtype='f')

file.close()
