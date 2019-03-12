###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durhama.ac.uk)
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
import numpy as np

# Generates a SWIFT IC file with a constant density and pressure

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1  # 1 kiloparsec
rho = 3.2e3  # Density in code units (3.2e6 is 0.1 hydrogen atoms per cm^3)
T = 4e3  # Initial Temperature
gamma = 5./3.  # Gas adiabatic index
fileName = "coolingBox.hdf5"
# ---------------------------------------------------

# defines some constants
# need to be changed in plotTemperature.py too
h_frac = 0.76
mu = 4. / (1. + 3. * h_frac)

m_h_cgs = 1.67e-24
k_b_cgs = 1.38e-16

# defines units
unit_length = 3.0857e21  # kpc
unit_mass = 2.0e33  # solar mass
unit_time = 3.0857e16  # ~ Gyr

# Read id, position and h from glass
glass = h5py.File("glassCube_32.hdf5", "r")
ids = glass["/PartType0/ParticleIDs"][:]
pos = glass["/PartType0/Coordinates"][:, :] * boxSize
h = glass["/PartType0/SmoothingLength"][:] * boxSize

# Compute basic properties
numPart = np.size(pos) // 3
mass = boxSize**3 * rho / numPart
internalEnergy = k_b_cgs * T * mu / ((gamma - 1.) * m_h_cgs)
internalEnergy *= (unit_time / unit_length)**2

# File
f = h5py.File(fileName, 'w')

# Header
grp = f.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0

# Runtime parameters
grp = f.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

# Units
grp = f.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass
grp.attrs["Unit time in cgs (U_t)"] = unit_time
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = f.create_group("/PartType0")

v = np.zeros((numPart, 3))
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = v

m = np.full((numPart, 1), mass)
ds = grp.create_dataset('Masses', (numPart, 1), 'f')
ds[()] = m

h = np.reshape(h, (numPart, 1))
ds = grp.create_dataset('SmoothingLength', (numPart, 1), 'f')
ds[()] = h

u = np.full((numPart, 1), internalEnergy)
ds = grp.create_dataset('InternalEnergy', (numPart, 1), 'f')
ds[()] = u

ids = np.reshape(ids, (numPart, 1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids

ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = pos

f.close()

print("Initial condition generated")
