###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a swift IC file for the Smoothed Metallicity test in a periodic cubic box

# Parameters
gamma = 5./3.      # Gas adiabatic index
rho0 = 1.          # Background density
P0 = 1.e-6         # Background pressure
Nelem = 9          # Gear: 9, EAGLE: 9
low_metal = -6     # Low iron fraction
high_metal = -5    # high iron fraction
sigma_metal = 0.1  # relative standard deviation for the metallicities
fileName = "smoothed_metallicity.hdf5"

# shift all metals in order to obtain nicer plots
low_metal = [low_metal] * Nelem + np.linspace(0, 3, Nelem)
high_metal = [high_metal] * Nelem + np.linspace(0, 3, Nelem)

# ---------------------------------------------------
glass = h5py.File("glassCube_32.hdf5", "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

numPart = np.size(h)
vol = 1.

# Generate extra arrays
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
u = np.zeros(numPart)
r = np.zeros(numPart)
Z = np.zeros((numPart, Nelem))

r = np.sqrt((pos[:, 0] - 0.5)**2 + (pos[:, 1] - 0.5)**2 + (pos[:, 2] - 0.5)**2)
m[:] = rho0 * vol / numPart
u[:] = P0 / (rho0 * (gamma - 1))

# set metallicities
select = pos[:, 0] < 0.5
nber = sum(select)
Z[select, :] = low_metal * (1 + np.random.normal(loc=0., scale=sigma_metal,
                                                 size=(nber, Nelem)))
nber = numPart - nber
Z[np.logical_not(select), :] = high_metal * (1 + np.random.normal(
    loc=0., scale=sigma_metal, size=(nber, Nelem)))

# --------------------------------------------------

# File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 1.]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')
grp.create_dataset('ElementAbundance', data=Z, dtype='f')


file.close()
