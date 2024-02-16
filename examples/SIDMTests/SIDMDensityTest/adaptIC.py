###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Camila Correa (camila.correa@cea.fr)
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

sim = h5py.File("EAGLE_ICs_6.hdf5", "r")
pos = sim["/PartType0/Coordinates"][:, :]
mass = sim["/PartType0/Masses"][:]
h = sim["/PartType0/SmoothingLength"][:]
vel = sim["/PartType0/Velocities"][:, :]
e = sim["/PartType0/InternalEnergy"][:]

numPart = len(h)

#File
file = h5py.File('ICs_Hydro.hdf5', 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = 4.235625
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [mass[0], 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Set units
unit_length_in_cgs = 3.085678e24
unit_mass_in_cgs = 1.98848e43
unit_time_in_cgs = 3.085678e19
unit_current_in_cgs = 1
unit_temperature_in_cgs = 1
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_time_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = unit_current_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = unit_temperature_in_cgs


#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', (numPart, 3), 'd', data=pos)
grp.create_dataset('Velocities', (numPart, 3), 'f', data=vel)
grp.create_dataset('Masses', (numPart,1), 'f', data=mass)
grp.create_dataset('InternalEnergy', (numPart, 1), 'd',data=e)
grp.create_dataset('SmoothingLength', (numPart, 1), 'd',data=h)
ids = np.linspace(0, numPart, numPart, endpoint=False).reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart, 1), 'L')
ds[()] = ids + 1

file.close()


#####

#File
file = h5py.File('ICs_DM.hdf5', 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = 4.235625
grp.attrs["NumPart_Total"] =  [0, numPart, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [0, numPart, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, mass[0], 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Set units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_time_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = unit_current_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = unit_temperature_in_cgs


#Particle group
grp = file.create_group("/PartType1")
grp.create_dataset('Coordinates', (numPart, 3), 'd', data=pos)
grp.create_dataset('Velocities', (numPart, 3), 'f', data=vel)
grp.create_dataset('Masses', (numPart,1), 'f', data=mass)
grp.create_dataset('SmoothingLength', (numPart, 1), 'd',data=h)
grp.create_dataset('ParticleIDs', (numPart, 1), 'L', data=ids+1)
file.close()


