################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

import h5py
import numpy as np

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
gamma = 5./3.     # Gas adiabatic index
P0    = 2.5       # Pressure
rho0  = 1.        # Density
d     = 0.0317    # Thickness of the transition layer
B     = 0.0005    # Amplitude of the seed velocity

fileOutputName = "kelvinHelmholtzGrowthRate.hdf5"

#---------------------------------------------------

glass = h5py.File("glassPlane_128.hdf5", 'r')
pos = glass["/PartType0/Coordinates"][:, :]
h = glass["/PartType0/SmoothingLength"][:]

N = len(h)
vol = 1.

# Generate extra arrays
v = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.))

v[pos[:, 1] <= 0.25 - d, 0] = -0.5
v[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 0] = \
  -0.5 + \
  0.5 * (pos[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 1] + d - 0.25) / d
v[(pos[:, 1] <= 0.75 - d) & (pos[:, 1] >= 0.25 + d), 0] = 0.5
v[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 0] = \
  0.5 - \
  0.5 * (pos[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 1] + d - 0.75) / d
v[pos[:, 1] >= 0.75 + d, 0] = -0.5

v[:, 1] = B * np.sin(4. * np.pi * pos[:, 0]) * \
          (np.exp(-(pos[:, 1] - 0.25)**2 / 32. / d**2) + \
           np.exp(-(pos[:, 1] - 0.75)**2 / 32. / d**2))
            
#File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 1.]
grp.attrs["NumPart_Total"] =  [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

#Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data = pos, dtype = 'd')
grp.create_dataset("Velocities", data = v, dtype = 'f')
grp.create_dataset("Masses", data = m, dtype = 'f')
grp.create_dataset("SmoothingLength", data = h, dtype = 'f')
grp.create_dataset("InternalEnergy", data = u, dtype = 'f')
grp.create_dataset("ParticleIDs", data = ids, dtype = 'L')

fileOutput.close()
