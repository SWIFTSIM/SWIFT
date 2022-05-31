################################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
gamma = 5.0 / 3.0  # Gas adiabatic index
P0 = 2.5  # Pressure
rho0 = 1.0  # Density
d = 0.0317  # Thickness of the transition layer
B = 0.0005  # Amplitude of the seed velocity
N1D = 128  # Number of particles in one dimension

fileOutputName = "kelvinHelmholtzGrowthRate.hdf5"

# ---------------------------------------------------

N = N1D ** 2
x = np.linspace(0.0, 1.0, N1D + 1)
x = 0.5 * (x[1:] + x[:-1])
y = x
xx, yy = np.meshgrid(x, y)
pos = np.zeros((N, 3))
pos[:, 0] = xx.reshape((N))
pos[:, 1] = yy.reshape((N))
h = np.ones(N) * 2.0 / N1D

vol = 1.0

# Generate extra arrays
v = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))

v[pos[:, 1] <= 0.25 - d, 0] = -0.5
v[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 0] = (
    -0.5
    + 0.5 * (pos[(pos[:, 1] < 0.25 + d) & (pos[:, 1] > 0.25 - d), 1] + d - 0.25) / d
)
v[(pos[:, 1] <= 0.75 - d) & (pos[:, 1] >= 0.25 + d), 0] = 0.5
v[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 0] = (
    0.5 - 0.5 * (pos[(pos[:, 1] < 0.75 + d) & (pos[:, 1] > 0.75 - d), 1] + d - 0.75) / d
)
v[pos[:, 1] >= 0.75 + d, 0] = -0.5

v[:, 1] = (
    B
    * np.sin(4.0 * np.pi * pos[:, 0])
    * (
        np.exp(-(pos[:, 1] - 0.25) ** 2 / 32.0 / d ** 2)
        + np.exp(-(pos[:, 1] - 0.75) ** 2 / 32.0 / d ** 2)
    )
)

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [1.0, 1.0, 1.0]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

fileOutput.close()
