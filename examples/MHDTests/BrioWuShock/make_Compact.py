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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
gamma = 5.0 / 3.0  # Gas adiabatic index
P0 = 2.5  # Pressure
rho0 = 1.0  # Density
N1D = 20  # Number of particles in one dimension

N1D = 34
fileOutputName = "HCP_low.hdf5"
# N1D = 50
# fileOutputName = "HCP__hi.hdf5"


# ---------------------------------------------------

xx, yy, zz = 0, 0, 0
i, j, k = 0, 0, 0
ax = [0.0]
ay = [0.0]
az = [0.0]
for i in range(1, 60):
    for j in range(0, 60):
        for k in range(0, 60):
            #            idx = k*N1D*N1D+j*N1D+i
            xx = (2.0 * i + (j + k) % 2) / N1D
            yy = (np.sqrt(3) * (j + 1.0 / 3.0 * (k % 2))) / N1D
            zz = (2.0 * np.sqrt(6.0) / 3.0 * k) / N1D
            if xx <= 1.0 and yy <= 1.0 and zz <= 1.0:
                ax.append(xx)
                ay.append(yy)
                az.append(zz)

ax.pop(0)
ay.pop(0)
az.pop(0)
N = len(ax[1:])
pos = np.zeros((N, 3))
pos[:, 0] = ax[1:]  # .reshape((N))
pos[:, 1] = ay[1:]  # .reshape((N))
pos[:, 2] = az[1:]  # .reshape((N))
fig, axx = pl.subplots(3, 1, sharex=True, figsize=(3, 9))
axx[0].plot(pos[:, 0], pos[:, 1], ".")
axx[0].set_ylabel("y")
axx[0].set_xlabel("x")
axx[1].plot(pos[:, 0], pos[:, 2], ".")
axx[1].set_ylabel("z")
axx[1].set_xlabel("x")
axx[2].plot(pos[:, 1], pos[:, 2], ".")
axx[2].set_ylabel("z")
axx[2].set_xlabel("y")
pl.savefig("a.png", dpi=300)
NN = N1D ** 3
print(N, NN)
print(min(pos[:, 0]), max(pos[:, 0]))
print(min(pos[:, 1]), max(pos[:, 1]))
print(min(pos[:, 2]), max(pos[:, 2]))

h = np.ones(N) * 2.0 / (N ** (1 / 3))

vol = 1.0

# Generate extra arrays
v = np.zeros((N, 3))
ids = np.linspace(1, N, N)
m = np.ones(N) * rho0 * vol / N
u = np.ones(N) * P0 / (rho0 * (gamma - 1.0))

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
grp.attrs["Dimension"] = 3

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
