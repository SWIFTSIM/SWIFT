###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Generates a swift IC file for the isentropic (Yee et al. 2000) vortex in a
# periodic box, in the form used by Pakmor et al. 2016 (MNRAS 455, 1134,
# appendix A): a stationary vortex with
#   T(r)   = T_inf - (gamma - 1) beta^2 / (8 gamma pi^2) exp(1 - r^2)
#   rho(r) = T^(1 / (gamma - 1)),  u = T / (gamma - 1),
#   v      = beta / (2 pi) exp((1 - r^2) / 2) (-y, x)
# centred on the box, with T_inf = 1, gamma = 7/5 and beta = 5. The paper's
# domain is [-5, 5]^2; ours is [0, 10]^2.
#
# Particles come from a 2D glass (glassPlane_32/64/128.hdf5, tiled for larger
# N) scaled to the box, with particle masses rho(x) L^2 / N^2 so that the
# particle spacing is uniform and the density is carried by the masses.
#
# Usage: python3 makeIC.py [N]   (N = 32, 64, 128 (default), 256, 512, ...)

import sys
import numpy as np
import h5py

# Parameters
gamma = 7.0 / 5.0  # Gas adiabatic index
beta = 5.0  # Vortex strength
T_inf = 1.0  # Temperature far from the vortex
L = 10.0  # Box size
N = int(sys.argv[1]) if len(sys.argv) > 1 else 128  # Particles per dimension
fileOutputName = "yeeVortex.hdf5"
# ---------------------------------------------------

base = min(N, 128)
tiles = N // base
if tiles * base != N:
    raise ValueError("N must be 32, 64 or a multiple of 128")
fileGlass = "glassPlane_%d.hdf5" % base

with h5py.File(fileGlass, "r") as f:
    glass = f["/PartType0/Coordinates"][:, :2] / f["/Header"].attrs["BoxSize"][0]
    h_glass = f["/PartType0/SmoothingLength"][:]

pos = np.concatenate(
    [(glass + [i, j]) / tiles for i in range(tiles) for j in range(tiles)]
)
h = np.tile(h_glass / tiles, tiles * tiles)
pos *= L
h *= L
numPart = len(pos)
if numPart != N * N:
    raise RuntimeError("Unexpected particle number")

x = pos[:, 0] - L / 2
y = pos[:, 1] - L / 2
r2 = x * x + y * y
T = T_inf - (gamma - 1.0) * beta**2 / (8.0 * gamma * np.pi**2) * np.exp(1.0 - r2)
rho = T ** (1.0 / (gamma - 1.0))
u = T / (gamma - 1.0)
vphi_over_r = beta / (2.0 * np.pi) * np.exp((1.0 - r2) / 2.0)
v = np.zeros((numPart, 3))
v[:, 0] = -y * vphi_over_r
v[:, 1] = x * vphi_over_r
m = rho * L * L / numPart

coords = np.zeros((numPart, 3))
coords[:, :2] = pos

# File
fileOutput = h5py.File(fileOutputName, "w")

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, 1.0]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
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
grp.create_dataset("Coordinates", data=coords, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=np.arange(1, numPart + 1), dtype="L")

fileOutput.close()

print(
    "%s: %d^2 particles, rho %.3f..%.3f, max |v| %.3f, c_s(inf) %.3f"
    % (
        fileOutputName,
        N,
        rho.min(),
        rho.max(),
        np.sqrt((v * v).sum(1)).max(),
        np.sqrt(gamma * T_inf),
    )
)
