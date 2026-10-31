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

"""Initial conditions for the driven subsonic turbulence test of
Bauer & Springel (2012, MNRAS 423, 2558; arXiv:1109.4413), section 3:
a unit periodic box of gas with rho = 1 and isothermal sound speed c_s = 1,
N^3 particles on a Cartesian lattice at rest. The gas is nearly isothermal
(gamma = 1.001, see the README), so the internal energy is u = c_s^2 / (gamma - 1).

Usage: python3 makeIC.py [N] [gamma] [eta]   (defaults: 64, 1.001, 1.2348)
Writes drivenTurbulence.hdf5.
"""

import sys
import numpy as np
import h5py

N = int(sys.argv[1]) if len(sys.argv) > 1 else 64
gamma = float(sys.argv[2]) if len(sys.argv) > 2 else 1.001
eta = float(sys.argv[3]) if len(sys.argv) > 3 else 1.2348
fileName = "drivenTurbulence.hdf5"

L = 1.0  # Box size
rho = 1.0  # Density
cs = 1.0  # Isothermal sound speed
u = cs * cs / (gamma - 1.0)  # P = (gamma - 1) rho u = rho c_s^2

numPart = N**3
g = (np.arange(N) + 0.5) * L / N
x, y, z = np.meshgrid(g, g, g, indexing="ij")
pos = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)

with h5py.File(fileName, "w") as f:
    grp = f.create_group("Header")
    grp.attrs["BoxSize"] = [L, L, L]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 3

    grp = f.create_group("Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    grp = f.create_group("PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=np.zeros((numPart, 3)), dtype="f")
    grp.create_dataset("Masses", data=np.full(numPart, rho * L**3 / numPart), dtype="f")
    grp.create_dataset("SmoothingLength", data=np.full(numPart, eta * L / N), dtype="f")
    grp.create_dataset("InternalEnergy", data=np.full(numPart, u), dtype="f")
    grp.create_dataset("ParticleIDs", data=np.arange(1, numPart + 1), dtype="L")

print(f"Wrote {fileName}: {N}^3 particles, u = {u:g}, h = {eta * L / N:g}")
