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
"""Strong blast wave (Rosswog 2020, section 3.4.3; Phantom paper, section
5.1.2): a planar Riemann problem with uniform density and a pressure jump of
10^4 (P = 1000 for x < L/2, 0.1 beyond), run in 3D on N unit glass cubes
along x. The shock Mach number is of order 100."""
import argparse

import h5py
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("-n", "--ncubes", type=int, default=4, help="number of glass cubes along x")
parser.add_argument("--glass", default="glassCube_32.hdf5", help="3D glass file (unit cube)")
parser.add_argument("--P_L", type=float, default=1000.0, help="left pressure")
parser.add_argument("--P_R", type=float, default=0.1, help="right pressure")
parser.add_argument("-o", "--output", default="blastWave.hdf5")
args = parser.parse_args()

gamma = 5.0 / 3.0
rho0 = 1.0

with h5py.File(args.glass, "r") as glass:
    pos_g = glass["/PartType0/Coordinates"][:, :]
    h_g = glass["/PartType0/SmoothingLength"][:]
n_g = len(h_g)
L = float(args.ncubes)
pos = np.concatenate([pos_g + np.array([i, 0.0, 0.0]) for i in range(args.ncubes)])
h = np.tile(h_g, args.ncubes)
N = len(h)
x = pos[:, 0]

P = np.where(x < 0.5 * L, args.P_L, args.P_R)
v = np.zeros((N, 3))
m = np.ones(N) * rho0 / n_g
u = P / (rho0 * (gamma - 1.0))
ids = np.arange(1, N + 1)

with h5py.File(args.output, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [L, 1.0, 1.0]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 3
    grp = f.create_group("/Units")
    for key in ["Unit length in cgs (U_L)", "Unit mass in cgs (U_M)", "Unit time in cgs (U_t)",
                "Unit current in cgs (U_I)", "Unit temperature in cgs (U_T)"]:
        grp.attrs[key] = 1.0
    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=v, dtype="f")
    grp.create_dataset("Masses", data=m, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=u, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
print(f"Strong blast wave: {N} particles in a {L:g} x 1 x 1 box, P_L / P_R = {args.P_L / args.P_R:g}")
