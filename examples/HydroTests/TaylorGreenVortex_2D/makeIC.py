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
"""Taylor-Green vortex (Taylor & Green 1937) in 2D on a glass.

v_x = v0 sin(kx) cos(ky), v_y = -v0 cos(kx) sin(ky), k = 2 pi, uniform density
and P = P0 + rho v0^2 / 4 (cos 2kx + cos 2ky), which makes the flow an exact
stationary solution of the Euler equations. With a kinematic viscosity nu the
kinetic energy decays as exp(-4 nu k^2 t) = exp(-16 pi^2 nu t) (Phantom paper,
section 5.4.1), so the decay measured without physical viscosity gives the
effective numerical viscosity of the scheme in a smooth subsonic flow.
"""
import argparse

import h5py
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--glass", default="glassPlane_128.hdf5", help="2D glass file")
parser.add_argument("--v0", type=float, default=0.1, help="velocity amplitude (Mach v0 / cs)")
parser.add_argument("--cs", type=float, default=1.0, help="sound speed of the unperturbed gas")
parser.add_argument("-o", "--output", default="taylorGreen.hdf5")
args = parser.parse_args()

gamma = 5.0 / 3.0
rho0 = 1.0
k = 2.0 * np.pi
P0 = rho0 * args.cs**2 / gamma

with h5py.File(args.glass, "r") as glass:
    pos = glass["/PartType0/Coordinates"][:, :]
    h = glass["/PartType0/SmoothingLength"][:]
N = len(h)
x = pos[:, 0]
y = pos[:, 1]

v = np.zeros((N, 3))
v[:, 0] = args.v0 * np.sin(k * x) * np.cos(k * y)
v[:, 1] = -args.v0 * np.cos(k * x) * np.sin(k * y)
P = P0 + 0.25 * rho0 * args.v0**2 * (np.cos(2.0 * k * x) + np.cos(2.0 * k * y))
m = np.ones(N) * rho0 / N
u = P / (rho0 * (gamma - 1.0))
ids = np.arange(1, N + 1)

with h5py.File(args.output, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [1.0, 1.0, 1.0]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2
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
print(f"Taylor-Green vortex: {N} particles, v0 = {args.v0}, cs = {args.cs} (Mach {args.v0 / args.cs:.2f})")
