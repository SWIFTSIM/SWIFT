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
"""Thermal (entropy) discontinuity in pressure equilibrium, 2D glass.

A periodic box of n copies of a 2D glass along x (box n x 1) with a dense,
cold slab in its middle half: density contrast C across the two interfaces,
the same pressure P0 = 1 everywhere and no velocity. The particles keep the
spacing of the glass and carry the density contrast in their masses
(m_cold = C m_hot), so that the interfaces are not resolution jumps; the
internal energies are set from the SPH density of the particles (Wendland C4
kernel, --eta) so that the pressure is uniform at the start to the accuracy of
the kernel-sum density.

The exact solution is static. What the run measures is the artificial thermal
conduction of the scheme: the growth of the width of the two u profiles
(an erf of width sigma grows as sigma^2 = sigma0^2 + 2 kappa_eff t), the
pressure blip that develops at the interfaces and the spurious velocities it
drives. A scheme with no conduction keeps the profiles and the blip; too much
conduction smears the contact.
"""
import argparse

import h5py
import numpy as np
from scipy.spatial import cKDTree

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--glass", default="glassPlane_128.hdf5", help="2D glass file")
parser.add_argument("--ncopies", type=int, default=2, help="copies of the glass along x")
parser.add_argument("--contrast", type=float, default=4.0, help="density contrast of the slab")
parser.add_argument("--eta", type=float, default=1.5, help="resolution eta for the SPH density of the ICs")
parser.add_argument("--nominal", action="store_true", help="internal energies from the nominal (step) density instead of the SPH density")
parser.add_argument("-o", "--output", default="thermalDiffusion.hdf5")
args = parser.parse_args()

gamma = 5.0 / 3.0
P0 = 1.0
rho_hot = 1.0
rho_cold = args.contrast * rho_hot

with h5py.File(args.glass, "r") as glass:
    pos1 = glass["/PartType0/Coordinates"][:, :]
    h1 = glass["/PartType0/SmoothingLength"][:]
    L1 = float(np.atleast_1d(glass["/Header"].attrs["BoxSize"])[0])
n1 = len(h1)
Lx = args.ncopies * L1
pos = np.concatenate([pos1 + np.array([i * L1, 0.0, 0.0]) for i in range(args.ncopies)])
pos[:, 2] = 0.0
# Keep the positions inside [0, L) for the periodic tree
pos[:, 0] = np.mod(pos[:, 0], Lx)
pos[:, 1] = np.mod(pos[:, 1], L1)
pos[:, 0][pos[:, 0] >= Lx] = 0.0
pos[:, 1][pos[:, 1] >= L1] = 0.0
h = np.tile(h1, args.ncopies)
N = len(h)
x = pos[:, 0]

# Dense, cold slab in the middle half of the box
inside = (x > 0.25 * Lx) & (x < 0.75 * Lx)
m_hot = rho_hot * L1 * L1 / n1
m = np.where(inside, args.contrast * m_hot, m_hot)
rho_nominal = np.where(inside, rho_cold, rho_hot)

# SPH density of this particle set with the Wendland C4 kernel in 2D
# (SWIFT convention: support H = gamma_k h, h = eta (m / rho)^(1/2))
kernel_gamma = 2.1712405933672376  # sqrt(33/7)
kernel_norm = 9.0 / np.pi


def w_c4(q):
    q = np.minimum(q, 1.0)
    return (1.0 - q) ** 6 * (1.0 + 6.0 * q + 35.0 / 3.0 * q * q)


def sph_density(h_ngb):
    H = kernel_gamma * h_ngb
    tree = cKDTree(pos[:, :2], boxsize=[Lx, L1])
    rho = np.zeros(N)
    pairs = tree.query_ball_point(pos[:, :2], r=H.max())
    for i, nb in enumerate(pairs):
        nb = np.asarray(nb)
        d = pos[nb, :2] - pos[i, :2]
        d -= np.array([Lx, L1]) * np.round(d / np.array([Lx, L1]))
        r = np.hypot(d[:, 0], d[:, 1])
        rho[i] = np.sum(m[nb] * w_c4(r / H[i])) * kernel_norm / H[i] ** 2
    return rho


if args.nominal:
    rho = rho_nominal
    h_ic = args.eta * np.sqrt(m / rho)
else:
    rho = rho_nominal.copy()
    for _ in range(3):
        h_ic = args.eta * np.sqrt(m / rho)
        rho = sph_density(h_ic)
u = P0 / ((gamma - 1.0) * rho)
v = np.zeros((N, 3))
ids = np.arange(1, N + 1)

with h5py.File(args.output, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [Lx, L1, 1.0]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2
    grp.attrs["ThermalDiffusionContrast"] = args.contrast
    grp.attrs["ThermalDiffusionPressure"] = P0
    grp.attrs["ThermalDiffusionInterfaces"] = [0.25 * Lx, 0.75 * Lx]
    grp = f.create_group("/Units")
    for key in ["Unit length in cgs (U_L)", "Unit mass in cgs (U_M)", "Unit time in cgs (U_t)",
                "Unit current in cgs (U_I)", "Unit temperature in cgs (U_T)"]:
        grp.attrs[key] = 1.0
    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=v, dtype="f")
    grp.create_dataset("Masses", data=m, dtype="f")
    grp.create_dataset("SmoothingLength", data=h_ic, dtype="f")
    grp.create_dataset("InternalEnergy", data=u, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
blip = np.abs((gamma - 1.0) * rho * u / P0 - 1.0).max()
print(f"thermal discontinuity: {N} particles in a {Lx} x {L1} box, contrast {args.contrast}, "
      f"u from the {'nominal' if args.nominal else 'SPH'} density (max |P / P0 - 1| = {blip:.1e})")
