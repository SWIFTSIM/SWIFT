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
"""Profiles of the strong blast wave against the exact Riemann solution."""
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "../")
from riemannSolver import RiemannSolver  # noqa: E402

snap = int(sys.argv[1]) if len(sys.argv) > 1 else 2
gamma = 5.0 / 3.0

with h5py.File(f"blastWave_{snap:04d}.hdf5", "r") as f:
    t = float(np.atleast_1d(f["/Header"].attrs["Time"])[0])
    L = float(np.atleast_1d(f["/Header"].attrs["BoxSize"])[0])
    x = f["/PartType0/Coordinates"][:, 0]
    v = f["/PartType0/Velocities"][:, 0]
    rho = f["/PartType0/Densities"][:]
    P = f["/PartType0/Pressures"][:]
    u = f["/PartType0/InternalEnergies"][:]
with h5py.File("blastWave.hdf5", "r") as f:
    u0 = f["/PartType0/InternalEnergy"][:]
P_L, P_R = (gamma - 1.0) * u0.max(), (gamma - 1.0) * u0.min()

solver = RiemannSolver(gamma)
xs = np.linspace(0, L, 4000)
rho_s, v_s, P_s, _ = solver.solve(1.0, 0.0, P_L, 1.0, 0.0, P_R, (xs - 0.5 * L) / t)
u_s = P_s / (rho_s * (gamma - 1.0))

fig, axes = plt.subplots(2, 2, figsize=(9, 6))
for ax, y, ys, lab, log in zip(axes.flat, [rho, v, P, u], [rho_s, v_s, P_s, u_s],
                               [r"$\rho$", r"$v_x$", r"$P$", r"$u$"], [False, False, True, True]):
    ax.plot(x, y, ".", ms=1, alpha=0.3, color="0.4")
    ax.plot(xs, ys, "r-", lw=1)
    ax.set_xlabel("x")
    ax.set_ylabel(lab)
    if log:
        ax.set_yscale("log")
fig.suptitle(f"Strong blast wave, t = {t:.4f}, P_L / P_R = {P_L / P_R:g}")
fig.tight_layout()
fig.savefig(f"blastWave_{snap:04d}.png", dpi=150)
