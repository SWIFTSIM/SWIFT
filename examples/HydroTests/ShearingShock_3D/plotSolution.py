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
"""Profiles of the shearing shock against the exact 1D solution of the
collision (central Riemann problem) and of the receding streams at the box
edge, and the change of each particle's v_y, which should be zero."""
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

with h5py.File(f"shearingShock_{snap:04d}.hdf5", "r") as f:
    t = float(np.atleast_1d(f["/Header"].attrs["Time"])[0])
    L = float(np.atleast_1d(f["/Header"].attrs["BoxSize"])[0])
    x = f["/PartType0/Coordinates"][:, 0]
    v = f["/PartType0/Velocities"][:, :]
    rho = f["/PartType0/Densities"][:]
    P = f["/PartType0/Pressures"][:]
    ids = f["/PartType0/ParticleIDs"][:]
with h5py.File("shearingShock.hdf5", "r") as f:
    v0 = f["/PartType0/Velocities"][:, :]
    ids0 = f["/PartType0/ParticleIDs"][:]
    u0 = f["/PartType0/InternalEnergy"][0]
dv = float(np.max(v0[:, 0]))
s = float(np.max(np.abs(v0[:, 1])))
P0 = (gamma - 1.0) * u0
vy0 = v0[np.argsort(ids0)][np.argsort(np.argsort(ids)), 1]

solver = RiemannSolver(gamma)
xs = np.linspace(0, L, 2000)
xi_c = (xs - 0.5 * L) / t
xi_e = np.where(xs < 0.5 * L, xs, xs - L) / t
rho_c, v_c, P_c, _ = solver.solve(1.0, dv, P0, 1.0, -dv, P0, xi_c)
rho_e, v_e, P_e, _ = solver.solve(1.0, -dv, P0, 1.0, dv, P0, xi_e)
central = np.abs(xs - 0.5 * L) < 0.25 * L
rho_s = np.where(central, rho_c, rho_e)
v_s = np.where(central, v_c, v_e)
P_s = np.where(central, P_c, P_e)

fig, axes = plt.subplots(2, 2, figsize=(9, 6))
for ax, y, ys, lab in zip(axes.flat, [rho, v[:, 0], P, v[:, 1] - vy0],
                          [rho_s, v_s, P_s, np.zeros_like(xs)],
                          [r"$\rho$", r"$v_x$", r"$P$", r"$v_y - v_{y,0}$"]):
    ax.plot(x, y, ".", ms=1, alpha=0.3, color="0.4")
    ax.plot(xs, ys, "r-", lw=1)
    ax.set_xlabel("x")
    ax.set_ylabel(lab)
fig.suptitle(f"Shearing shock, t = {t:.3f}, s = {s:g} dv")
fig.tight_layout()
fig.savefig(f"shearingShock_{snap:04d}.png", dpi=150)
print(f"mean |v_y - v_y0| / s = {np.mean(np.abs(v[:, 1] - vy0)) / s:.3e}")
