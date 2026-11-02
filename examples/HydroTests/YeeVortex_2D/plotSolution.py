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

# Plots the SPH answer for the Yee vortex against the analytic (stationary)
# solution and prints the L1 errors. Usage: python3 plotSolution.py <snapshot number>

# Parameters (must match makeIC.py)
gamma = 7.0 / 5.0
beta = 5.0
T_inf = 1.0

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import h5py
import sys

snap = int(sys.argv[1])
basename = sys.argv[2] if len(sys.argv) > 2 else "yeeVortex"
filename = "%s_%04d.hdf5" % (basename, snap)


def analytic(r2):
    """rho, u, P and v_phi of the vortex at r^2."""
    T = T_inf - (gamma - 1.0) * beta**2 / (8.0 * gamma * np.pi**2) * np.exp(1.0 - r2)
    rho = T ** (1.0 / (gamma - 1.0))
    vphi = beta / (2.0 * np.pi) * np.exp((1.0 - r2) / 2.0) * np.sqrt(r2)
    return rho, T / (gamma - 1.0), rho * T, vphi


with h5py.File(filename, "r") as f:
    L = f["/Header"].attrs["BoxSize"][0]
    time = float(np.ravel(f["/Header"].attrs["Time"])[0])
    scheme = f["/HydroScheme"].attrs["Scheme"].decode("utf-8")
    kernel = f["/HydroScheme"].attrs["Kernel function"].decode("utf-8")
    pos = f["/PartType0/Coordinates"][:, :2]
    vel = f["/PartType0/Velocities"][:, :2]
    m = f["/PartType0/Masses"][:]
    rho = f["/PartType0/Densities"][:]
    u = f["/PartType0/InternalEnergies"][:]
    P = f["/PartType0/Pressures"][:]

x = pos[:, 0] - L / 2
y = pos[:, 1] - L / 2
r = np.sqrt(x * x + y * y)
vphi = (-y * vel[:, 0] + x * vel[:, 1]) / np.maximum(r, 1e-12)
rho_a, u_a, P_a, vphi_a = analytic(r * r)
vx_a = -y * vphi_a / np.maximum(r, 1e-12)
vy_a = x * vphi_a / np.maximum(r, 1e-12)

# L1 errors, volume weighted with V_i = m_i / rho_i (Pakmor et al. 2016, eq. 13)
V = m / rho
Vtot = L * L
L1 = {
    "rho": np.sum(np.abs(rho - rho_a) * V) / Vtot,
    "v": np.sum(np.hypot(vel[:, 0] - vx_a, vel[:, 1] - vy_a) * V) / Vtot,
    "u": np.sum(np.abs(u - u_a) * V) / Vtot,
    "P": np.sum(np.abs(P - P_a) * V) / Vtot,
}
N = int(round(np.sqrt(len(m))))
print("%s (t = %.2f, %d^2 particles, %s, %s):" % (filename, time, N, scheme, kernel))
for k in ("rho", "v", "u", "P"):
    print("  L1(%s) = %.3e" % (k, L1[k]))

rr = np.linspace(0, 5, 400)
rho_r, u_r, P_r, vphi_r = analytic(rr * rr)

plt.rcParams.update({"font.size": 10, "figure.figsize": (9, 7)})
fig, axes = plt.subplots(2, 2, sharex=True)
for ax, q, qa, lab, l1 in zip(
    axes.ravel(),
    (rho, vphi, u, P),
    (rho_r, vphi_r, u_r, P_r),
    (
        "Density $\\rho$",
        "Azimuthal velocity $v_\\phi$",
        "Internal energy $u$",
        "Pressure $P$",
    ),
    (L1["rho"], L1["v"], L1["u"], L1["P"]),
):
    ax.scatter(r, q, s=0.5, c="C0", alpha=0.4, rasterized=True)
    ax.plot(rr, qa, "k-", lw=1.0, label="analytic")
    ax.set_ylabel(lab)
    ax.set_xlim(0, 5)
    ax.text(
        0.97, 0.95, "$L_1$ = %.2e" % l1, transform=ax.transAxes, ha="right", va="top"
    )
for ax in axes[1]:
    ax.set_xlabel("$r$")
axes[0, 0].legend(loc="lower right", frameon=False)
fig.suptitle(
    "Yee vortex, $t = %.1f$, $%d^2$ particles, %s (%s)" % (time, N, scheme, kernel),
    fontsize=10,
)
fig.tight_layout()
fig.savefig("YeeVortex_%04d.png" % snap, dpi=200)
