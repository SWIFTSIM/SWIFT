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
"""u and P profiles of every snapshot, width of the interfaces and the effective conductivity."""
import glob

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402
from scipy.special import erf  # noqa: E402

gamma = 5.0 / 3.0
with h5py.File("thermalDiffusion.hdf5", "r") as f:
    C = float(f["/Header"].attrs["ThermalDiffusionContrast"])
    P0 = float(f["/Header"].attrs["ThermalDiffusionPressure"])
    x1, x2 = f["/Header"].attrs["ThermalDiffusionInterfaces"]
    Lx = float(np.atleast_1d(f["/Header"].attrs["BoxSize"])[0])
u_hot = P0 / (gamma - 1.0)
u_cold = u_hot / C


def width(x, u, xi, sign):
    """Width sigma of the erf step of u at the interface xi (sign: +1 if u rises with x)."""
    sel = np.abs(x - xi) < 0.3
    theta = (u[sel] - u_cold) / (u_hot - u_cold)
    model = lambda xx, s, x0: 0.5 * (1.0 + sign * erf((xx - x0) / (np.sqrt(2.0) * s)))
    p, _ = curve_fit(model, x[sel], theta, p0=[0.02, xi], bounds=([1e-4, xi - 0.1], [0.3, xi + 0.1]))
    return p[0]


files = sorted(glob.glob("thermalDiffusion_[0-9][0-9][0-9][0-9].hdf5"))
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
t, sig, blip, vmax = [], [], [], []
for i, fn in enumerate(files):
    with h5py.File(fn, "r") as f:
        time = float(f["/Header"].attrs["Time"])
        x = f["/PartType0/Coordinates"][:, 0]
        u = f["/PartType0/InternalEnergies"][:]
        P = f["/PartType0/Pressures"][:]
        v = f["/PartType0/Velocities"][:]
    cs = np.sqrt(gamma * P0 / C)
    t.append(time)
    sig.append(0.5 * (width(x, u, x1, -1) + width(x, u, x2, +1)))
    blip.append(np.abs(P / P0 - 1.0).max())
    vmax.append(np.sqrt((v**2).sum(1)).max() / cs)
    c = plt.cm.viridis(i / max(len(files) - 1, 1))
    xb = np.linspace(0, Lx, 201)
    ub = np.array([np.mean(u[(x >= a) & (x < b)]) for a, b in zip(xb[:-1], xb[1:])])
    Pb = np.array([np.mean(P[(x >= a) & (x < b)]) for a, b in zip(xb[:-1], xb[1:])])
    ax[0].plot(0.5 * (xb[1:] + xb[:-1]), ub / u_hot, color=c, lw=1, label=f"t = {time:.2f}")
    ax[1].plot(0.5 * (xb[1:] + xb[:-1]), Pb / P0, color=c, lw=1)
t, sig = np.array(t), np.array(sig)
kappa = np.polyfit(t, sig**2, 1)[0] / 2.0
ax[0].set_xlabel("x")
ax[0].set_ylabel(r"$u / u_{\rm hot}$")
ax[0].legend(fontsize=6)
ax[1].set_xlabel("x")
ax[1].set_ylabel(r"$P / P_0$")
ax[2].plot(t, sig, "o-", label=rf"$\sigma$, $\kappa_{{\rm eff}}$ = {kappa:.2e}")
ax[2].plot(t, blip, "s-", label=r"max $|P/P_0 - 1|$")
ax[2].plot(t, vmax, "^-", label=r"max $|v| / c_{s,\rm cold}$")
ax[2].set_xlabel("t")
ax[2].set_yscale("log")
ax[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig("thermalDiffusion.png", dpi=200)
print(f"interface width {sig[0]:.4f} -> {sig[-1]:.4f}, kappa_eff = {kappa:.3e}, "
      f"pressure blip max {max(blip):.3f}, max |v| / cs_cold {max(vmax):.3f}")
