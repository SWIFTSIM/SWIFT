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

# Convergence of the Yee vortex: L1 errors at the final time against the
# resolution, compared with the four Arepo variants of Pakmor et al. 2016
# (MNRAS 455, 1134, figure 4; curves digitised from the paper), plus the
# kinetic-energy history of every run.
#
# Usage: python3 plotConvergence.py [--label NAME] run_N32/yeeVortex_0004.hdf5 run_N64/yeeVortex_0004.hdf5 ...
#
# Each argument is the final snapshot of one run; the t = 0 snapshot
# (<basename>_0000.hdf5) and statistics.txt are read from the same directory.
# The L1 norm is (1/V) sum_i |f_i - f(x_i)| m_i / rho_i (their equation 13).
# The "profile" error is the same norm applied to the azimuthally averaged
# density (100 radial bins to r = 5): it removes the particle noise of the
# SPH density estimate, which on a glass is a floor of about 2e-3 in L1(rho)
# independent of N, and isolates the systematic error.

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
import os
import sys

args = sys.argv[1:]
label = "SWIFT"
if "--label" in args:
    i = args.index("--label")
    label = args[i + 1]
    del args[i : i + 2]
here = os.path.dirname(os.path.abspath(__file__))


def analytic(pos, L):
    x = pos[:, 0] - L / 2
    y = pos[:, 1] - L / 2
    r2 = x * x + y * y
    T = T_inf - (gamma - 1.0) * beta**2 / (8.0 * gamma * np.pi**2) * np.exp(1.0 - r2)
    vf = beta / (2.0 * np.pi) * np.exp((1.0 - r2) / 2.0)
    return np.sqrt(r2), T ** (1.0 / (gamma - 1.0)), np.c_[-y * vf, x * vf]


def errors(filename):
    with h5py.File(filename, "r") as f:
        L = f["/Header"].attrs["BoxSize"][0]
        t = float(np.ravel(f["/Header"].attrs["Time"])[0])
        pos = f["/PartType0/Coordinates"][:, :2]
        v = f["/PartType0/Velocities"][:, :2]
        m = f["/PartType0/Masses"][:]
        rho = f["/PartType0/Densities"][:]
    r, rho_a, v_a = analytic(pos, L)
    V = m / rho
    Vtot = L * L
    rb = np.linspace(0, 5, 101)
    ib = np.clip(np.digitize(r, rb) - 1, 0, 99)
    cnt = np.bincount(ib, minlength=100)
    prof = lambda q: np.bincount(ib, weights=q, minlength=100) / np.maximum(cnt, 1)
    area = np.pi * (rb[1:] ** 2 - rb[:-1] ** 2)
    ok = cnt > 0
    return dict(
        t=t,
        N=int(round(np.sqrt(len(m)))),
        L1_rho=np.sum(np.abs(rho - rho_a) * V) / Vtot,
        L1_v=np.sum(np.sqrt(((v - v_a) ** 2).sum(1)) * V) / Vtot,
        L1_prof=np.sum((np.abs(prof(rho) - prof(rho_a)) * area)[ok]) / Vtot,
        Ekin=0.5 * np.sum(m * (v * v).sum(1)),
    )


runs = []
for fn in args:
    d = os.path.dirname(fn) or "."
    base = os.path.basename(fn).rsplit("_", 1)[0]
    e1 = errors(fn)
    e0 = errors(os.path.join(d, base + "_0000.hdf5"))
    stats = os.path.join(d, "statistics.txt")
    ek = np.loadtxt(stats)[:, [1, 13]] if os.path.exists(stats) else None
    runs.append((e0, e1, ek))
    print(
        "N = %4d, t = %.1f: L1(rho) %.3e (t = 0: %.3e), profile L1(rho) %.3e, L1(v) %.3e, Ekin/Ekin(0) %.4f"
        % (
            e1["N"],
            e1["t"],
            e1["L1_rho"],
            e0["L1_rho"],
            e1["L1_prof"],
            e1["L1_v"],
            e1["Ekin"] / e0["Ekin"],
        )
    )
runs.sort(key=lambda r: r[1]["N"])
N = np.array([e1["N"] for e0, e1, ek in runs])
L1 = np.array([e1["L1_rho"] for e0, e1, ek in runs])
L10 = np.array([e0["L1_rho"] for e0, e1, ek in runs])
L1p = np.array([e1["L1_prof"] for e0, e1, ek in runs])
L1v = np.array([e1["L1_v"] for e0, e1, ek in runs])
slope = lambda y: np.polyfit(np.log(N), np.log(y), 1)[0] if len(N) > 1 else np.nan
print(
    "slopes: L1(rho) %.2f, profile L1(rho) %.2f, L1(v) %.2f"
    % (slope(L1), slope(L1p), slope(L1v))
)

plt.rcParams.update({"font.size": 9, "figure.figsize": (11, 4.5)})
fig, (ax, bx) = plt.subplots(1, 2)
ax.set_xscale("log")
ax.set_yscale("log")
for name, c in (("MH-GG", "0.5"), ("RK-GG", "C0"), ("MH-LSF", "C1"), ("RK-LSF", "C2")):
    fn = os.path.join(here, "pakmor2016_fig4_%s.txt" % name)
    if os.path.exists(fn):
        d = np.loadtxt(fn)
        ax.plot(
            d[:, 0],
            d[:, 1],
            "o--",
            ms=3,
            lw=1,
            color=c,
            label="Arepo %s (Pakmor+16)" % name.replace("-", " "),
        )
ax.plot(
    N,
    L1,
    "s-",
    color="k",
    lw=2,
    ms=6,
    label="%s $L_1(\\rho)$ (slope %.2f)" % (label, slope(L1)),
)
ax.plot(
    N,
    L10,
    "s:",
    color="k",
    lw=1,
    ms=4,
    mfc="none",
    label="$L_1(\\rho)$ at $t=0$ (SPH estimate of the IC)",
)
ax.plot(
    N,
    L1p,
    "D-",
    color="C3",
    lw=1.5,
    ms=5,
    label="%s azimuthal-profile $L_1(\\rho)$ (slope %.2f)" % (label, slope(L1p)),
)
ax.plot(
    N,
    L1v,
    "^-",
    color="C4",
    lw=1.5,
    ms=5,
    label="%s $L_1(v)$ (slope %.2f)" % (label, slope(L1v)),
)
for p, lab in ((1, "$N^{-1}$"), (2, "$N^{-2}$")):
    ax.plot(
        [30, 3000],
        [L1[0] * (30.0 / N[0]) ** -p, L1[0] * (3000.0 / N[0]) ** -p],
        "k:",
        lw=0.7,
    )
    ax.text(
        2500, L1[0] * (2500.0 / N[0]) ** -p, lab, fontsize=8, ha="right", va="bottom"
    )
ax.set_xlabel("$N$ (particles or cells per dimension)")
ax.set_ylabel("$L_1$ at $t = %.0f$" % runs[0][1]["t"])
ax.set_xlim(25, 3000)
ax.set_ylim(1e-7, 5e-2)
ax.legend(fontsize=7, frameon=False, loc="lower left")
ax.set_title("Yee vortex: convergence (Pakmor et al. 2016, fig. 4)")
for e0, e1, ek in runs:
    if ek is not None:
        bx.plot(ek[:, 0], ek[:, 1] / ek[0, 1], lw=1.5, label="$N = %d$" % e1["N"])
bx.set_xlabel("$t$")
bx.set_ylabel("$E_{\\rm kin}(t) / E_{\\rm kin}(0)$")
bx.set_title("kinetic energy of the vortex (the loss is thermalised)")
bx.legend(fontsize=8, frameon=False)
fig.tight_layout()
fig.savefig("YeeVortex_convergence.png", dpi=200)
