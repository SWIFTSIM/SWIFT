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

"""Velocity power spectrum and Mach-number history of the driven-turbulence
box, compared with the curves of Bauer & Springel (2012). The spectrum of every
snapshot with t >= 10 is computed with power_spectrum.py (cached in spectra/)
and averaged, as in the paper. The reference curves in reference/ were
digitised from figures 2, 5 and 6 of the paper (bs12_*.txt; Arepo = moving mesh,
SPH = alpha = 1 artificial viscosity, SPH-tav = time-variable viscosity).

Usage: python3 plotSolution.py [snapshot basename] [t_start]
(defaults: drivenTurbulence, 10). Writes drivenTurbulence.png.
"""

import sys
import os
import glob
import subprocess
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

basename = sys.argv[1] if len(sys.argv) > 1 else "drivenTurbulence"
here = os.path.dirname(os.path.abspath(__file__))
t_start = (
    float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
)  # start of the averaging window (quasi-stationary state)

# Spectra of the snapshots in the averaging window
os.makedirs("spectra", exist_ok=True)
rows, N = [], None
for fn in sorted(glob.glob(f"{basename}_*.hdf5")):
    with h5py.File(fn, "r") as f:
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
        N = int(round(f["PartType0/Coordinates"].shape[0] ** (1.0 / 3.0)))
    if t < t_start - 1e-6:
        continue
    out = f"spectra/spectrum_{os.path.basename(fn)[len(basename) + 1:-5]}.txt"
    if not os.path.exists(out):
        subprocess.run(
            [sys.executable, f"{here}/power_spectrum.py", fn, out], check=True
        )
    rows.append(np.loadtxt(out))
if not rows:
    sys.exit(f"No snapshot with t >= {t_start} found")
k = rows[0][:, 0]
kE_all = np.array([r[:, 2] for r in rows])
kE = kE_all.mean(axis=0)

# Mach-number history from the statistics file (c_s = 1): M = sqrt(2 E_kin / M_gas)
stats = np.loadtxt("statistics.txt").T
t_stat, mach = stats[1], np.sqrt(2.0 * stats[13] / stats[4])

# Reference curves: same resolution if the paper has it (64^3: figure 6, 256^3: figure 5), else 256^3
ref = here
if N == 64:
    refs = [
        ("Arepo $64^3$ (A1)", "bs12_fig6top_A1.txt"),
        ("SPH, var. visc. $64^3$ (S1-tav)", "bs12_fig6bot_S1-tav.txt"),
        ("SPH, $\\alpha=1$ $64^3$ (S1)", "bs12_fig6bot_S1.txt"),
    ]
else:
    refs = [
        ("Arepo $256^3$ (A3)", "bs12_fig5_A3.txt"),
        ("SPH, var. visc. $256^3$ (S3-tav)", "bs12_fig5_S3-tav.txt"),
        ("SPH, $\\alpha=1$ $256^3$ (S3)", "bs12_fig5_S3.txt"),
    ]
colors = ["#1baf7a", "#6250d6", "#e87ba4"]

fig, (ax, axm) = plt.subplots(
    1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [1.4, 1]}
)
ax.axvspan(6.27, 12.57, color="#f0efec", zorder=0)
ax.axvline(np.pi * N, color="#52514e", lw=0.8, ls=":")
for (label, fn), c in zip(refs, colors):
    d = np.loadtxt(f"{ref}/{fn}")
    d = d[np.argsort(d[:, 0])]
    ax.plot(d[:, 0], d[:, 1], ls="--", lw=1.5, color=c, label=label)
kk = np.geomspace(7, 300, 10)
ax.plot(
    kk, 0.05 * (kk / 7.0) ** (-2.0 / 3.0), ls=":", lw=1, color="k", label="$k^{-5/3}$"
)
ax.fill_between(
    k, kE_all.min(axis=0), kE_all.max(axis=0), color="#2a78d6", alpha=0.2, lw=0
)
ax.plot(k, kE, lw=2.2, color="#2a78d6", label=f"SWIFT ${N}^3$ ({len(rows)} snapshots)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(5, 1.6 * np.pi * N)
ax.set_ylim(1e-6, 2e-1)
ax.set_xlabel("$k$")
ax.set_ylabel("$k\\,E(k)$")
ax.set_title(
    f"Velocity power spectrum, average over $t \\geq {t_start:g}$", fontsize=10
)
ax.legend(fontsize=8, frameon=False, loc="lower left")

axm.axvspan(t_start, t_stat.max(), color="#f0efec", zorder=0)
for (label, fn), c in zip(
    [
        ("Arepo (A3)", "bs12_fig2_A3.txt"),
        ("SPH, var. visc. (S3-tav)", "bs12_fig2_S3-tav.txt"),
        ("SPH, $\\alpha=1$ (S3)", "bs12_fig2_S3.txt"),
    ],
    colors,
):
    d = np.loadtxt(f"{ref}/{fn}")
    axm.plot(d[:, 0], d[:, 1], ls="--", lw=1.3, color=c, label=label)
axm.plot(t_stat, mach, lw=1.8, color="#2a78d6", label="SWIFT")
axm.set_xlim(0, t_stat.max())
axm.set_ylim(0, 0.5)
axm.set_xlabel("$t$")
axm.set_ylabel("rms Mach number")
axm.set_title("Mach number (dashed: paper, $256^3$)", fontsize=10)
axm.legend(fontsize=8, frameon=False, loc="lower right")

fig.tight_layout()
fig.savefig("drivenTurbulence.png", dpi=150)
sel = t_stat >= t_start
print(
    f"Mean rms Mach number for t >= {t_start:g}: {mach[sel].mean():.3f} (paper: Arepo 0.33, SPH-tav 0.30, SPH 0.25)"
)
