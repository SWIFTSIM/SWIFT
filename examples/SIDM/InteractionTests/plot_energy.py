###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
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
#################################################################################

# Energy conservation of the isolated SIDM halo (run.sh), compared with the
# same ICs run as CDM (run_cdm.sh).
#
# Scatterings are elastic, so they redistribute kinetic energy between the
# particles but must not change the total. The total energy
#
#   E = sum_i m_i |v_i|^2 / 2  +  sum_i m_i phi_i / 2
#
# is computed from each snapshot (the global statistics.txt does not include
# SIDM particles). Gravity alone does not conserve E exactly, so the drift of
# the SIDM run is compared with that of the CDM run: scattering should not add
# any drift on top of it while the core forms. The enclosed mass M(<R_CORE) is
# shown to confirm that the scatterings do change the halo over this time.
#
# Usage: python3 plot_energy.py [--sidm snap] [--cdm snap_cdm]

import os
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("--sidm", default="snap", help="Directory of the SIDM snapshots")
parser.add_argument("--cdm", default="snap_cdm", help="Directory of the CDM snapshots")
parser.add_argument("-o", "--output", default="energy_conservation.png")

R_CORE = 0.5  # radius of the enclosed mass [kpc]
TOL_DRIFT = 1e-3  # maximum allowed |E(t) - E(0)| / |E(0)| of the SIDM run

S_IN_GYR = 3.15576e16
KPC_IN_CM = 3.085677e21

C_SIDM, C_CDM = "#2a78d6", "#eb6834"  # categorical slots 1-2
INK, MUTED = "#3b3b38", "#8a8980"


def read_run(directory, ptype):
    rows = []
    for fname in sorted(glob.glob(os.path.join(directory, "snapshot_*.hdf5"))):
        with h5py.File(fname, "r") as f:
            if f["Cosmology"].attrs["Cosmological run"][0]:
                raise ValueError("This script assumes a non-cosmological run")
            u = f["Units"].attrs
            g = f[ptype]
            m = g["Masses"][:].astype(np.float64)
            v = g["Velocities"][:].astype(np.float64)
            phi = g["Potentials"][:].astype(np.float64)
            pos = g["Coordinates"][:].astype(np.float64)
            t = f["Header"].attrs["Time"][0] * u["Unit time in cgs (U_t)"][0] / S_IN_GYR
            r_core = R_CORE * KPC_IN_CM / u["Unit length in cgs (U_L)"][0]

        v_com = (m[:, None] * v).sum(axis=0) / m.sum()
        r = np.linalg.norm(pos - np.median(pos, axis=0), axis=1)
        rows.append(dict(
            t=t,
            K=0.5 * (m * ((v - v_com) ** 2).sum(axis=1)).sum(),
            W=0.5 * (m * phi).sum(),
            M_core=m[r < r_core].sum(),
        ))
    if not rows:
        return None
    run = {k: np.array([row[k] for row in rows]) for k in rows[0]}
    run["E"] = run["K"] + run["W"]
    run["dE"] = (run["E"] - run["E"][0]) / abs(run["E"][0])
    return run


def main():
    args = parser.parse_args()
    runs = {}
    for name, directory, ptype in (("SIDM", args.sidm, "SIDMParticles"),
                                   ("CDM", args.cdm, "DMParticles")):
        run = read_run(directory, ptype)
        if run is None:
            print(f"No {name} snapshots in {directory}/: skipped")
        else:
            runs[name] = run
    if "SIDM" not in runs:
        raise FileNotFoundError(f"No SIDM snapshots in {args.sidm}/")

    for name, run in runs.items():
        print(f"\n{name}")
        print(f"  {'t [Gyr]':>8} {'K':>11} {'W':>12} {'E':>12} {'dE/|E0|':>10} "
              f"{'2K/|W|':>7} {'M(<R)/M0':>9}")
        for i in range(len(run["t"])):
            print(f"  {run['t'][i]:8.4f} {run['K'][i]:11.5e} {run['W'][i]:12.5e} "
                  f"{run['E'][i]:12.5e} {run['dE'][i]:10.2e} "
                  f"{2 * run['K'][i] / abs(run['W'][i]):7.4f} "
                  f"{run['M_core'][i] / run['M_core'][0]:9.4f}")

    sidm = runs["SIDM"]
    drift = np.abs(sidm["dE"]).max()
    ok = drift < TOL_DRIFT
    print(f"\n[{'PASS' if ok else 'FAIL'}] max |dE/E0| (SIDM) = {drift:.2e} "
          f"(tolerance {TOL_DRIFT:g})")
    if "CDM" in runs:
        print(f"       max |dE/E0| (CDM)  = {np.abs(runs['CDM']['dE']).max():.2e}")
    print(f"       M(<{R_CORE:g} kpc) changed by a factor "
          f"{sidm['M_core'][-1] / sidm['M_core'][0]:.3f} (SIDM)"
          + (f", {runs['CDM']['M_core'][-1] / runs['CDM']['M_core'][0]:.3f} (CDM)"
             if "CDM" in runs else ""))

    # ------------------ Plot ------------------
    fig, (ax_e, ax_m) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    colours = {"SIDM": C_SIDM, "CDM": C_CDM}

    for name, run in runs.items():
        style = dict(color=colours[name], lw=2, marker="o", ms=8,
                     mec="white", mew=1.5, label=name)
        ax_e.plot(run["t"], run["dE"], **style)
        ax_m.plot(run["t"], run["M_core"] / run["M_core"][0], **style)

    ax_e.axhspan(-TOL_DRIFT, TOL_DRIFT, color=MUTED, alpha=0.12, lw=0,
                 label=f"$\\pm${TOL_DRIFT:g}")
    ax_e.axhline(0.0, color=INK, lw=0.8)
    ax_e.set_ylim(-2 * TOL_DRIFT, 2 * TOL_DRIFT)
    ax_e.set_ylabel(r"$(E - E_0) / |E_0|$")

    ax_m.axhline(1.0, color=INK, lw=0.8)
    ax_m.set_ylabel(rf"$M(<{R_CORE:g}\,{{\rm kpc}}) / M_0$")
    ax_m.set_xlabel("Time [Gyr]")

    for ax in (ax_e, ax_m):
        ax.legend(fontsize=8, frameon=False, labelcolor=INK)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
