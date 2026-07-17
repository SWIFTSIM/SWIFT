################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
################################################################################
"""
ID-tracked radial density profile of the clump: bin the original t=0
clump-region particles by their t=0 distance from the clump centre, and
plot the median n_H of those same particles at t=0 against the final
snapshot, for both schemes. Complements clump_erosion_matrix.py and
clump_hemisphere_erosion.py -- a radially/spatially aggregated diagnostic
like this one cannot resolve a directional, near-hemisphere-only effect
(see clump_hemisphere_erosion.py and the README for why).

Usage:
    python3 clump_radial_profile.py
    python3 clump_radial_profile.py --gas-mass 20 --metallicity z0231
"""

import argparse
import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import swiftsimio as sw

UnitLength_in_cgs = 3.0856775814913673e21  # kpc
PC_IN_CODE = (3.0856775814913673e18) / UnitLength_in_cgs

CLUMP_DISTANCE_PC = 20.0
CLUMP_RADIUS_PC = 10.0
DENSITY_THRESHOLD_NH = 1000.0  # cm^-3, the clump's own IC density
M_P_CGS = 1.67262192369e-24
N_BINS = 20
R_MAX_PC = 14.0


def load_gas(path: str) -> tuple:
    """Load gas IDs, positions, and n_H from a snapshot."""
    data = sw.load(path)
    pos = data.gas.coordinates.to("kpc").value
    ids = data.gas.particle_ids.value
    rho = data.gas.densities.to("g/cm**3").value
    boxsize = data.metadata.boxsize.to("kpc").value
    star_pos = data.stars.coordinates.to("kpc").value[0]
    n_H = rho / M_P_CGS
    return ids, pos, n_H, boxsize, star_pos


def profile_by_id(run_dir: str) -> tuple:
    """Median n_H at t=0 and at the final snapshot, binned by t=0 radius
    from the clump centre, for the same fixed set of particle IDs."""
    ids0, pos0, nH0, box0, star0 = load_gas(f"{run_dir}/snap/snapshot_0000.hdf5")
    clump_center = star0 + np.array([CLUMP_DISTANCE_PC * PC_IN_CODE, 0.0, 0.0])
    dx0 = pos0 - clump_center
    dx0 -= box0 * np.round(dx0 / box0)
    r0_pc = np.sqrt((dx0**2).sum(axis=1)) / PC_IN_CODE
    in_range = r0_pc < R_MAX_PC

    ids0_sel = ids0[in_range]
    r0_sel = r0_pc[in_range]
    nH0_sel = nH0[in_range]

    snap8 = sorted(glob.glob(f"{run_dir}/snap/snapshot_*.hdf5"))[-1]
    ids8, pos8, nH8, box8, star8 = load_gas(snap8)
    id_to_idx8 = {pid: i for i, pid in enumerate(ids8)}
    idx8 = np.array([id_to_idx8.get(pid, -1) for pid in ids0_sel])
    found = idx8 >= 0
    nH8_sel = np.full(len(ids0_sel), np.nan)
    nH8_sel[found] = nH8[idx8[found]]

    bins = np.linspace(0, R_MAX_PC, N_BINS + 1)
    bin_idx = np.digitize(r0_sel, bins) - 1
    med0 = np.full(N_BINS, np.nan)
    med8 = np.full(N_BINS, np.nan)
    for b in range(N_BINS):
        m = bin_idx == b
        if m.sum() > 0:
            med0[b] = np.median(nH0_sel[m])
            med8[b] = np.nanmedian(nH8_sel[m])
    centers = 0.5 * (bins[:-1] + bins[1:])
    return centers, med0, med8


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gas-mass", nargs="+", default=["20", "95"], help="gas_mass labels (default: 20 95)."
    )
    parser.add_argument(
        "--metallicity",
        nargs="+",
        default=["z0231", "z0"],
        help="Metallicity labels (default: z0231 z0).",
    )
    parser.add_argument(
        "--rebuild", type=str, default="0p01", help="Rebuild-cadence label (default: 0p01)."
    )
    parser.add_argument("-o", "--output", type=str, default="clump_radial_profile_matrix.png")
    args = parser.parse_args()

    panels = [
        (f"gas_mass={g}\n{z}", f"gas{g}_{z}") for g in args.gas_mass for z in args.metallicity
    ]

    fig, axes = plt.subplots(1, len(panels), figsize=(4 * len(panels), 4), sharey=True)
    if len(panels) == 1:
        axes = [axes]
    for ax, (title, prefix) in zip(axes, panels):
        c_u, m0_u, m8_u = profile_by_id(f"{prefix}_nside0_rebuild{args.rebuild}")
        c_a, m0_a, m8_a = profile_by_id(f"{prefix}_nside1_rebuild{args.rebuild}")
        ax.plot(c_u, m0_u, "k--", label="t=0 (IC)")
        ax.plot(c_u, m8_u, color="tab:orange", label="uniform, final")
        ax.plot(c_a, m8_a, color="tab:green", label="angular, final")
        ax.axvline(CLUMP_RADIUS_PC, color="grey", linewidth=0.7, linestyle=":")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("$r$ from clump centre [pc] (t=0 position)")
    axes[0].set_ylabel(r"median $n_H$ [cm$^{-3}$]")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle(
        "ID-tracked: binned by t=0 radius, density measured on the SAME particles at final time"
    )
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
