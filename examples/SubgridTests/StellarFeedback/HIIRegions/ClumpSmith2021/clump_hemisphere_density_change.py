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
Continuous companion to clump_hemisphere_erosion.py's binary "still above
the clump's own IC density threshold" metric. Uses the same ID-tracked,
near/far (star-facing/opposite) hemisphere split, but reports the ratio of
each hemisphere's median n_H at the final snapshot to its own median n_H
at t=0, for the same particles. This catches dilution the threshold-based
metric cannot: gas that drops from, e.g., 10,000 to 3,000 cm^-3 counts as
zero erosion under "still-dense" (n_H>1000 cm^-3) but is a real, sizeable
density change. A radially-binned profile (clump_radial_profile.py) cannot
show this either -- it averages over the whole sphere at a given radius,
diluting a near-side-only effect with the far side's near-total absence of
change. Reports median, not mean, and only over particles found in both
snapshots (no re-selection by position), for the same reasons given in
clump_hemisphere_erosion.py.

Usage:
    python3 clump_hemisphere_density_change.py
    python3 clump_hemisphere_density_change.py --gas-mass 20 --metallicity z0 --rebuild 0p01
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


def load_gas(path: str) -> tuple:
    """Load gas IDs, positions, and n_H from a snapshot."""
    data = sw.load(path)
    pos = data.gas.coordinates.to("kpc").value
    rho = data.gas.densities.to("g/cm**3").value
    ids = data.gas.particle_ids.value
    boxsize = data.metadata.boxsize.to("kpc").value
    star_pos = data.stars.coordinates.to("kpc").value[0]
    n_H = rho / M_P_CGS
    return ids, pos, n_H, boxsize, star_pos


def run_case(snap0: str, snap_final: str) -> dict:
    """Return near/far median n_H at t=0 and at the final snapshot, and
    their ratio, for the same ID-tracked clump particles."""
    ids0, pos0, nH0, box0, star0 = load_gas(snap0)
    clump_center = star0 + np.array([CLUMP_DISTANCE_PC * PC_IN_CODE, 0.0, 0.0])

    dx0 = pos0 - clump_center
    dx0 -= box0 * np.round(dx0 / box0)
    r0 = np.sqrt((dx0**2).sum(axis=1))
    in_clump0 = (r0 < CLUMP_RADIUS_PC * PC_IN_CODE) & (nH0 > DENSITY_THRESHOLD_NH)

    star_to_clump = clump_center - star0
    star_to_clump -= box0 * np.round(star_to_clump / box0)
    axis = star_to_clump / np.linalg.norm(star_to_clump)
    d_along_axis0 = dx0 @ axis
    near0 = in_clump0 & (d_along_axis0 < 0)
    far0 = in_clump0 & (d_along_axis0 >= 0)

    ids_final, _, nH_final, _, _ = load_gas(snap_final)
    id_to_idx_final = {pid: i for i, pid in enumerate(ids_final)}

    def median_change(mask0: np.array) -> tuple:
        ids_hemi = ids0[mask0]
        nH0_hemi = nH0[mask0]
        idx_final = np.array(
            [id_to_idx_final.get(pid, -1) for pid in ids_hemi], dtype=int
        )
        found = idx_final >= 0
        n_found = int(found.sum())
        if n_found == 0:
            return 0.0, 0.0, 0
        median_nH0 = np.median(nH0_hemi[found])
        median_nH_final = np.median(nH_final[idx_final[found]])
        return median_nH0, median_nH_final, n_found

    near_nH0, near_nH_final, n_near_found = median_change(near0)
    far_nH0, far_nH_final, n_far_found = median_change(far0)
    total_nH0, total_nH_final, n_total_found = median_change(in_clump0)

    return {
        "n_clump": int(in_clump0.sum()),
        "n_near": int(near0.sum()),
        "n_far": int(far0.sum()),
        "n_near_found": n_near_found,
        "n_far_found": n_far_found,
        "n_total_found": n_total_found,
        "near_ratio": 100 * near_nH_final / near_nH0 if near_nH0 > 0 else 0.0,
        "far_ratio": 100 * far_nH_final / far_nH0 if far_nH0 > 0 else 0.0,
        "total_ratio": 100 * total_nH_final / total_nH0 if total_nH0 > 0 else 0.0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gas-mass", nargs="+", default=["20", "95"], help="gas_mass labels."
    )
    parser.add_argument(
        "--metallicity", nargs="+", default=["z0"], help="Metallicity labels."
    )
    parser.add_argument(
        "--rebuild", nargs="+", default=["0p01"], help="Rebuild-cadence labels."
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument(
        "-o", "--output", type=str, default="clump_hemisphere_density_change.png"
    )
    args = parser.parse_args()

    schemes = [("uniform", "nside0"), ("angular", "nside1")]
    plot_data = {}

    for gas_label in args.gas_mass:
        for z_label in args.metallicity:
            for rb in args.rebuild:
                panel_key = f"gas_mass={gas_label}"
                plot_data.setdefault(panel_key, {})
                for scheme, nside in schemes:
                    run_dir = f"gas{gas_label}_{z_label}_{nside}_rebuild{rb}"
                    snap0 = f"{run_dir}/snap/snapshot_0000.hdf5"
                    snap_final = sorted(glob.glob(f"{run_dir}/snap/snapshot_*.hdf5"))[
                        -1
                    ]
                    r = run_case(snap0, snap_final)
                    plot_data[panel_key][scheme] = r
                    print(f"{run_dir}")
                    print(
                        f"  N_clump(t=0)={r['n_clump']}  N_near={r['n_near']} "
                        f"(found={r['n_near_found']})  N_far={r['n_far']} "
                        f"(found={r['n_far_found']})"
                    )
                    print(
                        f"  median n_H final/t=0: total={r['total_ratio']:5.1f}%  "
                        f"near={r['near_ratio']:5.1f}%  far={r['far_ratio']:5.1f}%"
                    )
                    print()

    if args.no_plot:
        return

    categories = ["total", "near\n(star-facing)", "far"]
    keys = ["total_ratio", "near_ratio", "far_ratio"]
    fig, axes = plt.subplots(
        1, len(plot_data), figsize=(4.5 * len(plot_data), 4.2), sharey=True
    )
    if len(plot_data) == 1:
        axes = [axes]
    for ax, (panel_label, schemes_data) in zip(axes, plot_data.items()):
        x = np.arange(len(keys))
        width = 0.35
        uni_vals = [schemes_data["uniform"][k] for k in keys]
        ang_vals = [schemes_data["angular"][k] for k in keys]
        ax.bar(x - width / 2, uni_vals, width, label="uniform", color="tab:orange")
        ax.bar(x + width / 2, ang_vals, width, label="angular", color="tab:green")
        ax.axhline(100, color="black", linewidth=1, linestyle=":")
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.set_title(panel_label)
        for xi, v in zip(x - width / 2, uni_vals):
            ax.text(xi, v + 0.5, f"{v:.1f}%", ha="center", fontsize=8)
        for xi, v in zip(x + width / 2, ang_vals):
            ax.text(xi, v + 0.5, f"{v:.1f}%", ha="center", fontsize=8)

    axes[0].set_ylabel("Median $n_H$, final / t=0 (%, ID-tracked)")
    axes[0].legend(loc="lower right")
    fig.suptitle("ID-tracked hemisphere-split clump density change (continuous)")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
