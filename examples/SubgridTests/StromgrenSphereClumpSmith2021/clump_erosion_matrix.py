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
ID-tracked clump mass-remaining measurement across the validation matrix
(gas_mass x metallicity x scheme x rebuild cadence). Membership is fixed
once, by particle ID, at t=0 (position + density threshold); at the final
snapshot the same IDs are looked up wherever they have since moved to, and
the fraction of their t=0 mass still above the clump's own IC density
threshold is reported. Position-based re-selection at t=8 Myr (an earlier
version of this measurement) is not used -- it implicitly assumes the
clump neither drifts nor has non-clump gas swept into its old footprint,
which produced internally inconsistent results (see README).

Usage:
    python3 clump_erosion_matrix.py
    python3 clump_erosion_matrix.py --gas-mass 20 --metallicity z0231
    python3 clump_erosion_matrix.py --no-plot
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
    """Load gas IDs, positions, masses, and n_H from a snapshot."""
    data = sw.load(path)
    pos = data.gas.coordinates.to("kpc").value
    mass = data.gas.masses.to("Msun").value
    rho = data.gas.densities.to("g/cm**3").value
    ids = data.gas.particle_ids.value
    boxsize = data.metadata.boxsize.to("kpc").value
    star_pos = data.stars.coordinates.to("kpc").value[0]
    n_H = rho / M_P_CGS
    return ids, pos, mass, n_H, boxsize, star_pos


def clump_ids_and_mass_t0(snap0: str) -> dict:
    """Identify clump-member particle IDs and their t=0 mass, once."""
    ids0, pos0, mass0, nH0, box0, star0 = load_gas(snap0)
    clump_center = star0 + np.array([CLUMP_DISTANCE_PC * PC_IN_CODE, 0.0, 0.0])
    dx0 = pos0 - clump_center
    dx0 -= box0 * np.round(dx0 / box0)
    r0 = np.sqrt((dx0**2).sum(axis=1))
    in_clump0 = (r0 < CLUMP_RADIUS_PC * PC_IN_CODE) & (nH0 > DENSITY_THRESHOLD_NH)

    return {
        "clump_ids": ids0[in_clump0],
        "mass0_total": mass0[in_clump0].sum(),
    }


def still_dense_fraction(
    id_list: np.array, ids8: np.array, nH8: np.array, mass8: np.array
) -> float:
    """Sum the t=8 Myr mass of the given IDs that is still above the
    clump's own density threshold; 0.0 for IDs no longer found at all."""
    id_to_idx8 = {pid: i for i, pid in enumerate(ids8)}
    idx = np.array([id_to_idx8[pid] for pid in id_list if pid in id_to_idx8])
    if len(idx) == 0:
        return 0.0
    still_dense = nH8[idx] > DENSITY_THRESHOLD_NH
    return mass8[idx][still_dense].sum()


def last_snapshot(run_dir: str) -> str:
    """Return the path of the highest-index snapshot in run_dir/snap/."""
    files = sorted(glob.glob(f"{run_dir}/snap/snapshot_*.hdf5"))
    if not files:
        raise FileNotFoundError(f"No snapshots found under {run_dir}/snap/")
    return files[-1]


def run_case(run_dir: str) -> float:
    """Return the % of the clump's t=0 mass still dense at the final
    snapshot of run_dir."""
    t0 = clump_ids_and_mass_t0(f"{run_dir}/snap/snapshot_0000.hdf5")
    ids8, pos8, mass8, nH8, box8, star8 = load_gas(last_snapshot(run_dir))
    dense8_total = still_dense_fraction(t0["clump_ids"], ids8, nH8, mass8)
    return 100 * dense8_total / t0["mass0_total"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gas-mass",
        nargs="+",
        default=["20", "95"],
        help="gas_mass labels to sweep (matches run_name's gas${gas_mass} "
        "prefix, default: 20 95).",
    )
    parser.add_argument(
        "--metallicity",
        nargs="+",
        default=["z0231", "z0"],
        help="Metallicity labels to sweep (matches run_name's z${z_label} "
        "prefix, default: z0231 z0).",
    )
    parser.add_argument(
        "--rebuild",
        nargs="+",
        default=["0p01", "0p05", "0p1"],
        help="Rebuild-cadence labels to sweep (default: 0p01 0p05 0p1).",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Only print the table, skip the bar-chart figure.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="clump_erosion_matrix.png",
        help="Figure output path (default: clump_erosion_matrix.png).",
    )
    args = parser.parse_args()

    schemes = [("uniform", "nside0"), ("angular", "nside1")]
    results = {}

    print("=" * 100)
    print("Clump mass remaining at t=8 Myr (% of t=0), ID-tracked, still n_H>1000")
    print("=" * 100)
    for gas_label in args.gas_mass:
        for z_label in args.metallicity:
            panel_key = f"gas_mass={gas_label}, {z_label}"
            results[panel_key] = {}
            for scheme, nside in schemes:
                row = []
                for rb in args.rebuild:
                    run_dir = f"gas{gas_label}_{z_label}_{nside}_rebuild{rb}"
                    row.append(run_case(run_dir))
                results[panel_key][scheme] = row
                print(
                    f"gas_mass={gas_label:5s} {z_label:6s} {scheme:8s} : "
                    + "  ".join(f"{v:5.1f}%" for v in row)
                )
            print()

    if args.no_plot:
        return

    fig, axes = plt.subplots(1, len(results), figsize=(4 * len(results), 4), sharey=True)
    if len(results) == 1:
        axes = [axes]
    x = np.arange(len(args.rebuild))
    width = 0.35
    for ax, (title, schemes_data) in zip(axes, results.items()):
        ax.bar(
            x - width / 2, schemes_data["uniform"], width, label="uniform", color="tab:orange"
        )
        ax.bar(
            x + width / 2, schemes_data["angular"], width, label="angular", color="tab:green"
        )
        ax.axhline(100, color="black", linewidth=1, linestyle=":")
        ax.set_xticks(x)
        ax.set_xticklabels([rb.replace("p", ".") + "\nMyr" for rb in args.rebuild])
        ax.set_title(title, fontsize=10)
        ax.set_ylim(90, 101)
        for xi, v in zip(x - width / 2, schemes_data["uniform"]):
            ax.text(xi, v + 0.15, f"{v:.1f}", ha="center", fontsize=7)
        for xi, v in zip(x + width / 2, schemes_data["angular"]):
            ax.text(xi, v + 0.15, f"{v:.1f}", ha="center", fontsize=7)

    axes[0].set_ylabel("Clump mass remaining at $t=8$ Myr (%, ID-tracked)")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("ID-tracked clump-erosion matrix")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
