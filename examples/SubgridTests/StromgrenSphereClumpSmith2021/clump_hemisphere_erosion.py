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
Split the ID-tracked clump-erosion measurement (see
clump_erosion_matrix.py) into the star-facing (near) and opposite (far)
hemisphere, split by the star->clump axis. Mass-biasing
(Section "Angular (HealPix) budget splitting" in the theory logbook) is
specifically a directional failure: an isotropic cause (bulk heating,
clump settling) would erode both hemispheres equally, a directional cause
would not. Also reports total mass still present (dense or not) as a
mass-conservation sanity check: this must read 100% for every
configuration, since the same fixed set of particle IDs is summed on both
sides -- any deviation would mean particles were lost, not eroded.

Usage:
    python3 clump_hemisphere_erosion.py
    python3 clump_hemisphere_erosion.py --gas-mass 20 --metallicity z0 --rebuild 0p01
"""

import argparse

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


def run_case(snap0: str, snap8: str) -> dict:
    """Return near/far/total still-dense and still-present percentages."""
    ids0, pos0, mass0, nH0, box0, star0 = load_gas(snap0)
    clump_center = star0 + np.array([CLUMP_DISTANCE_PC * PC_IN_CODE, 0.0, 0.0])

    # Membership is defined ONCE, at t=0, by position + density -- this is
    # the only time a spatial selection is unambiguous (nothing has had
    # time to drift or be swept in yet).
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

    clump_ids = ids0[in_clump0]
    near_ids = ids0[near0]
    far_ids = ids0[far0]
    mass0_total = mass0[in_clump0].sum()
    mass0_near = mass0[near0].sum()
    mass0_far = mass0[far0].sum()

    # At t=8 Myr, look the SAME IDs up wherever they now are, and ask
    # whether they are still above the clump's own IC density threshold --
    # no re-selection by position, so drift and any background gas swept
    # into the old spatial footprint cannot contaminate the count.
    ids8, pos8, mass8, nH8, box8, star8 = load_gas(snap8)
    id_to_idx8 = {pid: i for i, pid in enumerate(ids8)}

    def surviving_mass(id_list: np.array) -> tuple:
        idx = np.array([id_to_idx8[pid] for pid in id_list if pid in id_to_idx8])
        if len(idx) == 0:
            return 0.0, 0.0
        still_dense = nH8[idx] > DENSITY_THRESHOLD_NH
        return mass8[idx][still_dense].sum(), mass8[idx].sum()

    total_dense8, total_all8 = surviving_mass(clump_ids)
    near_dense8, near_all8 = surviving_mass(near_ids)
    far_dense8, far_all8 = surviving_mass(far_ids)

    return {
        "n_clump": len(clump_ids),
        "n_near": len(near_ids),
        "n_far": len(far_ids),
        "total": 100 * total_dense8 / mass0_total,
        "near": 100 * near_dense8 / mass0_near,
        "far": 100 * far_dense8 / mass0_far,
        "total_present": 100 * total_all8 / mass0_total,
        "near_present": 100 * near_all8 / mass0_near,
        "far_present": 100 * far_all8 / mass0_far,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gas-mass", nargs="+", default=["20", "95"], help="gas_mass labels (default: 20 95)."
    )
    parser.add_argument(
        "--metallicity", nargs="+", default=["z0"], help="Metallicity labels (default: z0)."
    )
    parser.add_argument(
        "--rebuild", nargs="+", default=["0p01"], help="Rebuild-cadence labels (default: 0p01)."
    )
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("-o", "--output", type=str, default="clump_hemisphere_erosion.png")
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
                    r = run_case(
                        f"{run_dir}/snap/snapshot_0000.hdf5",
                        sorted(
                            __import__("glob").glob(f"{run_dir}/snap/snapshot_*.hdf5")
                        )[-1],
                    )
                    plot_data[panel_key][scheme] = r
                    print(f"{run_dir}")
                    print(
                        f"  N_clump(t=0)={r['n_clump']}  N_near={r['n_near']}  "
                        f"N_far={r['n_far']}"
                    )
                    print(
                        f"  still-dense (n_H>{DENSITY_THRESHOLD_NH:.0f}): "
                        f"total={r['total']:5.1f}%  near={r['near']:5.1f}%  "
                        f"far={r['far']:5.1f}%"
                    )
                    print(
                        f"  still-present at all (mass conservation check): "
                        f"total={r['total_present']:5.1f}%  "
                        f"near={r['near_present']:5.1f}%  "
                        f"far={r['far_present']:5.1f}%"
                    )
                    print()

    if args.no_plot:
        return

    categories = ["total", "near\n(star-facing)", "far"]
    keys = ["total", "near", "far"]
    fig, axes = plt.subplots(1, len(plot_data), figsize=(4.5 * len(plot_data), 4.2), sharey=True)
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
        ax.set_ylim(90, 101)
        for xi, v in zip(x - width / 2, uni_vals):
            ax.text(xi, v + 0.15, f"{v:.1f}%", ha="center", fontsize=8)
        for xi, v in zip(x + width / 2, ang_vals):
            ax.text(xi, v + 0.55, f"{v:.1f}%", ha="center", fontsize=8)

    axes[0].set_ylabel("Clump mass remaining at $t=8$ Myr (%, ID-tracked)")
    axes[0].legend(loc="lower right")
    fig.suptitle("ID-tracked hemisphere-split clump erosion")
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
