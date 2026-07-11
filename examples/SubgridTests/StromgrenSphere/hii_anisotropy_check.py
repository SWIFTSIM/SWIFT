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
Check whether the HII ionization search is finding gas isotropically around
the star, or whether some directions are being systematically missed (e.g.
because a star sits near/at a top-level cell corner and the corner-direction
task splitting/linking is broken for that direction).

Methodology: within a radius safely inside the box (well below the box
half-width, to avoid periodic self-interaction contaminating the result) and
past at least one top-level cell width (so the check actually exercises
cell-boundary crossing), count what fraction of gas particles are NOT tagged
as ionized. Report this both overall and broken down by octant (sign of each
coordinate relative to the star) -- a real geometric/task-graph reach bug
shows up as a strong asymmetry between octants or axes; ordinary
photon-budget lag near the growing edge shows up as a small, isotropic
deficit.

Requires a run built with -DIONIZATION_FEEDBACK_DEBUG_NO_COOLING (so ionized
particles stay tagged forever and the deficit isn't confounded by cooling/
re-ionization) -- see the configure line in the project's memory/plan notes.

Usage:
    python3 hii_anisotropy_check.py [-s snap/snapshot_*.hdf5] [-r 0.02]
"""
import argparse
import glob

import numpy as np

try:
    import swiftsimio as sw
except ImportError:
    sw = None


def find_peak_hii_and_star_position(snapshot_glob):
    """Scan all snapshots for the largest HIIRegionRadii and the star
    position at that moment (h_hii resets to 0 once a star dies, so the
    final snapshot alone isn't enough)."""
    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    peak_r_hii = 0.0
    star_pos = None
    for f in files:
        data = sw.load(f)
        r = float(np.max(data.stars.hiiregion_radii).to("kpc").value)
        if r > peak_r_hii:
            peak_r_hii = r
            star_pos = np.array(
                data.stars.coordinates.to("kpc").value
            )[np.argmax(data.stars.hiiregion_radii)]

    return peak_r_hii, star_pos, files[-1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot-glob",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to read (default: snap/snapshot_*.hdf5)",
    )
    parser.add_argument(
        "-r",
        "--radius",
        type=float,
        default=None,
        help="Analysis radius in kpc. Default: min(peak r_hii, 0.4 * boxsize) "
        "to stay clear of periodic self-interaction.",
    )
    args = parser.parse_args()

    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    peak_r_hii, star_pos, last_snapshot = find_peak_hii_and_star_position(
        args.snapshot_glob
    )

    data = sw.load(last_snapshot)
    box = float(data.metadata.boxsize.to("kpc").value[0])

    radius = args.radius
    if radius is None:
        radius = min(peak_r_hii, 0.4 * box)

    print(f"Star position       : {star_pos} kpc")
    print(f"Peak r_hii reached   : {peak_r_hii:.6f} kpc")
    print(f"Box size             : {box:.6f} kpc (half-width {0.5*box:.6f} kpc)")
    print(f"Analysis radius      : {radius:.6f} kpc")
    if radius > 0.4 * box:
        print(
            "  WARNING: analysis radius exceeds 0.4x the box size -- periodic "
            "self-interaction may contaminate the result. Pass --radius "
            "explicitly to shrink it."
        )

    ionized = np.array(data.gas.is_ionized_flags, dtype=bool)
    pos = np.array(data.gas.coordinates.to("kpc").value)

    disp = pos - star_pos
    disp -= box * np.round(disp / box)  # periodic wrap
    r_dist = np.sqrt(np.sum(disp**2, axis=1))
    within = r_dist < radius

    n_within = int(np.sum(within))
    if n_within == 0:
        raise RuntimeError("No gas particles found within the analysis radius.")

    not_ionized_within = within & ~ionized
    n_missing = int(np.sum(not_ionized_within))
    print(
        f"\nOverall: {n_missing}/{n_within} particles within {radius:.4f} kpc "
        f"not ionized ({100 * n_missing / n_within:.2f}%)"
    )

    # Octant breakdown: sign of each coordinate relative to the star.
    print("\nBy octant (sign of displacement from star along each axis):")
    octant_signs = (disp[within] > 0).astype(int)
    octant_id = octant_signs[:, 0] * 4 + octant_signs[:, 1] * 2 + octant_signs[:, 2]
    not_ion_within = ~ionized[within]
    for oct_i in range(8):
        mask = octant_id == oct_i
        n_tot = int(np.sum(mask))
        n_miss = int(np.sum(mask & not_ion_within))
        signs = [("-", "+")[b] for b in [(oct_i >> 2) & 1, (oct_i >> 1) & 1, oct_i & 1]]
        pct = 100 * n_miss / n_tot if n_tot else float("nan")
        print(f"  ({signs[0]}x,{signs[1]}y,{signs[2]}z): {n_miss}/{n_tot} missing ({pct:.2f}%)")

    # Per-axis breakdown (collapses the other two axes) -- highlights a
    # single-axis asymmetry (e.g. +x vs -x) more directly than the octant
    # table above.
    print("\nBy single-axis sign (collapsing the other two axes):")
    for axis, name in enumerate(["x", "y", "z"]):
        for sign, label in [(False, "-"), (True, "+")]:
            mask = (disp[within][:, axis] > 0) == sign
            n_tot = int(np.sum(mask))
            n_miss = int(np.sum(mask & not_ion_within))
            pct = 100 * n_miss / n_tot if n_tot else float("nan")
            print(f"  {label}{name}: {n_miss}/{n_tot} missing ({pct:.2f}%)")

    # Directional anisotropy relative to the (1,1,1) body diagonal, for the
    # un-ionized particles specifically -- a nonzero mean indicates the
    # missing gas is clustered toward/away from that direction rather than
    # scattered randomly.
    if n_missing > 0:
        d_missing = disp[not_ionized_within]
        r_missing = np.sqrt(np.sum(d_missing**2, axis=1))
        unit = d_missing / r_missing[:, None]
        diag = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        cos_diag = unit @ diag
        print(
            f"\nMean cos(angle to (1,1,1) diagonal) of missing particles: "
            f"{np.mean(cos_diag):.4f} (isotropic expectation ~0, std ~0.577)"
        )
        print(f"std: {np.std(cos_diag):.4f}")


if __name__ == "__main__":
    main()
