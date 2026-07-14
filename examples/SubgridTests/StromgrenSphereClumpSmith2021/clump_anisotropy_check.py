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
Compare a nside=0 (spherical) run against a nside=1 (12 HEALPix pixels) run
of this example: bin gas particles by angle from the star to the clump's
pixel direction (+x, see makeIC_clump.py), and report the ionized fraction
and mass inside that pixel versus the other 11 combined, for both runs.

A real redistribution signal looks like: nside=0 shows the "other pixels"
bucket under-served relative to nside=1 (the shared budget got soaked up by
the clump's pixel); nside=1 then shows that corrected. Total ionized mass
may drop slightly at nside=1 -- hard pixel assignment can strand some budget
in an over-provisioned pixel; that alone isn't a failure.

Requires both runs built with -DIONIZATION_FEEDBACK_DEBUG_NO_COOLING (see
README) so IsIonizedFlags accumulates monotonically instead of resetting
every step.

Usage:
    python3 clump_anisotropy_check.py \
        --nside0-snapshot nside0/snap/snapshot_0020.hdf5 \
        --nside1-snapshot nside1/snap/snapshot_0020.hdf5
"""
import argparse

import numpy as np

try:
    import swiftsimio as sw
except ImportError:
    sw = None


def load_buckets(snapshot_path, clump_direction, pixel_halfangle_deg):
    data = sw.load(snapshot_path)

    star_pos = np.array(data.stars.coordinates.to("kpc").value)[0]
    ionized = np.array(data.gas.is_ionized_flags, dtype=bool)
    pos = np.array(data.gas.coordinates.to("kpc").value)
    mass = np.array(data.gas.masses.to("Msun").value)
    box = float(data.metadata.boxsize.to("kpc").value[0])

    disp = pos - star_pos
    disp -= box * np.round(disp / box)  # periodic wrap
    r = np.sqrt(np.sum(disp**2, axis=1))
    nonzero = r > 0
    unit = np.zeros_like(disp)
    unit[nonzero] = disp[nonzero] / r[nonzero, None]

    cos_angle = unit @ clump_direction
    cos_halfangle = np.cos(np.radians(pixel_halfangle_deg))
    in_clump_pixel = cos_angle >= cos_halfangle

    def stats(mask):
        n_tot = int(np.sum(mask))
        n_ion = int(np.sum(mask & ionized))
        m_ion = float(np.sum(mass[mask & ionized]))
        frac = n_ion / n_tot if n_tot else float("nan")
        return {"n_tot": n_tot, "n_ion": n_ion, "mass_ion": m_ion, "frac": frac}

    return {
        "clump_pixel": stats(in_clump_pixel),
        "other_pixels": stats(~in_clump_pixel),
    }


def print_bucket_table(label, buckets):
    print(f"\n--- {label} ---")
    for name in ("clump_pixel", "other_pixels"):
        b = buckets[name]
        print(
            f"  {name:14s}: {b['n_ion']:7d}/{b['n_tot']:7d} ionized "
            f"({100 * b['frac']:6.2f}%), ionized mass = {b['mass_ion']:.4g} Msun"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nside0-snapshot", required=True)
    parser.add_argument("--nside1-snapshot", required=True)
    parser.add_argument(
        "--clump-direction",
        type=float,
        nargs=3,
        default=[1.0, 0.0, 0.0],
        help="Unit direction from the star to the clump (default: +x, "
             "matching makeIC_clump.py's placement at HEALPix nside=1 "
             "pixel 4's center).",
    )
    parser.add_argument(
        "--pixel-halfangle",
        type=float,
        default=29.0,
        help="Half-angle (degrees) defining the 'clump pixel' bucket. "
             "Default 29 matches nside=1 pixel 4's inscribed-cap radius "
             "around its center direction.",
    )
    args = parser.parse_args()

    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    clump_dir = np.array(args.clump_direction)
    clump_dir = clump_dir / np.linalg.norm(clump_dir)

    b0 = load_buckets(args.nside0_snapshot, clump_dir, args.pixel_halfangle)
    b1 = load_buckets(args.nside1_snapshot, clump_dir, args.pixel_halfangle)

    print_bucket_table("nside=0 (spherical, shared budget)", b0)
    print_bucket_table("nside=1 (12 HEALPix pixels)", b1)

    other_frac_0 = b0["other_pixels"]["frac"]
    other_frac_1 = b1["other_pixels"]["frac"]
    clump_frac_0 = b0["clump_pixel"]["frac"]
    clump_frac_1 = b1["clump_pixel"]["frac"]

    total_mass_0 = b0["clump_pixel"]["mass_ion"] + b0["other_pixels"]["mass_ion"]
    total_mass_1 = b1["clump_pixel"]["mass_ion"] + b1["other_pixels"]["mass_ion"]

    print("\n--- Verdict ---")
    print(
        f"'other pixels' ionized fraction: nside=0={100*other_frac_0:.2f}%  "
        f"nside=1={100*other_frac_1:.2f}%  "
        f"(delta = {100*(other_frac_1 - other_frac_0):+.2f} pp)"
    )
    print(
        f"'clump pixel'  ionized fraction: nside=0={100*clump_frac_0:.2f}%  "
        f"nside=1={100*clump_frac_1:.2f}%  "
        f"(delta = {100*(clump_frac_1 - clump_frac_0):+.2f} pp)"
    )
    print(
        f"Total ionized mass             : nside=0={total_mass_0:.4g} Msun  "
        f"nside=1={total_mass_1:.4g} Msun"
    )

    redistribution_detected = (other_frac_1 > other_frac_0) and (
        clump_frac_1 <= clump_frac_0 + 1e-6
    )
    if redistribution_detected:
        print(
            "\nPASS: nside=1 ionizes a larger fraction of the diffuse "
            "gas outside the clump's pixel than nside=0, without the "
            "clump's own pixel doing better -- the angular split is "
            "redistributing the budget as intended."
        )
    else:
        print(
            "\nFAIL (or inconclusive): no clear redistribution signal. Try "
            "increasing makeIC_clump.py's --density_factor or "
            "--clump_radius_pc, or check the clump is within the front's "
            "reach (see README)."
        )


if __name__ == "__main__":
    main()
