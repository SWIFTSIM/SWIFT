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


def counting_sigma(bucket):
    """Binomial standard error on a bucket's ionized fraction.

    A LOWER BOUND on the real uncertainty only: the ionized field is
    spatially correlated (neighbouring particles are ionized together), so
    the effective number of independent samples is well below n_tot. It is
    used here as a floor a claimed effect must clear, never as the true
    error bar.
    """
    n_tot = bucket["n_tot"]
    if n_tot == 0:
        return float("nan")
    f = bucket["frac"]
    return np.sqrt(max(f * (1.0 - f), 0.0) / n_tot)


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
    parser.add_argument(
        "--min-effect-pp",
        type=float,
        default=1.0,
        help="Smallest 'other pixels' ionized-fraction gain (percentage "
        "points) that counts as redistribution (default: 1.0). Without a "
        "minimum effect size the verdict is a sign test between two "
        "independent simulations, which coin-flips to PASS about half the "
        "time when the true effect is zero.",
    )
    parser.add_argument(
        "--n-sigma",
        type=float,
        default=3.0,
        help="The gain must also exceed this many combined counting sigmas "
        "(default: 3.0). See counting_sigma(): that sigma is a lower bound "
        "on the real noise, so this is a floor, not a confidence level.",
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

    # The two runs are independent simulations, so a bare sign comparison
    # (other_frac_1 > other_frac_0) passes about half the time when the true
    # effect is zero. The gain has to clear both an explicit minimum effect
    # size and a counting-noise floor, and the clump pixel is allowed to move
    # up by the same noise before it counts as doing better.
    sigma_other = np.hypot(
        counting_sigma(b0["other_pixels"]), counting_sigma(b1["other_pixels"])
    )
    sigma_clump = np.hypot(
        counting_sigma(b0["clump_pixel"]), counting_sigma(b1["clump_pixel"])
    )
    other_gain = other_frac_1 - other_frac_0
    clump_gain = clump_frac_1 - clump_frac_0
    other_threshold = max(0.01 * args.min_effect_pp, args.n_sigma * sigma_other)
    clump_threshold = args.n_sigma * sigma_clump

    print(
        f"Redistribution threshold       : 'other pixels' gain must exceed "
        f"{100 * other_threshold:.2f} pp "
        f"(max of {args.min_effect_pp:.2f} pp minimum effect and "
        f"{args.n_sigma:.1f} x {100 * sigma_other:.3f} pp counting sigma, "
        f"itself a lower bound on the noise), while the 'clump pixel' gain "
        f"stays below {100 * clump_threshold:.2f} pp"
    )

    redistribution_detected = (other_gain > other_threshold) and (
        clump_gain <= clump_threshold
    )
    if redistribution_detected:
        print(
            f"\nPASS: 'other pixels' gained {100 * other_gain:+.2f} pp "
            f"(> {100 * other_threshold:.2f} pp) while the clump's own pixel "
            f"gained {100 * clump_gain:+.2f} pp (<= "
            f"{100 * clump_threshold:.2f} pp) -- nside=1 ionizes more of the "
            f"diffuse gas outside the clump's pixel without the clump's pixel "
            f"doing better, so the angular split is redistributing the budget "
            f"as intended."
        )
    else:
        print(
            f"\nFAIL (or inconclusive): 'other pixels' gained "
            f"{100 * other_gain:+.2f} pp against a "
            f"{100 * other_threshold:.2f} pp threshold, clump pixel gained "
            f"{100 * clump_gain:+.2f} pp against {100 * clump_threshold:.2f} "
            f"pp. No redistribution signal above the noise. Try increasing "
            f"makeIC_clump.py's --density_factor or --clump_radius_pc, or "
            f"check the clump is within the front's reach (see README)."
        )


if __name__ == "__main__":
    main()
