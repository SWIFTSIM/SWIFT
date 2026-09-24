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
Within-run bitwise identity check for the Lyman-Werner-band photon-number
moment (`ISRF_MOMENT_LW_PHOTON`): at H = 0, non-cosmological, it is injected
by direct assignment from the LW energy moment and shares its operator, so
the two must read bitwise identical on every snapshot.

Three independent field pairs, every particle, every snapshot, compared on
raw bit patterns (not float `==`): `LWSpecificEnergies` vs
`LWPhotonSpecificEnergies`, `LWSpecificFluxes` vs `LWPhotonSpecificFluxes`,
`LWSpecificFluxDivergences` vs `LWPhotonSpecificFluxDivergences`.
`LWArtificialDissipationCoefficients` vs
`LWPhotonArtificialDissipationCoefficients` is also compared but is
decorative: both getters index the same shared operator field by a literal
constant, so this pair cannot fail for any bug this check's sibling fields
would not already catch. The `grad_u` accumulator has no snapshot field in
any build; its own identity is checked in
`tests/testRadiationISRFGradientCache.c`, not here.

Every array is tested for finiteness before the bitwise comparison: a NaN
compares False against every equality test, so an undetected one would
silently pass rather than fail.
"""

import argparse
import glob
import sys

import h5py
import numpy as np

FIELD_PAIRS = [
    ("LWSpecificEnergies", "LWPhotonSpecificEnergies"),
    ("LWSpecificFluxes", "LWPhotonSpecificFluxes"),
    ("LWSpecificFluxDivergences", "LWPhotonSpecificFluxDivergences"),
]
DECORATIVE_PAIR = (
    "LWArtificialDissipationCoefficients",
    "LWPhotonArtificialDissipationCoefficients",
)


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s)",
    )
    return parser.parse_args()


def bits_equal(a, b):
    """Compare two float32 arrays on their raw bit patterns."""
    return np.array_equal(a.astype(np.float32).view(np.uint32),
                          b.astype(np.float32).view(np.uint32))


def check_snapshot(path, gate_pairs):
    """Run every field pair's bitwise check on one snapshot file.

    `gate_pairs` fail the check; `DECORATIVE_PAIR` is only reported.
    """
    ok = True
    with h5py.File(path, "r") as f:
        gas = f["PartType0"]
        for a_name, b_name in gate_pairs + [DECORATIVE_PAIR]:
            a = gas[a_name][()]
            b = gas[b_name][()]
            gated = (a_name, b_name) != DECORATIVE_PAIR
            if not (np.isfinite(a).all() and np.isfinite(b).all()):
                print(f"FAIL {path}: non-finite value in {a_name} or {b_name}")
                ok = ok and not gated
                continue
            if not bits_equal(a, b):
                n = a.astype(np.float32).view(np.uint32)
                m = b.astype(np.float32).view(np.uint32)
                n_diff = int(np.sum(n != m))
                print(
                    f"FAIL {path}: {a_name} vs {b_name} not bitwise identical "
                    f"({n_diff}/{n.size} elements differ)"
                )
                ok = ok and not gated
            elif gated and not np.any(a) and not np.any(b):
                print(f"NOTE {path}: {a_name}/{b_name} both all-zero this snapshot")
    return ok


def main():
    opt = parse_options()
    snapshots = sorted(glob.glob(opt.snapshot))
    if not snapshots:
        sys.exit(f"No snapshots matched {opt.snapshot!r}")

    all_ok = True
    any_nonzero = False
    for path in snapshots:
        with h5py.File(path, "r") as f:
            gas = f["PartType0"]
            if np.any(gas["LWSpecificEnergies"][()]):
                any_nonzero = True
        if not check_snapshot(path, FIELD_PAIRS):
            all_ok = False

    if not any_nonzero:
        sys.exit(
            "FAIL: LWSpecificEnergies is all-zero in every snapshot -- this "
            "fixture never illuminated any gas, so the identity check above "
            "is vacuous (both sides zero proves nothing)."
        )

    if not all_ok:
        sys.exit("RESULT: FAIL")
    print(f"RESULT: PASS ({len(snapshots)} snapshots, {len(FIELD_PAIRS)} gated "
          "field pairs, bitwise)")


if __name__ == "__main__":
    main()
