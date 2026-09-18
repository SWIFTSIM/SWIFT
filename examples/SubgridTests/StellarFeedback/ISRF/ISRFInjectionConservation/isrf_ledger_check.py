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
The ISRF energy ledger: the only valid running reference for the LW/FUV
propagation's carried band energy. It does not depend on any particular
`c_hyp` scheme (uniform pin, per-particle, kernel-local, ...): it is built
from `radiation_end_force_propagation`'s own exact-relaxation update
(`src/feedback/GEAR/radiation_isrf.c`), so it is valid for any arm that
runs that same update, and is meant to be re-run unchanged across
candidates.

Reads, per band (`FUV`, `LW`) and per snapshot:

    E   = sum(mass * <band>SpecificEnergies)
    Inj = sum(mass * <band>CumulativeInjectedSpecificEnergies)
    Abs = sum(mass * <band>CumulativeAbsorbedSpecificEnergies)

Inj and Abs are running totals since the particle's first init (see the
two fields' own doxygen, `src/feedback/GEAR_thermal/feedback_struct.h`),
so this is a global conservation check at each snapshot independently, not
a per-interval one. The gate is

    |E + Abs - Inj| / |Inj|  <=  --tol   (default 1e-3)

which isolates exactly the transport-and-dissipation residual that
`radiation_end_force_propagation`'s I/A split does not attribute to
either injection or the surviving field: the SPH divergence's kernel-sum
identity and the artificial-dissipation term's pairwise antisymmetry
drive this to ~0 when summed over every particle. A same-bin run (every
particle on one time bin, so every interacting pair is exactly mirrored)
closes it far below the bar; a run with active bin seams is expected to
sit closer to it.

The two cumulative fields only exist under `--enable-debugging-checks`
(zero otherwise); a run built without it makes every Inj/Abs column read
exactly 0, which this script reports as a hard error, not a silent 0/0.

Fails on any non-finite E, Inj or Abs (never lets a NaN/Inf compare false
against the bar and pass silently).
"""

import argparse
import glob
import sys

import h5py
import numpy as np

BANDS = ("FUV", "LW")


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to check (default: %(default)s). "
        "Every match is checked (a glob, not just the last one: the ledger "
        "is a running total, so it must close at every snapshot, not only "
        "at the end).",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-3,
        help="Max allowed |E + Abs - Inj| / |Inj| per snapshot per band "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--inj-floor",
        type=float,
        default=1e-30,
        help="Below this |Inj| (internal units), a snapshot is treated as "
        "pre-injection rather than divided by a near-zero denominator: "
        "PASS iff |E + Abs| is also below this floor, FAIL otherwise "
        "(default: %(default)s).",
    )
    return parser.parse_args()


def check_snapshot(path, tol, inj_floor):
    """Return True iff every band at this snapshot passes the ledger gate."""
    all_ok = True
    with h5py.File(path, "r") as f:
        time = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        gas = f["/PartType0"]
        mass = gas["Masses"][:].astype(np.float64)

        for band in BANDS:
            # Per-band: one band's failure must never suppress the other's check.
            band_ok = True

            u = gas[f"{band}SpecificEnergies"][:].astype(np.float64)
            inj = gas[f"{band}CumulativeInjectedSpecificEnergies"][:].astype(np.float64)
            absorbed = gas[f"{band}CumulativeAbsorbedSpecificEnergies"][:].astype(
                np.float64
            )

            for name, arr in (("u", u), ("Inj", inj), ("Abs", absorbed)):
                if not np.all(np.isfinite(arr)):
                    n_bad = int(np.sum(~np.isfinite(arr)))
                    print(
                        f"{path} t={time:.6e} {band}: {name} has {n_bad} "
                        "non-finite value(s) -> FAIL"
                    )
                    band_ok = False

            if not band_ok:
                all_ok = False
                continue

            E = float(np.sum(mass * u))
            Inj = float(np.sum(mass * inj))
            Abs = float(np.sum(mass * absorbed))

            if not all(np.isfinite(x) for x in (E, Inj, Abs)):
                print(
                    f"{path} t={time:.6e} {band}: mass-weighted sum is "
                    f"non-finite (E={E}, Inj={Inj}, Abs={Abs}) -> FAIL"
                )
                all_ok = False
                continue

            residual = E + Abs - Inj

            if abs(Inj) < inj_floor:
                # Pre-injection (or a run built without --enable-debugging-
                # checks, where Inj/Abs read exactly 0): only a genuine
                # all-zero state passes here, never a bare 0/0 divide.
                passed = abs(residual) < inj_floor
                metric = abs(residual)
                note = " (|Inj| below floor: checked |E+Abs-Inj| directly)"
            else:
                metric = abs(residual) / abs(Inj)
                passed = metric <= tol
                note = ""

            status = "PASS" if passed else "FAIL"
            print(
                f"{path} t={time:.6e} {band}: E={E:.10e} Inj={Inj:.10e} "
                f"Abs={Abs:.10e} metric={metric:.6e} (tol={tol:.3e}) "
                f"-> {status}{note}"
            )
            all_ok = all_ok and passed

    return all_ok


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    all_ok = True
    for path in files:
        all_ok = check_snapshot(path, opt.tol, opt.inj_floor) and all_ok

    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
