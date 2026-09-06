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
Checks the LW/FUV injection kernel-weight normalization directly: at Z=0
the receiver-side dust extinction factor is exp(0)=1 exactly (no dust), so
for one star-feedback pass, summing the injected energy over every gas
particle the star's kernel reaches must reproduce Delta_t * L_band exactly,
up to floating-point precision:

    sum_j(u_FUV_j * m_j) == Delta_t * L_FUV_star
    sum_j(u_LW_j  * m_j) == Delta_t * L_LW_star

This is an exact conservation identity (the injection formula's own
`sum_j weight_j = 1` normalization), not an approximate physics comparison
-- unlike Tier 1's Yukawa-profile decay-length fit, no loose tolerance is
expected here.

Delta_t (the star's own feedback-timestep at the checked snapshot) is read
from swift's own step-table log, not assumed to equal TimeIntegration:dt_max.
"""

import argparse
import glob
import sys

import h5py
import numpy as np


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s). "
        "The last one (highest time) is checked.",
    )
    parser.add_argument(
        "--log",
        default="output.log",
        help="Run log to read the star's own feedback Delta_t from "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-5,
        help="Max allowed relative error (default: %(default)s). This is an "
        "exact identity, so the only expected discrepancy is float32 "
        "snapshot-storage rounding (measured ~1e-8 for a ~50-particle "
        "kernel); the default leaves ample margin above that floor.",
    )
    return parser.parse_args()


def read_delta_t(log_path, time):
    """Read the star's feedback Delta_t from swift's own step-table log, by
    matching the row whose Time column equals the checked snapshot's time.
    Never assume Delta_t == TimeIntegration:dt_max -- individual
    time-stepping can settle on a smaller step (see README)."""
    with open(log_path, "r") as f:
        for line in f:
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                t = float(fields[1])
                dt = float(fields[4])
            except ValueError:
                continue
            if abs(t - time) < 1e-4 * abs(time) or (time == 0.0 and t == 0.0):
                return dt
    raise RuntimeError(
        f"Could not find a step in {log_path!r} matching snapshot time {time:.6e}."
    )


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    # Use the last snapshot. LW_FUV_propagation is off, so the field is
    # reset and fully re-injected on every star-feedback pass (see README):
    # any snapshot after at least one pass is a clean, self-contained check
    # of that pass's own injection, independent of how many earlier passes
    # ran.
    snap_path = files[-1]
    with h5py.File(snap_path, "r") as f:
        time = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        gas = f["/PartType0"]
        mass = gas["Masses"][:].astype(np.float64)
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        Z = gas["MetalMassFractions"][:, -1]

        star = f["/PartType4"]
        L_FUV = float(star["FUVLuminosities"][0])
        L_LW = float(star["LWLuminosities"][0])

    if np.any(Z != 0.0):
        raise RuntimeError(
            "This check requires GEARChemistry:initial_metallicity=0 (exact "
            "extinction=1); found nonzero metallicity in the snapshot."
        )

    Delta_t = read_delta_t(opt.log, time)

    n_illuminated = int(np.sum(u_fuv > 0))
    sum_fuv = float(np.sum(u_fuv * mass))
    sum_lw = float(np.sum(u_lw * mass))
    rhs_fuv = Delta_t * L_FUV
    rhs_lw = Delta_t * L_LW

    print(f"Snapshot: {snap_path} (t={time:.6e})")
    print(f"Delta_t (from {opt.log}): {Delta_t:.6e}")
    print(f"Star: L_FUV={L_FUV:.10e}, L_LW={L_LW:.10e}")
    print(f"Illuminated gas particles: {n_illuminated}")

    def report(band, lhs, rhs):
        rel_err = abs(lhs - rhs) / rhs
        status = "PASS" if rel_err < opt.tol else "FAIL"
        print(
            f"{band}: sum(u*m)={lhs:.10e}, Delta_t*L={rhs:.10e}, "
            f"rel_err={rel_err:.3e} -> {status}"
        )
        return rel_err < opt.tol

    ok_fuv = report("FUV", sum_fuv, rhs_fuv)
    ok_lw = report("LW", sum_lw, rhs_lw)

    if not (ok_fuv and ok_lw):
        sys.exit(1)


if __name__ == "__main__":
    main()
