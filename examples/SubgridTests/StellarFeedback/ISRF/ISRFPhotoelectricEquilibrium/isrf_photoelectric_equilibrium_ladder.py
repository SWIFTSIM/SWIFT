#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
"""Equilibrium curves P(n) and T(n) from a ladder of runs, and the two-phase shape check.

Run directories are named n<density>_hot and n<density>_cold (density in H atoms/cm^3),
each holding snap/. Each direction (hot start, cold start) gives one branch. The shape
gates follow from the existence of a thermally unstable density range, with no number
taken from a paper:

  hysteresis  at some density in 0.16 to 0.20 H/cm^3 (the window of the one-zone
              prediction, fixed before the runs) the hot-start and cold-start specific
              energies differ by a factor of at least 1.5
  maximum     along the hot-start branch the pressure has a first local maximum, and the
              cold-start branch a local minimum after it, as in a two-phase medium

The hysteresis gate asks for a factor of at least 1.5 between the two branches; the maximum gate asks for a
maximum of P(n) followed by a minimum.

Only runs whose check.txt holds no FAIL line and no SPLIT mark count; the others are listed
as excluded. A SPLIT run (p95/p5 of the specific energy above 1.5) mixes two phases, so its
median is no equilibrium, and the hot/cold pair of that density leaves both gates.
With --onezone-table the same gates run on the one-zone ladder (lines FINAL of onezone).
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

BOLTZMANN_CGS = 1.380649e-16
HYSTERESIS_RATIO = 1.5
HYSTERESIS_WINDOW = (0.16, 0.20)
SPLIT_MARK = "SPLIT into"


def run_verdict(check):
    """Return "ok", "no check", "a failed gate" or "SPLIT" for the check.txt of a run."""
    if not os.path.exists(check):
        return "no check"
    with open(check) as f:
        lines = f.read().splitlines()
    if any(line.startswith("FAIL") for line in lines):
        return "a failed gate"
    if any(SPLIT_MARK in line for line in lines):
        return "SPLIT"
    return "ok"


def from_runs(root):
    """Return {(n, start): (specific energy [erg/g], P/k)} from the last snapshots of a ladder."""
    table = {}
    for path in sorted(glob.glob(os.path.join(root, "n*_*"))):
        m = re.match(r"n([0-9.eE+-]+)_(hot|cold)$", os.path.basename(path))
        files = sorted(glob.glob(os.path.join(path, "snap", "*.hdf5")))
        if not m or not files:
            continue
        check = os.path.join(path, "check.txt")
        verdict = run_verdict(check)
        if verdict != "ok":
            print(f"excluded ({verdict}): {os.path.basename(path)}")
            continue
        with h5py.File(files[-1], "r") as f:
            units = f["Units"].attrs
            u_mass = float(np.ravel(units["Unit mass in cgs (U_M)"])[0])
            u_len = float(np.ravel(units["Unit length in cgs (U_L)"])[0])
            u_time = float(np.ravel(units["Unit time in cgs (U_t)"])[0])
            g = f["PartType0"]
            pressure = (
                np.median(g["Pressures"][:].astype(np.float64))
                * u_mass
                / u_len
                / u_time**2
            )
            energy = (
                np.median(g["InternalEnergies"][:].astype(np.float64))
                * (u_len / u_time) ** 2
            )
        table[(float(m.group(1)), m.group(2))] = (energy, pressure / BOLTZMANN_CGS)
    return table


def from_onezone(path):
    """Return the same table from the FINAL lines of the one-zone ladder (T0 below 1e3 K is cold)."""
    table = {}
    for line in open(path):
        if not line.startswith("FINAL"):
            continue
        d = {k: float(v) for k, v in (item.split("=") for item in line.split()[1:])}
        table[(d["n_H"], "cold" if d["T0"] < 1e3 else "hot")] = (d["T"], d["P_over_k"])
    return table


def first_extremum(values, sign, start=1):
    """Index of the first interior local maximum (sign 1) or minimum (sign -1), or -1."""
    for i in range(max(start, 1), len(values) - 1):
        if (
            sign * values[i] > sign * values[i - 1]
            and sign * values[i] >= sign * values[i + 1]
        ):
            return i
    return -1


def shape_gates(table):
    """Return (name, value, passed) for the hysteresis and the pressure extrema."""
    dens = sorted({n for n, _ in table})
    both = [n for n in dens if (n, "hot") in table and (n, "cold") in table]
    in_window = [
        n
        for n in both
        if HYSTERESIS_WINDOW[0] - 1e-9 <= n <= HYSTERESIS_WINDOW[1] + 1e-9
    ]
    ratio = [table[(n, "hot")][0] / table[(n, "cold")][0] for n in in_window]
    max_ratio = max(ratio) if ratio else float("nan")
    hot_p = np.array([table[(n, "hot")][1] for n in both])
    cold_p = np.array([table[(n, "cold")][1] for n in both])
    i_hot = first_extremum(hot_p, 1)
    i_cold = first_extremum(cold_p, -1, start=i_hot + 1) if i_hot >= 0 else -1
    interior_max = i_hot >= 0
    interior_min = i_cold >= 0
    return [
        (
            "hysteresis in 0.16 to 0.20: max hot/cold ratio",
            max_ratio,
            bool(np.isfinite(max_ratio) and max_ratio >= HYSTERESIS_RATIO),
        ),
        (
            "first hot-start pressure maximum is interior (P/k at n = %.3g)"
            % (both[i_hot] if i_hot >= 0 else float("nan")),
            float(hot_p[i_hot]) if i_hot >= 0 else float("nan"),
            bool(interior_max),
        ),
        (
            "first cold-start pressure minimum after it is interior (P/k at n = %.3g)"
            % (both[i_cold] if i_cold >= 0 else float("nan")),
            float(cold_p[i_cold]) if i_cold >= 0 else float("nan"),
            bool(interior_min),
        ),
    ]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "root", nargs="?", help="Ladder directory of n<density>_<hot|cold> runs"
    )
    parser.add_argument(
        "--onezone-table", help="Use the FINAL lines of a one-zone ladder instead"
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run the gates on a flat ladder, which must fail",
    )
    opt = parser.parse_args()
    if opt.onezone_table:
        table = from_onezone(opt.onezone_table)
    elif opt.root:
        table = from_runs(opt.root)
    else:
        parser.error("give a ladder directory or --onezone-table")
    if not table or not all(np.isfinite(v) for pair in table.values() for v in pair):
        print("FAIL empty or non-finite ladder")
        sys.exit(1)
    for (n, start), (t, p) in sorted(table.items()):
        print(f"n = {n:<8g} {start:<5} T-or-u = {t:<12.6g} P/k = {p:.6g}")
    gates = shape_gates(table)
    if opt.self_test:
        flat = {(n, s): (100.0 + 10.0 * np.log10(n), 100.0 * n) for (n, s) in table}
        gates += [
            (f"self-test must fail: {g[0]}", g[1], not g[2]) for g in shape_gates(flat)
        ]
    failed = 0
    for name, value, passed in gates:
        print(f"{'PASS' if passed else 'FAIL'}  {name}: {value:.4g}")
        failed += not passed
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
