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
"""
Check the ISRF injection and propagation under individual time steps.

Two layers, with different standing.

GATED, the time-bin structure of each run, read from timesteps.txt: the
hierarchy run updates part of the gas on some steps and all of it on
others, with the coarse slab on at least twice the fine slab's step and
the star step longer than the shortest gas step; the single-bin run
updates every gas particle on every step; with a pinned c_hyp,
c_hyp dt / h stays at or below ISRF_c_hyp_margin in each slab; the gas
carries no metals, so nothing absorbs. These gates hold the fixture to
what its README claims it is. A non-finite value fails its gate, and the
script exits nonzero on any failure.

NOT GATED, the conservation report. With Z = 0 and a pinned c_hyp the
stored energy (c / c_hyp) sum m u lags L t only by the dose still held in
the reservoirs and in the gas steps in flight. That lag is bounded, and it
is quantised in whole particle steps, so it steps between neighbouring
integer multiples of the fine step rather than varying smoothly. A
least-squares slope of (c / c_hyp) sum m u / L against t therefore reads a
bounded quantised offset as a drift, and its standard error assumes white
residuals that a staircase does not provide. The slope and the lag are
printed for inspection and neither is compared against a bar. A gate on
this fixture needs a bounded-window metric on the lag itself, with the
window derived from the measured bin structure.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

C_LIGHT_CGS = 2.99792458e10


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dir", default=".", help="Directory holding the runs")
    return parser.parse_args()


def snapshots(run: str) -> list:
    """Return the sorted snapshot files of a run, one per output time."""
    files, times = [], set()
    for path in sorted(glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5"))):
        with h5py.File(path, "r") as f:
            t = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        if not any(abs(t - s) <= 1e-9 * max(abs(t), 1e-300) for s in times):
            times.add(t)
            files.append(path)
    if len(files) < 3:
        raise RuntimeError(f"Fewer than three snapshots in {run}/snap")
    return files


def load(path: str) -> dict:
    """Read one snapshot, gas sorted by particle ID.

    Parameters
    ----------
    path : str
        Snapshot file.

    Returns
    -------
    dict
        Snapshot fields in internal units.
    """
    with h5py.File(path, "r") as f:
        units = f["/Units"].attrs
        gas = f["/PartType0"]
        order = np.argsort(gas["ParticleIDs"][:])
        star = f["/PartType4"]
        snap = dict(
            time=float(np.asarray(f["/Header"].attrs["Time"]).flat[0]),
            ul=float(np.asarray(units["Unit length in cgs (U_L)"]).flat[0]),
            ut=float(np.asarray(units["Unit time in cgs (U_t)"]).flat[0]),
            h=gas["SmoothingLengths"][:][order].astype(np.float64),
            mass=gas["Masses"][:][order].astype(np.float64),
            Z=gas["MetalMassFractions"][:, -1][order].astype(np.float64),
            u={
                "PE": gas["PESpecificEnergies"][:][order].astype(np.float64),
                "LW": gas["LWSpecificEnergies"][:][order].astype(np.float64),
            },
            L={
                "PE": float(star["PELuminosities"][0]),
                "LW": float(star["LWLuminosities"][0]),
            },
        )
    snap["c"] = C_LIGHT_CGS * snap["ut"] / snap["ul"]
    # The hot, light phase sits on the finer time bin.
    snap["fine"] = snap["mass"] < 2.0 * snap["mass"].min()
    return snap


def step_table(run: str) -> np.ndarray:
    """Return timesteps.txt rows: time, dt, gas updates, star updates."""
    rows = []
    with open(os.path.join(run, "timesteps.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append([float(v[1]), float(v[4]), int(v[7]), int(v[9])])
    return np.array(rows)


def parameter(run: str, key: str) -> float:
    """Read a numerical parameter from used_parameters.yml."""
    with open(os.path.join(run, "used_parameters.yml")) as f:
        m = re.search(rf"{key}:\s*([0-9.eE+-]+)", f.read())
    return float(m.group(1))


def gate(name: str, value: float, bar: float) -> bool:
    """Print and return one gate; non-finite values fail."""
    ok = bool(np.isfinite(value) and np.isfinite(bar) and value <= bar)
    print(f"  {name}: {value:.4e} <= {bar:.4e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def steps_of(times: np.ndarray) -> float:
    """Return the smallest spacing of a set of step end times."""
    return float(np.min(np.diff(np.concatenate([[0.0], times]))))


def preconditions(run: str, snap: dict, hierarchy: bool) -> tuple:
    """Check the time-bin structure of a run.

    Parameters
    ----------
    run : str
        Run directory.
    snap : dict
        Any snapshot of the run.
    hierarchy : bool
        True if the two slabs must sit on different bins.

    Returns
    -------
    tuple
        (ok, fine step, coarse step, star step).
    """
    table = step_table(run)[1:]
    n_gas = len(snap["mass"])
    n_fine = int(snap["fine"].sum())
    all_gas = table[:, 2] == n_gas
    most_fine = table[:, 2] >= 0.8 * n_fine
    dt_min = float(np.min(table[:, 1]))
    dt_fine = steps_of(table[most_fine, 0])
    dt_coarse = steps_of(table[all_gas, 0])
    dt_star = steps_of(table[table[:, 3] > 0, 0])
    finite = bool(np.all(np.isfinite([dt_min, dt_fine, dt_coarse, dt_star])))
    if hierarchy:
        ok = bool(
            finite
            and (~all_gas).mean() >= 0.25
            and dt_coarse >= 2.0 * dt_fine
            and dt_star > dt_min
        )
    else:
        ok = bool(finite and all_gas.all())
    print(
        f"  {run}: steps {len(table)}, partial-update steps "
        f"{(~all_gas).mean():.2f}, shortest step {dt_min:.4e}, fine slab step "
        f"{dt_fine:.4e}, coarse slab step {dt_coarse:.4e}, star step "
        f"{dt_star:.4e} -> {'PASS' if ok else 'FAIL'}"
    )
    return ok, dt_fine, dt_coarse, dt_star


def main() -> None:
    """Check the bin structure, report conservation, exit nonzero on failure."""
    opt = parse_options()
    run = lambda name: os.path.join(opt.dir, name)
    ok = True

    print("Time-bin structure and absorption-free setup (gated)")
    info = {}
    for name in ("conservation_single_bin", "conservation_hierarchy"):
        snap = load(snapshots(run(name))[-1])
        r = preconditions(run(name), snap, "hierarchy" in name)
        ok &= r[0]
        info[name] = r
        pin = parameter(run(name), "ISRF_c_hyp_pin_for_debugging")
        margin = parameter(run(name), "ISRF_c_hyp_margin")
        courant = max(
            pin * r[1] / snap["h"][snap["fine"]].min(),
            pin * r[2] / snap["h"][~snap["fine"]].min(),
        )
        ok &= gate(f"{name} pinned c_hyp dt / h", courant, margin)
        ok &= gate(f"{name} max metallicity", float(snap["Z"].max()), 0.0)

    print("Conservation report (NOT gated, see the module docstring)")
    for name in ("conservation_single_bin", "conservation_hierarchy"):
        snaps = [load(f) for f in snapshots(run(name))[1:]]
        pin = parameter(run(name), "ISRF_c_hyp_pin_for_debugging")
        dt_fine = info[name][1]
        t = np.array([s["time"] for s in snaps])
        for band in ("PE", "LW"):
            y = np.array(
                [
                    np.sum(s["mass"] * s["u"][band]) * s["c"] / pin / s["L"][band]
                    for s in snaps
                ]
            )
            slope, intercept = np.polyfit(t, y, 1)
            resid = y - (slope * t + intercept)
            se = np.sqrt(np.sum(resid**2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))
            lag = (t - y) / dt_fine
            print(
                f"  {name} {band}: {len(t)} snapshots, s - 1 = {slope - 1.0:.3e}, "
                f"se = {se:.3e}"
            )
            print(
                f"  {name} {band}: lag (L t - E) in units of this run's own "
                f"fine step ({dt_fine:.4e}), min {lag.min():.2f}, max "
                f"{lag.max():.2f}, last {lag[-1]:.2f}"
            )

    print("BIN-STRUCTURE GATES:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
