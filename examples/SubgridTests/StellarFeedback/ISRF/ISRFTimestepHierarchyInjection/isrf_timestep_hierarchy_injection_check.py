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

Gated, the time-bin structure of each run, read from timesteps.txt: the
hierarchy run updates part of the gas on some steps and all of it on
others, with the coarse slab on at least twice the fine slab's step and
the star step longer than the shortest gas step; the single-bin run
updates every gas particle on every step; with a pinned c_hyp,
c_hyp dt / h stays at or below ISRF_c_hyp_margin in each slab; the gas
carries no metals, so nothing absorbs.

Gated, the injection budget. With Z = 0 and a pinned c_hyp the stored
energy E = (c / c_hyp) sum m u trails the emitted energy L t only by the
dose still in flight between the star and the gas. That lag is bounded,
and the bound follows from the cadences the run itself ran at, see
lag_window. The lag is gated against that window at every snapshot after
the first. It is NOT fitted: a least-squares slope of E against L t reads
the lag's step quantisation as a drift and reports a spurious violation
even on a single-bin control, so the slope is printed as a diagnostic and
compared against nothing.

Every gated value is tested for finiteness before it is compared, and the
script exits nonzero on any failure.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

C_LIGHT_CGS = 2.99792458e10

# Relative rounding of the float32 specific energies the snapshots carry.
FLT_EPSILON = float(np.finfo(np.float32).eps)


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
    """Print and return one upper-bound gate; non-finite values fail."""
    ok = bool(np.isfinite(value) and np.isfinite(bar) and value <= bar)
    print(f"  {name}: {value:.4e} <= {bar:.4e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate_window(name: str, value: float, low: float, high: float) -> bool:
    """Print and return one two-sided gate; non-finite values fail.

    Parameters
    ----------
    name : str
        Label to print.
    value : float
        Measured value.
    low, high : float
        Window edges, inclusive.

    Returns
    -------
    bool
        True if the value is finite and inside a finite window.
    """
    finite = bool(np.isfinite(value) and np.isfinite(low) and np.isfinite(high))
    ok = bool(finite and low <= value <= high)
    print(
        f"  {name}: {value:.4e} in [{low:.4e}, {high:.4e}] -> "
        f"{'PASS' if ok else 'FAIL'}"
    )
    return ok


def steps_of(times: np.ndarray, longest: bool = False) -> float:
    """Return the shortest, or longest, spacing of a set of step end times.

    Parameters
    ----------
    times : numpy.ndarray
        End times of the steps on which some set of particles was updated.
    longest : bool
        Return the longest spacing instead of the shortest.

    Returns
    -------
    float
        Spacing, or NaN if the set is empty.
    """
    if times.size == 0:
        return float("nan")
    gaps = np.diff(np.concatenate([[0.0], times]))
    return float(np.max(gaps) if longest else np.min(gaps))


def lag_window(dt_gas: float, dt_star: float) -> tuple:
    """Return the window the in-flight injection budget must stay inside.

    The emitted energy L t reaches the gas through three stages, each of
    which holds at most one cadence's worth of dose, so the lag L t - E is
    bounded by their sum.

    1. Emission awaiting hand-off. The star hands its dose out once per
       star step, in one lump covering its whole step
       (`Delta_t = get_timestep(si->time_bin)`, radiation_iact.h). Between
       two hand-offs the emission clock advances while no dose moves, so up
       to `L dt_star` of emission sits in no reservoir.
    2. Dose awaiting drawdown. A lump lands in each illuminated particle's
       `u_dose_reservoir` and is drawn down only at that particle's next
       `feedback_reset_part`, for active particles only
       (radiation_snapshot_part_propagation, radiation_isrf.c). That wait
       is at most one gas step, so up to `L dt_gas`.
    3. Drawn dose awaiting the u update. The drawdown turns the reservoir
       into a rate released over a window ending at
       `ISRF_reservoir_end_ti`, the end of the star step the dose belongs
       to, and `radiation_end_force_propagation` adds one step of it per
       particle step. The residue is at most one release window, so up to
       `L max(dt_gas, dt_star)`.

    All three hold energy the gas has not yet stored, so they only add. A
    lump is handed out for a star step that has not elapsed yet, however,
    so a particle that releases it in one go can store energy the emission
    clock has not reached: a lead of at most one star step.

    The window is therefore `[-dt_star, dt_star + dt_gas + max(dt_gas,
    dt_star)]`, in units of the luminosity, that is, a time.

    Parameters
    ----------
    dt_gas : float
        Longest interval between two steps on which every gas particle was
        updated, an upper bound on any particle's own step.
    dt_star : float
        Longest interval between two star updates.

    Returns
    -------
    tuple
        (lower edge, upper edge), both times.
    """
    return -dt_star, dt_star + dt_gas + max(dt_gas, dt_star)


def preconditions(run: str, snap: dict, hierarchy: bool) -> tuple:
    """Check the time-bin structure of a run and measure its cadences.

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
        (ok, fine step, coarse step, star step, longest all-gas interval,
        longest star interval).
    """
    table = step_table(run)[1:]
    n_gas = len(snap["mass"])
    n_fine = int(snap["fine"].sum())
    all_gas = table[:, 2] == n_gas
    most_fine = table[:, 2] >= 0.8 * n_fine
    star_on = table[:, 3] > 0
    dt_min = float(np.min(table[:, 1]))
    dt_fine = steps_of(table[most_fine, 0])
    dt_coarse = steps_of(table[all_gas, 0])
    dt_star = steps_of(table[star_on, 0])
    dt_gas_longest = steps_of(table[all_gas, 0], longest=True)
    dt_star_longest = steps_of(table[star_on, 0], longest=True)
    finite = bool(
        np.all(
            np.isfinite(
                [dt_min, dt_fine, dt_coarse, dt_star, dt_gas_longest, dt_star_longest]
            )
        )
    )
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
        f"{dt_star:.4e}, longest all-gas interval {dt_gas_longest:.4e}, "
        f"longest star interval {dt_star_longest:.4e} -> "
        f"{'PASS' if ok else 'FAIL'}"
    )
    return ok, dt_fine, dt_coarse, dt_star, dt_gas_longest, dt_star_longest


def main() -> None:
    """Run every gate and exit nonzero on any failure."""
    opt = parse_options()
    run = lambda name: os.path.join(opt.dir, name)
    ok = True

    print("Time-bin structure and absorption-free setup")
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

    print("Injection budget in flight, L t - E against its derived window")
    for name in ("conservation_single_bin", "conservation_hierarchy"):
        snaps = [load(f) for f in snapshots(run(name))[1:]]
        pin = parameter(run(name), "ISRF_c_hyp_pin_for_debugging")
        low, high = lag_window(info[name][4], info[name][5])
        n_gas = len(snaps[0]["mass"])
        t = np.array([s["time"] for s in snaps])
        print(
            f"  {name}: window [{low:.4e}, {high:.4e}] from a longest all-gas "
            f"interval of {info[name][4]:.4e} and a longest star interval of "
            f"{info[name][5]:.4e}"
        )
        for band in ("PE", "LW"):
            L = np.array([s["L"][band] for s in snaps])
            # L t is the emitted energy only for a constant L. A float32
            # luminosity that did not move is bit-identical, so any spread
            # above one rounding unit means it did.
            positive = bool(np.all(np.isfinite(L)) and np.all(L > 0.0))
            print(
                f"  {name} {band} luminosity finite and positive: "
                f"{'PASS' if positive else 'FAIL'}"
            )
            ok &= positive
            spread = float(np.max(np.abs(L - L[0])) / max(abs(L[0]), 1e-300))
            ok &= gate(f"{name} {band} luminosity spread", spread, FLT_EPSILON)
            E = np.array(
                [np.sum(s["mass"] * s["u"][band]) * s["c"] / pin for s in snaps]
            )
            lag = t - E / L
            # E is a sum of N float32 specific energies, so its rounding
            # error grows as sqrt(N) eps E; carried to the lag through L,
            # that is sqrt(N) eps t.
            tol = np.sqrt(n_gas) * FLT_EPSILON * t
            # These reductions must propagate a non-finite lag to the gate,
            # so keep np.min/np.max: a nan-ignoring variant fails open.
            ok &= gate_window(
                f"{name} {band} min lag", float(np.min(lag + tol)), low, high
            )
            ok &= gate_window(
                f"{name} {band} max lag", float(np.max(lag - tol)), low, high
            )
            slope, intercept = np.polyfit(t, E / L, 1)
            resid = E / L - (slope * t + intercept)
            se = np.sqrt(np.sum(resid**2) / (len(t) - 2) / np.sum((t - t.mean()) ** 2))
            print(
                f"  {name} {band}: {len(t)} snapshots, lag / window "
                f"{np.min(lag) / high:.3f} to {np.max(lag) / high:.3f}, lag / "
                f"longest all-gas interval {np.min(lag) / info[name][4]:.3f} to "
                f"{np.max(lag) / info[name][4]:.3f}"
            )
            print(
                f"  {name} {band}: diagnostic fit, s - 1 = {slope - 1.0:.3e}, "
                f"se = {se:.3e} (quantised residual, not a gate)"
            )

    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
