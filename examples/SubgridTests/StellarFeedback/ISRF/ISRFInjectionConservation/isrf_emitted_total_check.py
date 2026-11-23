#!/usr/bin/env python3
"""Check that a star emitted the light it owes, with the packet transport on.

A star of constant luminosity L, born at time t_b, emits L * (t - t_b) of light in each band by
time t. The code books this light step by step into PacketEmittedEnergies. A lost or doubled
step, or a start-up that books the first step as zero, changes the booked total and leaves the
light budget of ISRFCosmology/compareLegs.py closed, because that budget only compares the
light booked with the light accounted for. This check compares the booked total with the clock:

    R = PacketEmittedEnergies / (L * (t - t_b)),   target 1, per star and per band (PE, LW)

L is the luminosity of the star in the snapshot and t is the time of the last snapshot. Every
star step ends at that time, because the run stops at TimeIntegration:time_end and the star
steps are bounded by it. A snapshot inside a star step reads low by the light of the part of
the step not yet booked: use the last snapshot only.

The test fails when |R - 1| exceeds the bar in any star and band, or on a non-finite value. It
refuses a run it cannot judge: no transport, a cosmological run, a last snapshot not at time_end,
a luminosity that changes between the last two snapshots (the formula needs a constant L, as
for a single star), a star count that changes, or a star with no age.

The bar 1e-6 is that of the light budget. L is stored in float32 (relative error 6e-8, one
rounding). The booked total gathers one float32 sum per star update, 2^-24 = 6e-8 each. The most
updates in a cadence of the example is 32, which gives 3.4e-7 for independent roundings (sqrt(32)
x 6e-8). The error measured at the three cadences is 2.5e-8. The smallest real defect, the first
step lost at the finest step limit of 2.5e-5 over 4e-4, is 6.25e-2. The check excludes the way
the light is spread in time, and a wrong L, which cancels in R.
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np
import yaml


def read_snapshot(path: str) -> dict:
    """Return the time, the star ids, the birth times, L and the booked light of a snapshot.

    Parameters
    ----------
    path : str
        Snapshot file.

    Returns
    -------
    dict
        Keys ``t``, ``ids``, ``birth``, ``L`` (stars by 2: PE, LW) and ``emitted`` (same shape),
        the stars sorted by id, and ``cosmological``.
    """
    with h5py.File(path, "r") as f:
        star = f["PartType4"]
        try:
            emitted = np.asarray(star["PacketEmittedEnergies"][:], dtype=np.float64)
        except KeyError:
            sys.exit(
                f"{path}: PacketEmittedEnergies missing: configure with --enable-isrf-debugging-checks."
            )
        ids = np.asarray(star["ParticleIDs"][:])
        order = np.argsort(ids, kind="stable")
        luminosity = np.stack(
            [
                np.asarray(star["PELuminosities"][:], dtype=np.float64),
                np.asarray(star["LWLuminosities"][:], dtype=np.float64),
            ],
            axis=1,
        )
        return {
            "t": float(np.ravel(f["Header"].attrs["Time"])[0]),
            "ids": ids[order],
            "birth": np.asarray(star["BirthTimes"][:], dtype=np.float64)[order],
            "L": luminosity[order],
            "emitted": emitted[order],
            "cosmological": int(np.ravel(f["Cosmology"].attrs["Cosmological run"])[0])
            != 0,
        }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("run", nargs="?", default=".", help="Run directory (default: .)")
    p.add_argument("--bar", type=float, default=1e-6, help="Largest |R - 1| allowed")
    opt = p.parse_args()
    if not (np.isfinite(opt.bar) and opt.bar > 0.0):
        sys.exit("The bar must be finite and positive.")

    with open(os.path.join(opt.run, "used_parameters.yml")) as handle:
        used = yaml.safe_load(handle)
    if int(used["GEARFeedback"]["ISRF_propagation"]) != 1:
        sys.exit(
            "Refused: GEARFeedback:ISRF_propagation is 0, so no light is booked into PacketEmittedEnergies."
        )
    time_end = float(used["TimeIntegration"]["time_end"])

    files = sorted(glob.glob(os.path.join(opt.run, "snap", "snapshot_*.hdf5")))
    if len(files) < 2:
        sys.exit("Need at least two snapshots.")
    before, last = read_snapshot(files[-2]), read_snapshot(files[-1])

    for name, snap in (("last", last), ("second to last", before)):
        for key in ("t", "birth", "L", "emitted"):
            if not np.all(np.isfinite(snap[key])):
                sys.exit(f"FAIL: non-finite {key} in the {name} snapshot.")
    if last["cosmological"]:
        sys.exit("Refused: the run is cosmological, t is not the time.")
    if abs(last["t"] - time_end) > 1e-9 * time_end:
        sys.exit(
            f"Refused: the last snapshot is at t={last['t']:.6e}, not at time_end={time_end:.6e}."
        )
    if last["ids"].size == 0 or not np.array_equal(last["ids"], before["ids"]):
        sys.exit(
            "Refused: no star, or the star set differs between the last two snapshots."
        )
    if np.any(np.abs(last["L"] - before["L"]) > 1e-6 * np.abs(last["L"])):
        sys.exit("Refused: the luminosity changes between the last two snapshots.")
    age = last["t"] - last["birth"]
    if not np.all(age > 0.0) or not np.all(last["L"] > 0.0):
        sys.exit("Refused: a star has no age or no luminosity.")

    ratio = last["emitted"] / (last["L"] * age[:, None])
    deviation = np.abs(ratio - 1.0)
    print(f"t={last['t']:.6e} stars={last['ids'].size} bar={opt.bar:.0e}")
    for i, star_id in enumerate(last["ids"]):
        print(
            f"star {star_id}: emitted/(L age) PE={ratio[i, 0]:.9f} LW={ratio[i, 1]:.9f}"
        )
    worst = float(deviation.max())
    verdict = "pass" if worst <= opt.bar else "FAIL"
    print(f"worst |R-1|={worst:.3e} {verdict}")
    return 0 if verdict == "pass" else 1


if __name__ == "__main__":
    sys.exit(main())
