#!/usr/bin/env python3
"""Compare the light field of a rest run and a boosted run, snapshot by snapshot.

D = max |u_boost - u_rest| / max u_rest, with u the PacketSpecificEnergies (both bands) matched by particle id.
The second column is the largest displacement of a particle between the two runs once the bulk motion is
removed, in median smoothing lengths. Compare two identical rest runs to get the repeat floor of D.
The default bar is 1.9e-6 = 32 x 2^-24: one float32 rounding step of D for each of the 32 time steps. It applies
to the lattice pairs. The glass pair gives D of about 7.4e-5 by design and fails it.
"""

import argparse
import glob
import sys

import h5py
import numpy as np


def load(path):
    """Return time, box, ids, positions, smoothing lengths and light field, sorted by particle id."""
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        o = np.argsort(g["ParticleIDs"][:], kind="stable")
        try:
            u = g["PacketSpecificEnergies"][:][o].astype(np.float64)
        except KeyError:
            sys.exit(
                "PacketSpecificEnergies missing: configure with --enable-isrf-debugging-checks."
            )
        return dict(
            t=float(np.ravel(f["Header"].attrs["Time"])[0]),
            box=float(np.ravel(f["Header"].attrs["BoxSize"])[0]),
            ids=g["ParticleIDs"][:][o],
            x=g["Coordinates"][:][o],
            h=g["SmoothingLengths"][:][o],
            u=u,
        )


def displacement(a, b):
    """Return the largest displacement between the runs after removing the mean shift along x, in median h."""
    box = a["box"]
    d = b["x"] - a["x"]
    d = (d + 0.5 * box) % box - 0.5 * box
    d[:, 0] -= np.angle(np.mean(np.exp(2j * np.pi * d[:, 0] / box))) / (2 * np.pi) * box
    d = (d + 0.5 * box) % box - 0.5 * box
    return float(np.abs(d).max() / np.median(a["h"]))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("rest", help="Run directory of the rest box")
    p.add_argument(
        "boost",
        help="Run directory of the boosted box (or of a repeat of the rest box)",
    )
    p.add_argument(
        "--bar",
        type=float,
        default=1.9e-6,
        help="Exit with status 1 if D of the last snapshot exceeds this (default 32 x 2^-24)",
    )
    opt = p.parse_args()

    fa = sorted(glob.glob(f"{opt.rest}/snap/snapshot_*.hdf5"))
    fb = sorted(glob.glob(f"{opt.boost}/snap/snapshot_*.hdf5"))
    if len(fa) < 2 or len(fa) != len(fb):
        sys.exit(
            f"Need the same number (at least 2) of snapshots: {len(fa)} and {len(fb)}."
        )
    d_last = np.nan
    print("snapshot  time        D           displacement [h]")
    for i, (pa, pb) in enumerate(zip(fa, fb)):
        if i == 0:
            continue
        a, b = load(pa), load(pb)
        if not (
            np.array_equal(a["ids"], b["ids"]) and np.isclose(a["t"], b["t"], rtol=1e-9)
        ):
            sys.exit(f"Snapshot {i}: the runs differ in particle set or time.")
        if not (np.all(np.isfinite(a["u"])) and np.all(np.isfinite(b["u"]))):
            sys.exit(f"Snapshot {i}: non-finite light field.")
        if not a["u"].max() > 0.0:
            sys.exit(f"Snapshot {i}: no light in the rest run.")
        d_last = float(np.abs(a["u"] - b["u"]).max() / a["u"].max())
        print(f"{i:8d}  {a['t']:.4e}  {d_last:.4e}  {displacement(a, b):.3e}")
    print(f"D at the last snapshot: {d_last:.3e}")
    if not d_last <= opt.bar:
        sys.exit(1)


if __name__ == "__main__":
    main()
