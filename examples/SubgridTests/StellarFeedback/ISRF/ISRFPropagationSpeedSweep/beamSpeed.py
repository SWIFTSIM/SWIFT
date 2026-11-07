#!/usr/bin/env python3
"""Measure the speed of the pulse of light in one or more runs, as a fraction of the set light speed.

The pulse position is the energy centroid of the light in flight (inbox, reservoir and held energy of the
PE band) along the beam direction. The speed is the advance of the centroid between two snapshots over their
time difference, divided by the set speed (the fraction of c of the run). With two runs the script also prints
the ratio of the second to the first.
"""

import argparse
import glob
import json
import sys

import h5py
import numpy as np

C_KMS = 299792.458
FIELDS = ("PacketInboxEnergies", "PacketReservoirEnergies", "PacketHeldEnergies")


def centroid(path, geometry):
    """Return the time and the centroid of the light along the beam, in kpc from the box centre."""
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        try:
            light = sum(np.asarray(g[k][:, 0], dtype=np.float64) for k in FIELDS)
        except KeyError as err:
            sys.exit(
                f"{path}: field {err} missing: configure with --enable-isrf-debugging-checks."
            )
        x = np.asarray(g["Coordinates"][:], dtype=np.float64)
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
        box = float(np.ravel(f["Header"].attrs["BoxSize"])[0])
    if not (np.all(np.isfinite(light)) and light.sum() > 0.0):
        sys.exit(f"{path}: no finite light in flight.")
    rel = x - np.asarray(geometry["centre_kpc"])
    rel -= box * np.round(rel / box)
    return t, float(
        (light * (rel @ np.asarray(geometry["direction"]))).sum() / light.sum()
    )


def speed(run, fraction, first, last):
    """Return the pulse speed of a run over the set speed, between two snapshot numbers."""
    geometry = json.load(open(f"{run}/geometry.json"))
    files = sorted(glob.glob(f"{run}/snap/snapshot_*.hdf5"))
    if len(files) <= last:
        sys.exit(f"{run}: snapshot {last} missing ({len(files)} found).")
    t0, c0 = centroid(files[first], geometry)
    t1, c1 = centroid(files[last], geometry)
    return (c1 - c0) / (t1 - t0) / (fraction * C_KMS)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("runs", nargs="+", help="Run directories")
    p.add_argument(
        "--fraction",
        type=float,
        default=3e-5,
        help="Set light speed as a fraction of c",
    )
    p.add_argument(
        "--first", type=int, default=2, help="First snapshot number of the window"
    )
    p.add_argument(
        "--last", type=int, default=5, help="Last snapshot number of the window"
    )
    p.add_argument(
        "--bar",
        type=float,
        default=None,
        help="With two runs, exit with status 1 if the ratio differs from 1 by more than this",
    )
    opt = p.parse_args()
    speeds = [speed(r, opt.fraction, opt.first, opt.last) for r in opt.runs]
    for r, s in zip(opt.runs, speeds):
        print(f"{r}: speed over set speed {s:.4f}")
    if len(speeds) == 2:
        print(f"ratio {opt.runs[1]} over {opt.runs[0]}: {speeds[1] / speeds[0]:.3f}")
        if opt.bar is not None and not abs(speeds[1] / speeds[0] - 1.0) <= opt.bar:
            sys.exit(1)


if __name__ == "__main__":
    main()
