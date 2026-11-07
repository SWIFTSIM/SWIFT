#!/usr/bin/env python3
"""Check that the transport of the light is linear in its energy, snapshot by snapshot.

Run the pulse at energy x1 and at energy x SCALE. The light in flight of the scaled run divided by SCALE must equal
the light in flight of the x1 run, up to float32 rounding:

    L = sum_i |E_scaled,i / SCALE - E_x1,i| / sum_i |E_x1,i|

E_i is the light in flight of particle i in the PE band (inbox + reservoir + held energy). The test fails when L
exceeds the bar in any snapshot, or on a non-finite energy, a different snapshot set, time or particle set, or no
light. The bar 1e-5 is a few float32 roundings of 6e-8 per hop over about 30 hops. The light-free first snapshot is
not compared.
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

FIELDS = ("PacketInboxEnergies", "PacketReservoirEnergies", "PacketHeldEnergies")


def read(path):
    """Return the particle ids, the light in flight per particle (sorted by id) and the time."""
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        ids = np.asarray(g["ParticleIDs"][:])
        o = np.argsort(ids, kind="stable")
        try:
            e = sum(np.asarray(g[k][:, 0], dtype=np.float64) for k in FIELDS)[o]
        except KeyError as err:
            sys.exit(
                f"{path}: field {err} missing: configure with --enable-isrf-debugging-checks."
            )
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
    return ids[o], e, t


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("x1", help="Run directory at energy x1")
    p.add_argument("scaled", help="Run directory at energy x SCALE")
    p.add_argument("scale", type=float, help="Energy factor of the scaled run")
    p.add_argument("--bar", type=float, default=1e-5, help="Largest L allowed")
    opt = p.parse_args()
    if not (
        np.isfinite(opt.scale)
        and opt.scale > 0.0
        and np.isfinite(opt.bar)
        and opt.bar > 0.0
    ):
        sys.exit("The scale and the bar must be finite and positive.")

    names = [
        os.path.basename(f) for f in sorted(glob.glob(f"{opt.x1}/snap/snapshot_*.hdf5"))
    ]
    if names != [
        os.path.basename(f)
        for f in sorted(glob.glob(f"{opt.scaled}/snap/snapshot_*.hdf5"))
    ]:
        sys.exit("The runs hold different snapshot sets.")
    if len(names) < 2:
        sys.exit("Need at least two snapshots.")

    worst = 0.0
    for name in names[1:]:
        ix, ex, tx = read(f"{opt.x1}/snap/{name}")
        iw, ew, tw = read(f"{opt.scaled}/snap/{name}")
        if not (np.all(np.isfinite(ex)) and np.all(np.isfinite(ew))):
            sys.exit(f"{name}: non-finite energy.")
        if tx != tw or not np.array_equal(ix, iw):
            sys.exit(f"{name}: the runs differ in time or particle set.")
        light = float(np.abs(ex).sum())
        if not light > 0.0:
            sys.exit(f"{name}: no light in the x1 run.")
        loss = float(np.abs(ew / opt.scale - ex).sum()) / light
        worst = max(worst, loss)
        print(f"{name} t={tx:.6e} L={loss:.3e} light={light:.4e}")
    verdict = "pass" if worst <= opt.bar else "FAIL"
    print(f"worst L={worst:.3e} bar={opt.bar:.0e} {verdict}")
    if verdict == "FAIL":
        sys.exit(1)


if __name__ == "__main__":
    main()
