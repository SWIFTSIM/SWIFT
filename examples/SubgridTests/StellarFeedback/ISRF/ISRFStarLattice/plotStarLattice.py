#!/usr/bin/env python3
"""Plot the light of a 4 x 4 x 4 lattice of stars in the plane of one layer of stars."""

import glob
import json
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def snapshots(pattern="snap/snapshot_*.hdf5"):
    """Return the sorted snapshot file names."""
    files = sorted(glob.glob(pattern))
    if len(files) < 2:
        sys.exit(f"Need at least two snapshots matching {pattern}.")
    return files


def load(path):
    """Return time, gas positions, gas ids and the PE-band specific energy of the gas."""
    with h5py.File(path, "r") as f:
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
        x = f["PartType0/Coordinates"][:]
        ids = f["PartType0/ParticleIDs"][:]
        q = np.ravel(f["PartType0/PESpecificEnergies"][:])
    return t, x, ids, q


def main():
    t, x, _, q = load(snapshots()[-1])
    with h5py.File(snapshots()[-1], "r") as f:
        box = float(f["Header"].attrs["BoxSize"][0])
        stars = f["PartType4/Coordinates"][:]
        h0 = float(f["PartType0/SmoothingLengths"][:].mean())
    z_layer = 3.0 * box / 8.0
    ok = (np.abs(x[:, 2] - z_layer) < 1.5 * h0) & (q > 0)
    fig, ax = plt.subplots(figsize=(5.2, 4.4), constrained_layout=True)
    im = ax.scatter(
        x[ok, 0] / box,
        x[ok, 1] / box,
        c=np.log10(q[ok]),
        s=3,
        marker="s",
        cmap="inferno",
        rasterized=True,
    )
    layer = np.abs(stars[:, 2] - z_layer) < 0.5 * h0
    ax.plot(stars[layer, 0] / box, stars[layer, 1] / box, "w*", ms=8, mec="k")
    fig.colorbar(im, label="log10 PE specific energy (internal units)")
    ax.set(aspect="equal", xlabel="x / box", ylabel="y / box", title=f"t = {t:.3g}")
    fig.savefig("star_lattice.png", dpi=150)
    print("star_lattice.png written")


if __name__ == "__main__":
    main()
