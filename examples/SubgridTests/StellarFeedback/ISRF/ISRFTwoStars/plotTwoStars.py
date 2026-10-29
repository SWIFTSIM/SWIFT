#!/usr/bin/env python3
"""Plot the light of two stars: a map of the plane of the stars and a cut along the line through them."""

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


def load(path, field="PacketSpecificEnergies"):
    """Return time, gas positions, gas ids and the PE-band value of a packet field."""
    with h5py.File(path, "r") as f:
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
        x = f["PartType0/Coordinates"][:]
        ids = f["PartType0/ParticleIDs"][:]
        q = f[f"PartType0/{field}"][:][:, 0]
    return t, x, ids, q


def main():
    t, x, _, q = load(snapshots()[-1], "PacketKernelMeanSpecificEnergies")
    with h5py.File(snapshots()[-1], "r") as f:
        box = float(f["Header"].attrs["BoxSize"][0])
        h0 = float(f["PartType0/SmoothingLengths"][:].mean())
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0), constrained_layout=True)
    c = x - 0.5 * box
    slab = np.abs(c[:, 2]) < 1.5 * h0
    ok = slab & (q > 0)
    im = axes[0].scatter(c[ok, 0] / box, c[ok, 1] / box, c=np.log10(q[ok]), s=1, marker="s", cmap="inferno", rasterized=True)
    fig.colorbar(im, ax=axes[0], label="log10 PE specific energy (internal units)")
    axes[0].set(aspect="equal", xlabel="x / box", ylabel="y / box", title=f"t = {t:.3g}")
    cut = slab & (np.abs(c[:, 1]) < 1.5 * h0) & (q > 0)
    axes[1].semilogy(c[cut, 0] / box, q[cut], ".", ms=2)
    axes[1].set(xlabel="x / box", ylabel="PE specific energy", title="cut through both stars")
    fig.savefig("two_stars.png", dpi=150)
    print("two_stars.png written")


if __name__ == "__main__":
    main()
