#!/usr/bin/env python3
"""Plot the light of one star in a nearly transparent box: the field times r^2 against r, at four times."""

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
    files = snapshots()[1:]
    files = [files[i] for i in np.unique(np.linspace(0, len(files) - 1, 4).astype(int))]
    fig, ax = plt.subplots(figsize=(5.0, 4.0), constrained_layout=True)
    for fn in files:
        t, x, _, q = load(fn, "PacketKernelMeanSpecificEnergies")
        with h5py.File(fn, "r") as f:
            box = float(f["Header"].attrs["BoxSize"][0])
            h0 = float(f["PartType0/SmoothingLengths"][:].mean())
        r = x - 0.5 * box
        r -= box * np.round(r / box)
        r = np.linalg.norm(r, axis=1) / h0
        edges = np.arange(2.0, 0.5 * box / h0, 1.0)
        idx = np.digitize(r, edges)
        prof = np.array([np.median(q[idx == i] * r[idx == i] ** 2) if np.any(idx == i) else np.nan for i in range(1, len(edges))])
        prof[prof <= 0] = np.nan
        ref = np.nanmedian(prof[:3])
        ax.plot(0.5 * (edges[1:] + edges[:-1]), prof / ref, label=f"t = {t:.3g}")
    ax.axhline(1.0, color="k", lw=0.5, ls="--")
    ax.set(xlabel="r [h]", ylabel="field x r^2 (normalised at r < 5 h)", title="one star, nearly transparent gas")
    ax.legend()
    fig.savefig("thin_box.png", dpi=150)
    print("thin_box.png written")


if __name__ == "__main__":
    main()
