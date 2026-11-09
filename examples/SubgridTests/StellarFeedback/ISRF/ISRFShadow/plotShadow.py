#!/usr/bin/env python3
"""Plot the ratio of the PE field with the clump to the field without it, in the plane of the stars and the clump."""

import glob
import json
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    """Return the time, and the gas ids, positions and PE-band specific energy."""
    with h5py.File(path, "r") as f:
        t = float(np.ravel(f["Header"].attrs["Time"])[0])
        x = f["PartType0/Coordinates"][:]
        ids = f["PartType0/ParticleIDs"][:]
        q = np.ravel(f["PartType0/PESpecificEnergies"][:])
    return t, x, ids, q


def main():
    geo = json.load(open("geometry.json"))
    h0 = geo["h0_kpc"]
    centre = np.array(geo["centre_kpc"])
    files = sorted(glob.glob("snap/snapshot_*.hdf5"))
    if not files or not glob.glob("../run_clump0/snap/snapshot_*.hdf5"):
        sys.exit("Need the snapshots of this run and of ../run_clump0.")
    t, x, ids, q = load(files[-1])
    _, _, ids0, q0 = load("../run_clump0/snap/" + files[-1].split("/")[-1])
    _, i1, i0 = np.intersect1d(
        ids, ids0, return_indices=True
    )  # gas that left the box is dropped
    r = (x[i1] - centre) / h0
    ratio = q[i1] / np.where(q0[i0] > 0, q0[i0], np.nan)
    keep = np.abs(r[:, 2]) < 1.5
    fig, ax = plt.subplots(figsize=(5.4, 4.6), constrained_layout=True)
    im = ax.scatter(
        r[keep, 0],
        r[keep, 1],
        c=ratio[keep],
        s=2,
        marker="s",
        vmin=0,
        vmax=1.2,
        cmap="viridis",
        rasterized=True,
    )
    fig.colorbar(im, ax=ax, label="field with clump / without")
    ax.add_patch(plt.Circle((0, 0), geo["clump_radius_h"], fill=False, ec="w", ls="--"))
    for s in (np.array(geo["stars_kpc"]) - centre) / h0:
        ax.plot(s[0], s[1], "w*", ms=8)
    ax.set(
        aspect="equal",
        xlabel="x from the clump centre [h]",
        ylabel="y from the clump centre [h]",
        title=f"t = {t:.3g}",
    )
    fig.savefig("shadow.png", dpi=150)
    print("shadow.png written")


if __name__ == "__main__":
    main()
