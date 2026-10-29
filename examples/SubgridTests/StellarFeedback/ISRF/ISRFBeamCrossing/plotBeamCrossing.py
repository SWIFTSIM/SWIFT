#!/usr/bin/env python3
"""Plot the light of the two beams in the plane that holds both directions, at four times."""

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
    geo = json.load(open("geometry.json"))
    h0, box = geo["h0_kpc"], geo["box_kpc"]
    da, db = (np.array(b["direction"]) for b in geo["beams"])
    e1 = da
    n = np.cross(da, db)
    if np.linalg.norm(n) < 1e-6:  # collinear beams: any plane that holds the axis
        n = np.cross(da, np.eye(3)[np.argmin(np.abs(da))])
    n /= np.linalg.norm(n)
    e2 = np.cross(n, e1)
    files = snapshots()[1:]
    files = [files[i] for i in np.unique(np.linspace(0, len(files) - 1, 4).astype(int))]
    fig, axes = plt.subplots(1, len(files), figsize=(3.2 * len(files), 3.6), squeeze=False, constrained_layout=True)
    data = []
    for fn in files:
        t, x, _, q = load(fn)
        r = x - np.array(geo["crossing_kpc"])
        r -= box * np.round(r / box)
        r /= h0
        keep = np.abs(r @ n) < 1.5
        data.append((t, r[keep] @ e1, r[keep] @ e2, q[keep]))
    vmax = max(d[3].max() for d in data)
    for ax, (t, u, v, q) in zip(axes[0], data):
        im = ax.scatter(u, v, c=q, s=2, marker="s", vmin=0, vmax=vmax, cmap="inferno", rasterized=True)
        ax.set(aspect="equal", title=f"t = {t:.3g}", xlabel="along beam 1 [h]")
    axes[0][0].set_ylabel("in the beam plane [h]")
    fig.colorbar(im, ax=axes[0].tolist(), label="PE specific energy of the light (internal units)")
    fig.savefig("beam_crossing.png", dpi=150)
    print("beam_crossing.png written")


if __name__ == "__main__":
    main()
