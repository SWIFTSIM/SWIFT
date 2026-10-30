#!/usr/bin/env python3
"""Plot the light behind the absorbing clump. With the control run (clump=0) it also plots the ratio of the two fields."""

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


def plane(geo):
    axis = np.array(geo["axis"])
    e2 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
    n = np.cross(axis, e2)
    return axis, e2, n / np.linalg.norm(n)


def main():
    geo = json.load(open("geometry.json"))
    h0 = geo["h0_kpc"]
    axis, e2, n = plane(geo)
    star = np.array(geo["star_kpc"])
    last = snapshots()[-1]
    t, x, ids, q = load(last)
    q0 = None
    try:  # control run: same gas, matched by particle id; gas that left the box is dropped
        _, _, ids0, c0 = load("../run_clump0/snap/" + last.split("/")[-1])
        common, i1, i0 = np.intersect1d(ids, ids0, return_indices=True)
        x, q, q0 = x[i1], q[i1], c0[i0]
    except OSError:
        print("No control run ../run_clump0: only the first panel is plotted.")
    r = (x - star) / h0
    dist = np.linalg.norm(r, axis=1)
    keep = (np.abs(r @ n) < 1.5) & (dist > 2.0)
    u, v = r[keep] @ axis, r[keep] @ e2
    flux = q * dist**2  # removes the 1/r^2 dilution
    ncol = 1 if q0 is None else 2
    fig, axes = plt.subplots(
        1, ncol, figsize=(4.6 * ncol, 4.2), squeeze=False, constrained_layout=True
    )
    ok = flux[keep] > 0
    im = axes[0][0].scatter(
        u[ok],
        v[ok],
        c=np.log10(flux[keep][ok]),
        s=2,
        marker="s",
        cmap="inferno",
        rasterized=True,
    )
    im.set_clim(np.log10(flux[keep][ok]).max() - 4.0, None)
    fig.colorbar(im, ax=axes[0][0], label="log10 (field x r^2), with clump")
    if q0 is not None:
        ratio = q[keep] / np.where(q0[keep] > 0, q0[keep], np.nan)
        im = axes[0][1].scatter(
            u,
            v,
            c=ratio,
            s=2,
            marker="s",
            vmin=0,
            vmax=1.2,
            cmap="viridis",
            rasterized=True,
        )
        fig.colorbar(im, ax=axes[0][1], label="field with clump / without")
    clump = (np.array(geo["clump_centre_kpc"]) - star) / h0
    for ax in axes[0]:
        ax.add_patch(
            plt.Circle(
                (clump @ axis, clump @ e2),
                geo["clump_radius_kpc"] / h0,
                fill=False,
                ec="w",
                ls="--",
            )
        )
        ax.plot(0, 0, "w*", ms=8)
        ax.set(
            aspect="equal",
            xlabel="along the star-clump axis [h]",
            ylabel="across [h]",
            title=f"t = {t:.3g}",
        )
    fig.savefig("shadow.png", dpi=150)
    print("shadow.png written")


if __name__ == "__main__":
    main()
