#!/usr/bin/env python3
"""Plot each run's x-axis metal profile against the exact solution of the
Cattaneo (telegrapher) equation.

A 3D run's x-profile is a line density: the 3D solution marginalised over
y and z -- which IS exactly the 1D Green's function (integrating the PDE
over y,z kills the transverse Laplacian), summed over the actual seed
particles. The front appears as edge deltas of weight e^{-at}/2 per side,
not as a flat box. Run with no arguments to sweep every run under
hyperbolic.limiter_dev_*/.
"""
import os
import sys
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import swiftsimio as sw
from scipy.special import iv
from scipy.integrate import quad


def u_3d(r, t, tau, kappa):
    """Smooth interior of the exact 3D point-source solution, unit mass."""
    c = np.sqrt(kappa / tau)
    a = 0.5 / tau
    r = np.atleast_1d(np.asarray(r, float))
    u = np.zeros_like(r)
    ins = r < c * t
    s = np.sqrt(np.maximum(t**2 - (r[ins] / c) ** 2, 0.0))
    z = a * s
    u[ins] = (
        np.exp(-a * t)
        * (a / (4 * np.pi * c**3))
        * (t * (a * s * iv(0, z) - 2 * iv(1, z)) / s**3 + a * iv(1, z) / s)
    )
    return u


def u_1d(x, t, tau, kappa):
    """Smooth interior of the exact 1D point-source solution, unit mass."""
    c = np.sqrt(kappa / tau)
    a = 0.5 / tau
    x = np.atleast_1d(np.asarray(x, float))
    u = np.zeros_like(x)
    ins = np.abs(x) < c * t
    s = np.sqrt(np.maximum(t**2 - (x[ins] / c) ** 2, 0.0))
    z = a * s
    u[ins] = a * np.exp(-a * t) / (2 * c) * (iv(0, z) + t * iv(1, z) / s)
    return u


def front_fraction(t, tau, ndim):
    """Mass fraction in the ballistic front. In 3D, differentiating the
    Heaviside cutoff and the 1/r geometry adds two front deltas beyond the
    1D one, giving the (1 + at + (at)^2/2) factor."""
    at = t * 0.5 / tau
    return np.exp(-at) * (1.0 + at + 0.5 * at**2) if ndim == 3 else np.exp(-at)


def exact_line_binned(edges, t, seed_x, seed_m, tau, kappa):
    """Exact Fe mass per x-bin for any 3D or 1D run. The y,z-marginal of
    the 3D Green's function is exactly the 1D Green's function (integrating
    the PDE over y,z kills the transverse Laplacian), so both geometries
    reduce to a sum of 1D solutions over the actual seed particles; the
    front's edge deltas (weight e^{-at}/2 each) land in their bins."""
    c = np.sqrt(kappa / tau)
    e_at = np.exp(-t / (2 * tau))
    ct = c * t
    M = seed_m.sum()
    xg = np.linspace(edges[0], edges[-1], 8001)
    lam = np.zeros_like(xg)
    mass = np.zeros(len(edges) - 1)
    for x0, m in zip(seed_x, seed_m):
        w = m / M
        lam += w * u_1d(xg - x0, t, tau, kappa)
        for xf in (x0 + ct, x0 - ct):
            if edges[0] < xf < edges[-1]:
                mass[np.searchsorted(edges, xf, side="right") - 1] += 0.5 * w * e_at
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * np.diff(xg))])
    mass += np.diff(np.interp(edges, xg, cum))
    return M * mass


def metal_mass(data):
    mfe = data.gas.metal_mass_fractions
    if hasattr(mfe, "fe"):
        return (mfe.fe * data.gas.masses).value
    return (mfe[:, 0] * data.gas.masses).value


def load_seed(snapdir):
    """Actual seed particles (x-offsets, masses) and an effective radius
    from snapshot_0000 -- the seed changed (0.02 -> 0.04, and discreteness
    varies with level) partway through these runs, so never assume it."""
    f = os.path.join(snapdir, "snapshot_0000.hdf5")
    d = sw.load(f)
    m = metal_mass(d)
    box = d.metadata.boxsize.value
    pos = d.gas.coordinates.value - box / 2.0
    ndim = 3 if (box[1] > 1e-3 and box[2] > 1e-3) else 1
    r = np.linalg.norm(pos, axis=1) if ndim == 3 else np.abs(pos[:, 0])
    sel = m > 0
    spacing = box[0] / round(len(m) ** (1.0 / ndim))
    eps = max(float(r[sel].max()) + 0.5 * spacing, spacing) if sel.any() else 0.04
    return pos[sel, 0], m[sel], eps, ndim


def seed_radius(snapdir):
    return load_seed(snapdir)[2]


def plot_run(rundir, snap_indices=(5, 10, 20, 50)):
    snapdir = os.path.join(rundir, "snap")
    snaps = {
        int(os.path.basename(f).split("_")[1].split(".")[0]): f
        for f in glob.glob(os.path.join(snapdir, "snapshot_*.hdf5"))
    }
    idxs = [i for i in snap_indices if i in snaps]
    if not idxs:
        return None
    p = sw.load(snaps[idxs[0]]).metadata.parameters
    kappa = float(p["GEARChemistry:diffusion_coefficient"])
    tau = float(p["GEARChemistry:tau"])
    c = np.sqrt(kappa / tau)
    seed_x, seed_m, eps, _ = load_seed(snapdir)

    fig, axes = plt.subplots(
        1, len(idxs), figsize=(4.0 * len(idxs), 3.8), squeeze=False
    )
    for ax, idx in zip(axes[0], idxs):
        d = sw.load(snaps[idx])
        t = float(d.metadata.time.value)
        m = metal_mass(d)
        box = d.metadata.boxsize.value
        ndim = 3 if (box[1] > 1e-3 and box[2] > 1e-3) else 1
        x = d.gas.coordinates.value[:, 0] - box[0] / 2.0
        total = m.sum()
        rmax = min(max((c * t + eps) * 2.0, 3 * eps), box[0] / 2)
        npart = int((np.abs(x) < rmax).sum())
        nb = max(8, min(80, npart // 4))
        edges = np.linspace(-rmax, rmax, nb + 1)
        mb, _ = np.histogram(x, bins=edges, weights=m)
        dens = mb / np.diff(edges)
        ex_mass = exact_line_binned(edges, t, seed_x, seed_m, tau, kappa)
        ax.stairs(dens, edges, lw=1.5, label="simulation")
        ax.stairs(
            ex_mass / np.diff(edges), edges, lw=1.5, color="k", label="exact (binned)"
        )
        ax.axvline(c * t + eps, color="k", ls="--", lw=0.9)
        ax.axvline(-c * t - eps, color="k", ls="--", lw=0.9)
        ax.set_xlabel("x [kpc]")
        ax.set_title(f"t={t:.2f}")
        if ax is axes[0][0]:
            ax.set_ylabel("Fe line density")
            ax.legend(fontsize=8)
    fig.suptitle(
        f"{os.path.basename(rundir)}   "
        f"{ndim}D, tau={tau:g}, kappa={kappa:.3e}, eps={eps:.3f}"
    )
    plt.tight_layout()
    out = os.path.join(rundir, "profile_vs_exact.png")
    plt.savefig(out, dpi=115)
    plt.close(fig)
    return out


if __name__ == "__main__":
    targets = sys.argv[1:] or sorted(
        os.path.dirname(s) for s in glob.glob("hyperbolic.limiter_dev_*/*/*/snap")
    )
    for rundir in targets:
        try:
            out = plot_run(rundir)
            print(("wrote " + out) if out else ("no snapshots: " + rundir))
        except Exception as exc:
            print(f"FAILED {rundir}: {type(exc).__name__}: {exc}")
