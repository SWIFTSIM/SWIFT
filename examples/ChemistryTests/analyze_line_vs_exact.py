#!/usr/bin/env python3
"""x-line-density comparison of any ChemistryTests diffusion run against the
exact solution, for arbitrary seed geometry (one peak, two peaks, half-box).

The exact reference is the periodic circular convolution of the run's OWN
initial binned line density with the exact 1D Green's function of the
equation the scheme solves: the telegrapher kernel (smooth Bessel interior
plus the two e^{-at}/2 ballistic front deltas) for hyperbolic runs -- valid
for 3D runs because the y,z-marginal of the 3D solution is exactly the 1D
Green's function -- or the heat kernel for parabolic runs (detected by tau
being absent from the used parameters).

Usage: analyze_line_vs_exact.py <rundir> [<rundir> ...]
Writes <rundir>/line_profile_vs_exact.png and prints per-snapshot metrics.
"""
import os
import sys
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import swiftsimio as sw

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "MetalDiffusionOnePeak")
)
from profile_vs_exact import u_1d, metal_mass

NBIN = 512


def green_kernel_bins(t, tau, kappa, L, nbin, parabolic):
    """Mass-per-displacement-bin of the exact 1D Green's function, wrapped;
    bin j covers displacement [j*dx-dx/2, j*dx+dx/2) so the kernel is
    symmetric about zero (no half-bin shift of the reference)."""
    dx = L / nbin
    edges = (np.arange(nbin + 1) - nbin / 2 - 0.5) * dx
    fine = np.linspace(edges[0], edges[-1], 8 * nbin + 1)
    if parabolic:
        sig = np.sqrt(2 * kappa * t)
        lam = np.exp(-0.5 * (fine / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    else:
        lam = u_1d(fine, t, tau, kappa)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (lam[1:] + lam[:-1]) * np.diff(fine))])
    g = np.diff(np.interp(edges, fine, cum))
    if not parabolic:
        c = np.sqrt(kappa / tau)
        e_at = np.exp(-t / (2 * tau))
        for xf in (c * t, -c * t):
            xw = (xf - edges[0]) % L + edges[0]
            g[np.searchsorted(edges, xw, side="right") - 1] += 0.5 * e_at
    return np.roll(g, -(nbin // 2))  # bin j <- displacement j*dx


def scan_run(rundir, snap_indices=(5, 10, 20, 50)):
    snapdir = os.path.join(rundir, "snap")
    snaps = {
        int(os.path.basename(f).split("_")[1].split(".")[0]): f
        for f in sorted(glob.glob(os.path.join(snapdir, "snapshot_*.hdf5")))
    }
    if 0 not in snaps:
        print("no snapshot_0000:", rundir)
        return
    d0 = sw.load(snaps[0])
    p = d0.metadata.parameters
    kappa = float(p["GEARChemistry:diffusion_coefficient"])
    tau_raw = p.get("GEARChemistry:tau")
    parabolic = tau_raw is None
    tau = float(tau_raw) if tau_raw is not None else np.inf
    L = float(d0.metadata.boxsize.value[0])
    edges = np.linspace(0.0, L, NBIN + 1)
    dx = L / NBIN

    m0 = metal_mass(d0)
    x0 = d0.gas.coordinates.value[:, 0]
    lam0, _ = np.histogram(x0 % L, bins=edges, weights=m0)
    M = m0.sum()
    if M <= 0 or lam0.max() <= 0:
        print("no seeded metal mass, skipping:", rundir)
        return
    seeded = lam0 > 1e-8 * lam0.max()

    idxs = [i for i in snap_indices if i in snaps]
    fig, axes = plt.subplots(
        1, len(idxs), figsize=(4.2 * len(idxs), 3.8), squeeze=False
    )
    rows = []
    for ax, idx in zip(axes[0], idxs):
        d = sw.load(snaps[idx])
        t = float(d.metadata.time.value)
        m = metal_mass(d)
        x = d.gas.coordinates.value[:, 0] % L
        sim, _ = np.histogram(x, bins=edges, weights=m)
        ker = green_kernel_bins(t, tau, kappa, L, NBIN, parabolic)
        ex = np.real(np.fft.ifft(np.fft.fft(lam0) * np.fft.fft(ker)))
        # circular W1: subtract the median of the cumulative difference so
        # the value doesn't depend on where the periodic box is cut
        cdiff = np.cumsum(sim - ex)
        w1 = np.abs(cdiff - np.median(cdiff)).sum() * dx / M
        # causality: mass in bins farther than c*t + one bin from any seed
        if parabolic:
            excess = np.nan
        else:
            c = np.sqrt(kappa / tau)
            sidx = np.where(seeded)[0]
            dist = (
                np.min(
                    np.abs(
                        (np.arange(NBIN)[:, None] - sidx[None, :] + NBIN / 2) % NBIN
                        - NBIN / 2
                    ),
                    axis=1,
                )
                * dx
            )
            excess = sim[dist > c * t + dx].sum() / M
        rows.append((t, w1, excess, int((m < 0).sum())))
        ax.stairs(sim / dx, edges, lw=1.2, label="simulation")
        ax.stairs(ex / dx, edges, color="k", lw=1.2, label="exact (binned)")
        ax.set_xlabel("x [kpc]")
        ax.set_title(f"t={t:.2f}")
        if ax is axes[0][0]:
            ax.set_ylabel("Fe line density")
            ax.legend(fontsize=8)
    kind = "parabolic (heat kernel)" if parabolic else f"hyperbolic tau={tau:g}"
    fig.suptitle(
        f"{os.path.basename(os.path.dirname(rundir))}/"
        f"{os.path.basename(rundir)}   {kind}, kappa={kappa:.3e}"
    )
    plt.tight_layout()
    out = os.path.join(rundir, "line_profile_vs_exact.png")
    plt.savefig(out, dpi=125)
    plt.close(fig)
    print(f"wrote {out}")
    for t, w1, ex, nn in rows:
        print(
            f"  t={t:5.2f}  W1={w1:.4f} kpc  excess_beyond_front="
            f"{ex if np.isnan(ex) else round(ex, 4)}  n_neg={nn}"
        )


if __name__ == "__main__":
    for rundir in sys.argv[1:]:
        try:
            scan_run(rundir)
        except Exception as exc:
            print(f"FAILED {rundir}: {type(exc).__name__}: {exc}")
