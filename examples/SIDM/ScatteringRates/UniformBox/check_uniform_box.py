###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################

# Checks of the SIDM scattering rate in a uniform-density periodic box, run at
# several scale factors from the same ICs (see run.sh).
#
# The expected rate of every particle is
#
#   Gamma = rho_phys (sigma/m) <v_rel>,   rho_phys = rho_comoving / a^3,
#
# with rho_comoving = total mass / box volume and <v_rel> the mean peculiar
# relative speed to the N_VREL nearest neighbours (IC velocities).
#
# 1. Scale-factor dependence. The comoving positions and peculiar velocities
#    are the same at every a, so <Gamma_SPH> must scale exactly as a^-3 and
#    Gamma_SPH / Gamma must be the same at every a. A wrong a-factor in the
#    rate would show up as an extra power of a.
#
# 2. Density bias. The kernel density written by the code includes the
#    particle's own contribution m W(0, h). Removing it leaves the
#    neighbour-only density that enters the rate. For Poisson-distributed
#    particles with h adapted to the local particle number this is still
#    biased high (~6% for 48 neighbours); the expected bias is computed here
#    with a Monte Carlo of the ghost's h iteration. Gamma_SPH / Gamma should
#    show the same bias.
#
# Everything is in the internal units of the snapshots. The scripts assume the
# default cubic-spline SIDM kernel.
#
# Usage: python3 check_uniform_box.py snap_a*/snapshot_0000.hdf5

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.optimize import brentq
from scipy.spatial import cKDTree

parser = argparse.ArgumentParser()
parser.add_argument("snaps", nargs="+", help="snapshot_0000 of each run (one per a)")
parser.add_argument("--ics", default="uniform_box.hdf5", help="Initial conditions")
parser.add_argument("-o", "--output", default="uniform_box.png")

N_VREL = 32  # neighbours used for the measured mean relative speed
N_SAMPLE = 100_000  # particles used for the measured mean relative speed

# Tolerances of the checks
TOL_SLOPE = 1e-3  # |n + 3| for <Gamma_SPH> ~ a^n
TOL_SPREAD = 1e-3  # max/min - 1 of Gamma_SPH / Gamma over a
N_SIGMA = 3  # agreement with the expected density bias

# Colours (categorical slots 1-3) and ink
C_SPH, C_NGB, C_RATE = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED = "#3b3b38", "#8a8980"

# Cubic spline (M4) in 3D, as in kernel_sidm.h
KERNEL_GAMMA = 1.825742
KERNEL_CONSTANT = 16.0 / np.pi


# ------------------ Kernel ------------------
def kernel_w(q):
    """Cubic spline W(r, h) h^3 as a function of q = r / h"""
    x = q / KERNEL_GAMMA
    w = np.where(x < 0.5, 3 * x**3 - 3 * x**2 + 0.5,
                 np.where(x < 1, (1 - x) ** 3, 0.0))
    return KERNEL_CONSTANT / KERNEL_GAMMA**3 * w


def n_neighbours(eta):
    return 4 / 3 * np.pi * (KERNEL_GAMMA * eta) ** 3


def expected_density_bias(eta, n_part=400_000, n_sample=10_000, seed=0):
    """
    Mean neighbour-only kernel density in units of the true density for
    Poisson particles, with h solved as in the SIDM ghost:
    h^3 sum_j W(r_ij, h) (including self) = eta^3. Returns mean and its error.
    """
    rng = np.random.default_rng(seed)
    pos = rng.random((n_part, 3))
    h_mean = eta * n_part ** (-1 / 3)

    # Enough neighbours to cover the kernel support for h up to 2 h_mean
    k = int(n_neighbours(2 * eta)) + 1
    dists, _ = cKDTree(pos, boxsize=1.0).query(
        pos[rng.choice(n_part, n_sample, replace=False)], k=k)

    ratio = np.empty(n_sample)
    for i, d in enumerate(dists):
        h = brentq(lambda h: kernel_w(d / h).sum() - eta**3, 0.3 * h_mean, 2 * h_mean)
        # d[0] = 0 is the particle itself; the true density is n_part
        ratio[i] = kernel_w(d[1:] / h).sum() / h**3 / n_part
    return ratio.mean(), ratio.std(ddof=1) / np.sqrt(n_sample)


def check_cubic_spline(f):
    config = f["Code"].attrs["Configuration options"].decode()
    for opt in config.split():
        opt = opt.strip("'")
        if opt.startswith(("--with-sidm-kernel=", "--with-kernel=")):
            if opt.split("=")[1] != "cubic-spline":
                raise ValueError(f"Run used {opt}; this script assumes cubic-spline")


# ------------------ Data ------------------
def read_snapshot(fname):
    with h5py.File(fname, "r") as f:
        check_cubic_spline(f)
        u = f["Units"].attrs
        params = f["Parameters"].attrs
        g = f["SIDMParticles"]
        return dict(
            a=float(f["Header"].attrs["Scale-factor"][0]),
            L=f["Header"].attrs["BoxSize"][0],
            eta=float(params["BasicSIDM:resolution_eta"]),
            sigma_over_m=float(params["BasicSIDM:sigma_over_m_cgs"])
            * u["Unit mass in cgs (U_M)"][0] / u["Unit length in cgs (U_L)"][0] ** 2,
            pos=g["Coordinates"][:],
            ids=g["ParticleIDs"][:],
            mass=g["Masses"][:].astype(np.float64),
            h=g["SmoothingLengths"][:].astype(np.float64),
            rho=g["Densities"][:].astype(np.float64),  # comoving
            rate=g["Rates"][:].astype(np.float64),  # physical
        )


def mean_vrel(pos, vel, L):
    """Mean relative speed to the N_VREL nearest neighbours"""
    rng = np.random.default_rng(0)
    sample = rng.choice(len(pos), min(N_SAMPLE, len(pos)), replace=False)
    _, nb = cKDTree(pos, boxsize=L).query(pos[sample], k=N_VREL + 1)
    return np.linalg.norm(vel[nb[:, 1:]] - vel[sample, None], axis=2).mean()


def mean_err(x):
    return x.mean(), x.std(ddof=1) / np.sqrt(len(x))


def report(name, ok, detail):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    return ok


# ------------------ Main ------------------
def main():
    args = parser.parse_args()
    snaps = sorted((read_snapshot(s) for s in args.snaps), key=lambda s: s["a"])

    # IC (peculiar) velocities; the rates in snapshot_0000 are computed
    # before the first kick, so these are the velocities that were used
    with h5py.File(args.ics, "r") as f:
        ic_ids = f["PartType7"]["ParticleIDs"][:]
        ic_vel = f["PartType7"]["Velocities"][:].astype(np.float64)
    order = np.argsort(ic_ids)

    eta = snaps[0]["eta"]
    n_part = len(snaps[0]["pos"])
    print(f"N = {n_part}, eta = {eta:g} (N_ngb ~ {n_neighbours(eta):.0f})")
    print("Computing the expected neighbour-only density bias...")
    bias, bias_err = expected_density_bias(eta)

    for s in snaps:
        vel = ic_vel[order[np.searchsorted(ic_ids, s["ids"], sorter=order)]]
        vel -= vel.mean(axis=0)
        s["rho_mean"] = s["mass"].sum() / s["L"] ** 3  # comoving
        s["gamma"] = (s["rho_mean"] / s["a"] ** 3 * s["sigma_over_m"]
                      * mean_vrel(s["pos"], vel, s["L"]))
        s["ratio_rate"] = s["rate"] / s["gamma"]
        s["ratio_rho"] = s["rho"] / s["rho_mean"]
        s["ratio_rho_ngb"] = (s["rho"] - s["mass"] * kernel_w(0.0) / s["h"] ** 3) / s["rho_mean"]

    a = np.array([s["a"] for s in snaps])
    gamma = np.array([s["gamma"] for s in snaps])
    rate_mean = np.array([s["rate"].mean() for s in snaps])
    ratio, ratio_err = np.array([mean_err(s["ratio_rate"]) for s in snaps]).T

    # ------------------ Test 1: a-dependence ------------------
    print("\n1. Scale-factor dependence of the rate")
    print(f"  {'a':>6} {'z':>6} {'Gamma':>12} {'<Gamma_SPH>':>12} {'<Gamma_SPH/Gamma>':>20}")
    for s, r, re in zip(snaps, ratio, ratio_err):
        print(f"  {s['a']:6.3f} {1 / s['a'] - 1:6.2f} {s['gamma']:12.5e} "
              f"{s['rate'].mean():12.5e} {r:10.5f} +/- {re:.5f}")
    ok = True
    if len(snaps) > 1:
        slope = np.polyfit(np.log(a), np.log(rate_mean), 1)[0]
        spread = ratio.max() / ratio.min() - 1
        ok &= report("<Gamma_SPH> ~ a^n", abs(slope + 3) < TOL_SLOPE,
                     f"n = {slope:.5f} (expected -3)")
        ok &= report("Gamma_SPH / Gamma independent of a", spread < TOL_SPREAD,
                     f"max/min - 1 = {spread:.1e}")
    else:
        slope = np.nan
        print("  (only one scale factor: skipped)")

    # ------------------ Test 2: density bias (at the largest a) ------------------
    s = snaps[-1]
    rho_m, rho_e = mean_err(s["ratio_rho"])
    ngb_m, ngb_e = mean_err(s["ratio_rho_ngb"])
    rate_m, rate_e = mean_err(s["ratio_rate"])
    print(f"\n2. Density and rate bias (a = {s['a']:g}), in units of the true values")
    print(f"  rho_SPH (code output, incl. self)   : {rho_m:.4f} +/- {rho_e:.4f}")
    print(f"  rho_SPH - m W(0, h) (neighbours)    : {ngb_m:.4f} +/- {ngb_e:.4f}")
    print(f"  Gamma_SPH                           : {rate_m:.4f} +/- {rate_e:.4f}")
    print(f"  expected neighbour-only bias        : {bias:.4f} +/- {bias_err:.4f}")
    for name, m, e in (("neighbour-only density", ngb_m, ngb_e), ("rate", rate_m, rate_e)):
        tol = N_SIGMA * np.hypot(e, bias_err)
        ok &= report(f"{name} matches expected bias", abs(m - bias) < tol,
                     f"|{m:.4f} - {bias:.4f}| < {tol:.4f}")
    print(f"\n{'All checks passed' if ok else 'SOME CHECKS FAILED'}")

    # ------------------ Plot ------------------
    fig, ((ax_a, ax_r), (ax_rho, ax_hist)) = plt.subplots(2, 2, figsize=(10, 7.5))

    a_line = np.geomspace(a.min() / 1.3, a.max() * 1.3, 50)
    ax_a.loglog(a_line, gamma[-1] * (a_line / a[-1]) ** -3, color=MUTED, lw=1.2,
                label=r"$\Gamma = \bar\rho\,a^{-3}(\sigma/m)\langle v_{\rm rel}\rangle$")
    ax_a.loglog(a, rate_mean, "o", color=C_RATE, ms=8,
                label=r"$\langle\Gamma_{\rm SPH}\rangle$")
    ax_a.set_xlabel("Scale factor $a$")
    ax_a.set_ylabel("Scattering rate [internal units]")

    ax_r.axhspan(bias - bias_err, bias + bias_err, color=MUTED, alpha=0.2, lw=0)
    ax_r.axhline(bias, color=MUTED, ls=":", lw=1.2, label="Expected density bias")
    ax_r.axhline(1.0, color=INK, lw=0.8)
    ax_r.errorbar(a, ratio, yerr=ratio_err, fmt="o", color=C_RATE, ms=8,
                  label=r"$\langle\Gamma_{\rm SPH}/\Gamma\rangle$")
    ax_r.set_xscale("log")
    ax_r.set_ylim(0.95, 1.12)
    ax_r.set_xlabel("Scale factor $a$")
    ax_r.set_ylabel(r"$\Gamma_{\rm SPH} / \Gamma$")

    bins = np.linspace(0.2, 3.0, 141)
    for x, c, name in ((s["ratio_rho"], C_SPH, r"$\rho_{\rm SPH}$ (code output)"),
                       (s["ratio_rho_ngb"], C_NGB, r"$\rho_{\rm SPH} - m\,W(0,h)$")):
        ax_rho.hist(x, bins=bins, density=True, histtype="step", color=c, lw=1.6,
                    label=f"{name}, mean {x.mean():.3f}")
        ax_rho.axvline(x.mean(), color=c, ls="--", lw=1.2)
    ax_rho.set_xlabel(r"$\rho / \bar\rho$")

    ax_hist.hist(s["ratio_rate"], bins=bins, density=True, histtype="step",
                 color=C_RATE, lw=1.6,
                 label=rf"$\Gamma_{{\rm SPH}}$, mean {rate_m:.3f}")
    ax_hist.axvline(rate_m, color=C_RATE, ls="--", lw=1.2)
    ax_hist.set_xlabel(r"$\Gamma_{\rm SPH} / \Gamma$")

    for ax in (ax_rho, ax_hist):
        ax.axvline(1.0, color=INK, lw=0.8)
        ax.axvline(bias, color=MUTED, ls=":", lw=1.2,
                   label=f"Expected bias ({bias:.3f})")
        ax.set_ylabel("PDF")
        ax.set_xlim(bins[0], bins[-1])

    for ax in (ax_a, ax_r, ax_rho, ax_hist):
        ax.legend(fontsize=8, frameon=False, labelcolor=INK)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
