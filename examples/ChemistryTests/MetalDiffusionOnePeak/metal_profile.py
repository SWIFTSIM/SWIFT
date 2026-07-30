#!/usr/bin/env python3
# -*- coding: utf-8 -*-
################################################################################
# This file is part of SWIFT.
# Copyright (c)  2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
################################################################################
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ive as scaled_bessel
import argparse
from tqdm import tqdm
import swiftsimio as sw

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from chemistry_tests_common import (
    get_fe_metal_mass,
    radial_profile,
    hyperbolic_diffusion_solution,
    hyperbolic_diffusion_solution_convolved,
)

# %%


def gaussian(r, t, q_0, r_0, kappa, epsilon):
    return (
        q_0
        * (2 * np.pi) ** (-1.5)
        / (epsilon**2 + 2 * kappa * t) ** 1.5
        * np.exp(-0.5 * ((r - r_0) ** 2) / (epsilon**2 + 2 * kappa * t))
    )


def hyperbolic_diffusion_solution_3d(r, t, total_mass, tau, kappa):
    """Smooth interior of the exact 3D point-source solution of the Cattaneo
    equation, for our IC (U=delta^3, zero initial flux). Substituting
    U=exp(-a t) w with a=1/(2 tau) turns the radial equation into a
    Klein-Gordon equation whose radial Green's function follows from the 1D
    one by w_3D = -(1/(2 pi r)) d(w_1D)/dr; our IC is (d_t+2a) applied to
    that momentum-kick propagator. Strictly positive; the remaining mass
    sits in the ballistic front at r=c*t (hyperbolic_3d_front_mass_fraction).
    Verified against numerical Laplace inversion to ~1e-14."""
    c = np.sqrt(kappa / tau)
    a = 0.5 / tau
    r = np.atleast_1d(np.asarray(r, dtype=float))
    u = np.zeros_like(r)
    inside = r < c * t
    s = np.sqrt(np.maximum(t**2 - (r[inside] / c) ** 2, 0.0))
    z = a * s
    scale = np.exp(z - a * t)  # ive-scaled Bessels: no overflow for t >> tau
    I0 = scaled_bessel(0, z) * scale
    I1 = scaled_bessel(1, z) * scale
    u[inside] = (
        total_mass
        * (a / (4 * np.pi * c**3))
        * (t * (a * s * I0 - 2.0 * I1) / s**3 + a * I1 / s)
    )
    return u


def hyperbolic_3d_front_mass_fraction(t, tau):
    """Fraction of the metal mass carried by the 3D ballistic front at r=c*t.
    Unlike 1D's plain exp(-t/2tau), differentiating the Heaviside cutoff and
    the 1/r geometry contributes two further front deltas, giving
    exp(-a t)(1 + a t + (a t)^2/2), a=1/(2 tau). Mass budget (this plus the
    integrated smooth interior) closes to 4e-16 over tau in [0.1,100],
    t in [0.1,5]."""
    at = t * 0.5 / tau
    return np.exp(-at) * (1.0 + at + 0.5 * at**2)


def hyperbolic_diffusion_cumulative_mass_3d(r_grid, t, total_mass, epsilon, tau, kappa):
    """Cumulative Fe mass M(<r,t) for the exact 3D solution: integrate the
    smooth interior, then add the whole front mass once r passes the front.
    Point-source front sits exactly at c*t; a finite seed of radius epsilon
    smears it over [c*t-epsilon, c*t+epsilon], which is applied here as a
    linear ramp, so the curve is comparable to a simulation's cumulative
    mass but should not be read as resolving the front's internal shape."""
    c = np.sqrt(kappa / tau)
    u = hyperbolic_diffusion_solution_3d(r_grid, t, total_mass, tau, kappa)
    integrand = 4 * np.pi * r_grid**2 * u
    cum = np.concatenate(
        [[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r_grid))]
    )
    front = total_mass * hyperbolic_3d_front_mass_fraction(t, tau)
    ramp = np.clip((r_grid - (c * t - epsilon)) / (2.0 * epsilon), 0.0, 1.0)
    return cum + front * ramp


# %%


def find_q0_reference_file(files):
    """snapshot_0001 (t~0.1): fitting q_0 there instead of per-snapshot
    avoids each run's own diffusion error biasing its analytic reference."""
    candidate = os.path.join(os.path.dirname(files[0]), "snapshot_0001.hdf5")
    if os.path.exists(candidate):
        return candidate
    return min(files, key=lambda f: sw.load(f).metadata.time.value)


def parse_option():
    description = """"
Plot the Fe 1D density profile
    """
    epilog = """
Examples:
--------
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_max 1.1
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_max 1.1 --log
python3 metal_profile.py snap/snapshot_*.hdf5 --track_edge --threshold 1e-7
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")

    parser.add_argument(
        "--epsilon",
        action="store",
        type=float,
        default=0.04,
        help="Size of the initial homogeneous sphere seeded with metals",
    )

    parser.add_argument(
        "--n_bins", action="store", type=int, default=40, help="Number bins"
    )

    parser.add_argument(
        "--r_max",
        action="store",
        type=float,
        default=None,
        help="Maximal r. Defaults to the box's actual half-extent.",
    )

    parser.add_argument(
        "--log", default=False, action="store_true", help="Density plot in log."
    )

    parser.add_argument(
        "--per_snapshot_q0",
        default=False,
        action="store_true",
        help="Refit q_0 independently for every snapshot instead of "
        "once from an early reference snapshot. Only for "
        "debugging/comparison -- biases cross-run comparisons.",
    )

    parser.add_argument(
        "--track_edge",
        default=False,
        action="store_true",
        help="Instead of per-snapshot profile plots, track the "
        "outer edge of the populated (above --threshold) "
        "region vs time and compare to the causal front "
        "c_hyp*t. Geometry/amplitude-independent shape "
        "diagnostic.",
    )

    parser.add_argument(
        "--convolve_source",
        default=True,
        action="store_true",
        help="(Now always on; kept for compatibility.) The "
        "hyperbolic reference is the exact response to "
        "the finite seed (radius --epsilon), using the "
        "dimension-appropriate seed x-marginal.",
    )

    parser.add_argument(
        "--with_front_term",
        default=True,
        action="store_true",
        help="(Now always on; kept for compatibility.) The "
        "ballistic front term (Masoliver, Porra & Weiss "
        "1993) is included in the reference curve.",
    )

    parser.add_argument(
        "--threshold",
        default=1e-4,
        type=float,
        help="Fe density [Msun/kpc^3] threshold defining the "
        "'populated' region for --track_edge.",
    )

    parser.add_argument(
        "--track_moments",
        default=False,
        action="store_true",
        help="Instead of per-snapshot profile plots, track "
        "continuous, mass-weighted statistics (RMS "
        "radius, r99, r999) vs time and compare to the "
        "causal front c_hyp*t. Unlike --track_edge "
        "(a binned threshold-crossing radius, quantized "
        "to about one bin/particle spacing), these are "
        "computed directly from unbinned per-particle "
        "data, so they resolve effects smaller than one "
        "particle spacing -- use this, not --track_edge, "
        "when comparing runs whose difference is "
        "expected to be a few percent (e.g. tuning "
        "hll_riemann_solver_psi or resolution_eta).",
    )

    parser.parse_args()
    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError("You need to provide one file")

    return args, files


# %%
# Parse the arguments
args, files = parse_option()
n_bins = args.n_bins
epsilon = args.epsilon
log = args.log
r_0 = 0.0  # Peak location for the analytical solution (center of our coords)

# Define the figure size
figsize = (6.4, 4.8)

# Open the data in the first snapshot to grab some information
data_init = sw.load(files[0])
boxsize = data_init.metadata.boxsize
r_max = args.r_max if args.r_max is not None else float(boxsize[0].value) / 2.0

# Read kappa from the parameter file
try:
    kappa = float(data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
except KeyError:
    kappa = args.kappa

print(f"Using kappa = {kappa:.3e}")

# Read tau from the parameter file
try:
    tau = float(data_init.metadata.parameters["GEARChemistry:tau"])
except KeyError:
    tau = 0.001

print(f"Using tau = {tau:.3e}")
print(f"Using r_max = {r_max:.4f}")

cross_section_area = float(boxsize[1].value) * float(boxsize[2].value)
ndim = 3 if (boxsize[1].value > 1e-3 and boxsize[2].value > 1e-3) else 1


def fit_q0_linear(func_at_q0_one, density):
    """Both gaussian() and hyperbolic_diffusion_solution() are exactly
    linear in q_0 (q_0 * f(r,t)), so the least-squares q_0 has a closed
    form -- more robust than curve_fit, which can stall (flat gradient,
    wrong scale) when the causal region is narrower than the bin width."""
    denom = np.sum(func_at_q0_one**2)
    return np.sum(density * func_at_q0_one) / denom if denom > 0 else 0.0


def density_at(filename, n_bins):
    data = sw.load(filename)
    t = data.metadata.time.value
    r = data.gas.coordinates[:, 0] - boxsize[0] / 2.0
    r_c, dens = radial_profile(
        get_fe_metal_mass(data), r, r_max, cross_section_area, n_bins
    )
    return t, r_c, dens


def find_hyperbolic_reference_file(files, n_fit_bins):
    """The hyperbolic causal front (c_hyp*t) can be narrower than even a
    fine bin at very early times, especially for large tau -- pick the
    earliest snapshot where it spans enough bins to give the linear q_0
    fit an actual signal to fit, rather than always using snapshot_0001."""
    c_hyp = np.sqrt(kappa / tau)
    bin_width = 2 * r_max / n_fit_bins
    for f in sorted(files, key=lambda f: sw.load(f).metadata.time.value):
        t = sw.load(f).metadata.time.value
        if t > 0 and c_hyp * t > 5 * bin_width:
            return f
    return sorted(files, key=lambda f: sw.load(f).metadata.time.value)[-1]


# Default reference curves are the exact, self-normalised (total-mass)
# x-marginals -- no amplitude fit. This is exact for 3D runs too: the
# y,z-marginal of the 3D Green's function is the 1D Green's function
# (integrating the PDE over y,z kills the transverse Laplacian), so only
# the seed's x-marginal shape (ball vs segment) differs by dimension.
# --per_snapshot_q0 keeps the legacy amplitude-fit curves for debugging.

edge_times, edge_radii, causal_front = [], [], []
moment_times, moment_rms, moment_r4, moment_r99, moment_r999 = [], [], [], [], []

# Now that we have all we need, do the work on everyone!
for filename in tqdm(files):
    data = sw.load(filename)

    if not args.track_edge and not args.track_moments:
        snapshot_number = int(os.path.basename(filename).split("_")[1].split(".")[0])
        output_name = "metal_profile" + str(snapshot_number)
        if log:
            output_name = "log_" + output_name

    # Get data
    m_fe = get_fe_metal_mass(data)
    t = data.metadata.time.value

    # Calculate signed distance from box center (using x-coordinate as the "radial" line)
    center = boxsize[0] / 2.0
    r_signed = data.gas.coordinates[:, 0] - center

    if args.track_moments:
        # Mass-weighted moments -- genuinely continuous (a weighted sum
        # over continuous per-particle masses), unlike a binned threshold
        # (--track_edge) or a cumulative-mass percentile radius (r99/r999
        # below): on a lattice, abs(r) only takes discrete site values, so
        # a percentile-of-mass radius snaps to lattice sites and is just as
        # quantized as --track_edge. r99/r999 are kept for a coarse sanity
        # check only -- rms and r4 are the statistics to trust for small
        # (few-percent) differences between runs.
        m_fe_v = m_fe.value
        r_v = r_signed.value
        total = m_fe_v.sum()
        rms = np.sqrt(np.sum(m_fe_v * r_v**2) / total)
        r4 = (np.sum(m_fe_v * r_v**4) / total) ** 0.25
        abs_r = np.abs(r_v)
        order = np.argsort(abs_r)
        cum = np.cumsum(m_fe_v[order]) / total
        r99 = abs_r[order][np.searchsorted(cum, 0.99)]
        r999 = abs_r[order][np.searchsorted(cum, 0.999)]
        moment_times.append(t)
        moment_rms.append(rms)
        moment_r4.append(r4)
        moment_r99.append(r99)
        moment_r999.append(r999)
        tqdm.write(
            f"t={t:.3f}  total_mass={total:.6e}  rms={rms:.5e}  "
            f"r4={r4:.5e}  r99={r99:.5e}  r999={r999:.5e}"
        )
        continue

    # Compute the profile (density: mass / bin volume, resolution-independent)
    r_centers, fe_bin = radial_profile(
        m_fe, r_signed, r_max, cross_section_area, n_bins=n_bins
    )

    if args.track_edge:
        above = np.abs(r_centers)[fe_bin > args.threshold]
        edge_times.append(t)
        edge_radii.append(above.max() if len(above) else 0.0)
        causal_front.append(np.sqrt(kappa / tau) * t)
        continue

    # Compute the analytical solutions
    r_sol = np.linspace(-r_max, r_max, 400)
    total_mass = float(m_fe.value.sum())
    if args.per_snapshot_q0:
        # legacy amplitude-fit references (debugging only): interior-only
        # hyperbolic curve, 3D-kernel-shaped parabolic curve
        q_0 = fit_q0_linear(gaussian(r_centers, t, 1.0, r_0, kappa, epsilon), fe_bin)
        q_0_hyperbolic = fit_q0_linear(
            hyperbolic_diffusion_solution(r_centers, t, 1.0, r_0, tau, kappa), fe_bin
        )
        fe_sol = gaussian(r_sol, t, q_0, r_0, kappa, epsilon)
        fe_sol_hyperbolic = hyperbolic_diffusion_solution(
            r_sol, t, q_0_hyperbolic, r_0, tau, kappa
        )
        parabolic_label = "Parabolic solution (fitted amplitude)"
        hyperbolic_label = (
            "Hyperbolic solution (interior only, front "
            "term missing, fitted amplitude)"
        )
    else:
        # exact self-normalised x-marginal references; the seed's marginal
        # variance is eps^2/5 (ball) or eps^2/3 (segment)
        seed_var = epsilon**2 / 5 if ndim == 3 else epsilon**2 / 3
        sigma2 = 2 * kappa * t + seed_var
        fe_sol = (
            total_mass
            / cross_section_area
            * np.exp(-0.5 * (r_sol - r_0) ** 2 / sigma2)
            / np.sqrt(2 * np.pi * sigma2)
        )
        fe_sol_hyperbolic = (
            hyperbolic_diffusion_solution_convolved(
                r_sol,
                t,
                total_mass,
                epsilon,
                r_0,
                tau,
                kappa,
                with_front_term=True,
                source_shape="ball" if ndim == 3 else "segment",
            )
            / cross_section_area
        )
        parabolic_label = "Exact parabolic (heat-kernel marginal)"
        hyperbolic_label = "Exact hyperbolic (front+interior marginal)"

    ###########
    # Now plot
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.clear()
    ax.plot(r_centers, fe_bin, label="Fe mass profile")
    ax.plot(r_sol, fe_sol, label=parabolic_label)
    ax.plot(r_sol, fe_sol_hyperbolic, label=hyperbolic_label)
    ax.set_xlabel("$r$ [kpc]")
    ax.set_ylabel(r"$Fe$ density [M$_\odot$ kpc$^{-3}$]")
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name + ".png", format="png", bbox_inches="tight", dpi=300)
    plt.close()

if args.track_edge:
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.plot(
        edge_times,
        edge_radii,
        "o-",
        label=f"Populated edge (Fe > {args.threshold:.0e})",
    )
    ax.plot(edge_times, causal_front, "--", label="Causal front $c_{hyp} t$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$r$ [kpc]")
    ax.legend()
    plt.savefig("metal_profile_edge.png", format="png", bbox_inches="tight", dpi=300)
    plt.close()

if args.track_moments:
    causal_front_moments = [np.sqrt(kappa / tau) * t for t in moment_times]
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.plot(moment_times, moment_rms, "o-", label="RMS radius")
    ax.plot(moment_times, moment_r4, "d-", label="$\\langle r^4\\rangle^{1/4}$")
    ax.plot(moment_times, moment_r99, "s-", label="r99 (lattice-quantized)")
    ax.plot(moment_times, moment_r999, "^-", label="r999 (lattice-quantized)")
    ax.plot(moment_times, causal_front_moments, "--", label="Causal front $c_{hyp} t$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$r$ [kpc]")
    ax.legend()
    plt.savefig("metal_profile_moments.png", format="png", bbox_inches="tight", dpi=300)
    plt.close()
