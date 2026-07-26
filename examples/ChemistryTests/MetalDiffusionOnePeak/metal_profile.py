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
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv as modified_bessel
import argparse
from tqdm import tqdm
import swiftsimio as sw

# %%


def get_fe_metal_mass(data):
    """Get the Fe metal mass (in code units) from swiftsimio data"""
    # Get the metal mass
    if hasattr(data.gas.metal_mass_fractions, 'fe'):
        m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else: # This case happens when we run the simulation without --feedback
        m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses
    return m_fe

def radial_profile(value, r, r_max, cross_section_area, n_bins=30):
    """
    Compute a 1D density profile (mass / bin volume) across -r_max to
    r_max. Summing mass and dividing by bin volume (rather than averaging
    per-particle values) makes the result resolution-independent -- a
    per-particle mean shrinks by construction as particle count grows,
    even when the underlying physical density is converged.
    """
    # Define a symmetric linear range
    r_min = -r_max

    # Create linear bin edges
    r_bins = np.linspace(r_min, r_max, n_bins + 1)
    bin_width = (r_max - r_min) / n_bins
    bin_volume = bin_width * cross_section_area

    # Get bin indices
    bin_indices = np.digitize(r, bins=r_bins) - 1

    # Calculate centers
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    densities = np.zeros(n_bins)
    for i in range(n_bins):
        in_bin = (bin_indices == i)
        if np.any(in_bin):
            densities[i] = np.sum(value[in_bin]) / bin_volume

    return r_centers, densities

def gaussian(r, t, q_0, r_0, kappa, epsilon):
    return q_0*(2*np.pi)**(-1.5) / (epsilon**2 + 2*kappa*t)**1.5 * np.exp(-0.5*((r-r_0)**2)/(epsilon**2 + 2*kappa*t))


def hyperbolic_diffusion_solution(x, t, q_0, x_0, tau, kappa):
    """
    Compute the solution u(x, t) of the hyperbolic diffusion equation.

    Parameters:
        x (array): Spatial positions.
        t (float): Time.
        tau (float): Relaxation time constant.
        kappa (float): Diffusion coefficient.

    Returns:
        u (array): Solution values at positions x and time t.
    """
    # Wave speed
    c = np.sqrt(kappa / tau)

    Delta_x = x - x_0

    # Exponential decay factor
    decay_factor = q_0*np.exp(-c**2 * t / (2 * kappa))

    # Compute radial term and modified Bessel functions for valid x values
    within_causal_region = np.abs(Delta_x) <= c * t
    radial_term = np.zeros_like(Delta_x)
    radial_term[within_causal_region] = np.sqrt(
        c**2 * t**2 - Delta_x[within_causal_region]**2)

    I0 = modified_bessel(0, (c / (2 * kappa)) *
                         radial_term[within_causal_region])
    I1 = modified_bessel(1, (c / (2 * kappa)) *
                         radial_term[within_causal_region])

    # Compute the solution
    u = np.zeros_like(Delta_x)
    if t > 0:
        u[within_causal_region] = 0.5 * decay_factor * (
            (c / (2 * kappa)) * I0 +
            (c**2 / (2 * kappa)) * t * I1 / radial_term[within_causal_region]
        )

    return u

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

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument("--epsilon",
                        action="store",
                        type=float,
                        default=0.04,
                        help="Size of the initial homogeneous sphere seeded with metals")

    parser.add_argument("--n_bins",
                        action="store",
                        type=int,
                        default=40,
                        help="Number bins")

    parser.add_argument("--r_max",
                        action="store",
                        type=float,
                        default=None,
                        help="Maximal r. Defaults to the box's actual half-extent.")

    parser.add_argument('--log', default=False, action="store_true",
                        help="Density plot in log.")

    parser.add_argument('--per_snapshot_q0', default=False, action="store_true",
                        help="Refit q_0 independently for every snapshot instead of "
                             "once from an early reference snapshot. Only for "
                             "debugging/comparison -- biases cross-run comparisons.")

    parser.add_argument('--track_edge', default=False, action="store_true",
                        help="Instead of per-snapshot profile plots, track the "
                             "outer edge of the populated (above --threshold) "
                             "region vs time and compare to the causal front "
                             "c_hyp*t. Geometry/amplitude-independent shape "
                             "diagnostic.")

    parser.add_argument('--threshold', default=1e-4, type=float,
                        help="Fe density [Msun/kpc^3] threshold defining the "
                             "'populated' region for --track_edge.")

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
    kappa = float(
        data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
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
    r_c, dens = radial_profile(get_fe_metal_mass(data), r, r_max, cross_section_area, n_bins)
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


# Fit q_0 once from a reference snapshot (not per-file/per-run) so
# amplitude comparisons across differently-tuned runs are apples-to-apples.
# Fit against the binned *density* profile (not raw per-particle mass) so
# the fit is in the same units as both the plotted profile and what the
# analytic q(x,t) formulas represent (a density, per the code's own U=rho*Z
# convention) -- fitting a density function to raw per-particle mass values
# would be off by the bin-volume normalization. The parabolic and
# hyperbolic curves use different reference times (see
# find_hyperbolic_reference_file): the parabolic solution has no sharp
# causal edge, so the early, cross-run-stable snapshot works directly.
FIT_N_BINS = 200
q_0_fixed = q_0_hyperbolic_fixed = None
if not args.per_snapshot_q0:
    ref_file = find_q0_reference_file(files)
    ref_t, ref_r_centers, ref_density = density_at(ref_file, FIT_N_BINS)
    q_0_fixed = fit_q0_linear(
        gaussian(ref_r_centers, ref_t, 1.0, r_0, kappa, epsilon), ref_density)

    hyp_ref_file = find_hyperbolic_reference_file(files, FIT_N_BINS)
    hyp_ref_t, hyp_ref_r_centers, hyp_ref_density = density_at(hyp_ref_file, FIT_N_BINS)
    q_0_hyperbolic_fixed = fit_q0_linear(
        hyperbolic_diffusion_solution(hyp_ref_r_centers, hyp_ref_t, 1.0, r_0, tau, kappa),
        hyp_ref_density)

    print(f"Fixed q_0: parabolic={q_0_fixed:.4e} (from {ref_file}, t={ref_t:.3f}), "
          f"hyperbolic={q_0_hyperbolic_fixed:.4e} (from {hyp_ref_file}, t={hyp_ref_t:.3f})")

edge_times, edge_radii, causal_front = [], [], []

# Now that we have all we need, do the work on everyone!
for filename in tqdm(files):
    data = sw.load(filename)

    if not args.track_edge:
        snapshot_number = int(os.path.basename(filename).split('_')[1].split('.')[0])
        output_name = "metal_profile" + str(snapshot_number)
        if log:
            output_name = "log_" + output_name

    # Get data
    m_fe = get_fe_metal_mass(data)
    t = data.metadata.time.value

    # Calculate signed distance from box center (using x-coordinate as the "radial" line)
    center = boxsize[0] / 2.0
    r_signed = data.gas.coordinates[:, 0] - center

    # Compute the profile (density: mass / bin volume, resolution-independent)
    r_centers, fe_bin = radial_profile(m_fe, r_signed, r_max, cross_section_area, n_bins=n_bins)

    if args.track_edge:
        above = np.abs(r_centers)[fe_bin > args.threshold]
        edge_times.append(t)
        edge_radii.append(above.max() if len(above) else 0.0)
        causal_front.append(np.sqrt(kappa / tau) * t)
        continue

    # Perform the fit on all the data (against the density profile, see note above)
    if args.per_snapshot_q0:
        q_0 = fit_q0_linear(
            gaussian(r_centers, t, 1.0, r_0, kappa, epsilon), fe_bin)
        q_0_hyperbolic = fit_q0_linear(
            hyperbolic_diffusion_solution(r_centers, t, 1.0, r_0, tau, kappa), fe_bin)
    else:
        q_0 = q_0_fixed
        q_0_hyperbolic = q_0_hyperbolic_fixed

    # Compute the analytical solutions
    r_sol = np.linspace(-r_max, r_max, 100)
    fe_sol = gaussian(r_sol, t, q_0, r_0, kappa, epsilon)
    fe_sol_hyperbolic = hyperbolic_diffusion_solution(
        r_sol, t, q_0_hyperbolic, r_0, tau, kappa)

    ###########
    # Now plot
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1,
                           figsize=figsize, layout="tight")
    ax.clear()
    ax.plot(r_centers, fe_bin, label='Fe mass profile')
    ax.plot(r_sol, fe_sol, label='Parabolic diffusion solution')
    ax.plot(r_sol, fe_sol_hyperbolic, label='Hyperbolic diffusion solution')
    ax.set_xlabel("$r$ [kpc]")
    ax.set_ylabel(r"$Fe$ density [M$_\odot$ kpc$^{-3}$]")
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name+".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()

if args.track_edge:
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.plot(edge_times, edge_radii, "o-", label=f"Populated edge (Fe > {args.threshold:.0e})")
    ax.plot(edge_times, causal_front, "--", label="Causal front $c_{hyp} t$")
    ax.set_xlabel("$t$")
    ax.set_ylabel("$r$ [kpc]")
    ax.legend()
    plt.savefig("metal_profile_edge.png", format='png', bbox_inches='tight', dpi=300)
    plt.close()
