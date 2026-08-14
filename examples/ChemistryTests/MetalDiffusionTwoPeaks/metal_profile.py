#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:17:37 2024

@author: darwinr
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse
from tqdm import tqdm
import swiftsimio as sw

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from chemistry_tests_common import (
    get_fe_metal_mass,
    radial_profile,
    lattice_bin_count,
    hyperbolic_diffusion_solution_convolved,
)

# %%


def periodic_hyperbolic(r, t, total_mass, epsilon, x_0, tau, kappa, L, source_shape, n_images=2):
    """Sum the exact single-seed hyperbolic response over periodic images of the source -- exact since the causal front is 0 beyond c*t until it wraps around the box."""
    u = np.zeros_like(r)
    for k in range(-n_images, n_images + 1):
        u += hyperbolic_diffusion_solution_convolved(
            r, t, total_mass, epsilon, x_0 + k * L, tau, kappa,
            with_front_term=True, source_shape=source_shape,
        )
    return u


def periodic_parabolic(r, t, total_mass, sigma2, x_0, L, n_images=2):
    """Sum the parabolic heat-kernel response over periodic images of the source."""
    u = np.zeros_like(r)
    for k in range(-n_images, n_images + 1):
        xi = x_0 + k * L
        u += np.exp(-0.5 * (r - xi) ** 2 / sigma2) / np.sqrt(2 * np.pi * sigma2)
    return total_mass * u


def find_t0_file(files):
    return min(files, key=lambda f: sw.load(f).metadata.time.value)


def wrap_centered(x, L):
    """Box-centred x-coordinate, wrapped into [-L/2, L/2)."""
    centered = x - 0.5 * L
    return np.fmod(centered + 0.5 * L, L) - 0.5 * L


# %%


def parse_option():
    description = """"
Plot the Fe 1D density profile
    """
    epilog = """
Examples:
--------
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --log
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")

    parser.add_argument(
        "--epsilon",
        action="store",
        type=float,
        default=0.04,
        help="Radius of the homogeneous sphere seeded with metal mass",
    )

    parser.add_argument(
        "--n_bins",
        action="store",
        type=int,
        default=None,
        help="Number of bins. Defaults to the count that lines up with the "
        "IC's particle lattice spacing (a mismatched fixed count aliases "
        "into spurious spikes where a bin boundary lands mid-spacing).",
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
epsilon = args.epsilon
log = args.log

# Define the figure size
figsize = (6.4, 4.8)

# Open the data in the first snapshot to grab some information
data_init = sw.load(files[0])
boxsize = data_init.metadata.boxsize
L = float(boxsize[0].value)
r_max = args.r_max if args.r_max is not None else L / 2.0
n_bins = args.n_bins
if n_bins is None:
    n_bins = lattice_bin_count(data_init.gas.coordinates.value[:, 0], L)
    print(f"Auto-selected n_bins = {n_bins} (matches the IC lattice spacing)")

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

ndim = int(data_init.metadata.dimension)
cross_section_area = float(boxsize[1].value) * float(boxsize[2].value) if ndim == 3 else 1.0
source_shape = "ball" if ndim == 3 else "segment"

# The two seeds sit at L/4 and 3L/4 in [0, L); -L/4 and +L/4 once box-centred.
x_1 = -L / 4.0
x_2 = L / 4.0

# Each seed's own total mass, measured at t=0 (before mixing): epsilon << L/2 separation, so splitting by box half is exact.
d0 = sw.load(find_t0_file(files))
m0 = get_fe_metal_mass(d0)
r0_centered = wrap_centered(d0.gas.coordinates.value[:, 0], L)
total_mass_1 = float(m0.value[r0_centered < 0.0].sum())
total_mass_2 = float(m0.value[r0_centered >= 0.0].sum())
print(f"Seed masses: M1={total_mass_1:.3e}  M2={total_mass_2:.3e}")

for filename in tqdm(files):
    snapshot_number = int(os.path.basename(filename).split("_")[1].split(".")[0])
    output_name = "metal_profile" + str(snapshot_number)
    data = sw.load(filename)

    if log:
        output_name = "log_" + output_name

    m_fe = get_fe_metal_mass(data)
    t = data.metadata.time.value
    r_signed = wrap_centered(data.gas.coordinates.value[:, 0], L)

    # Compute the profile (density: mass / bin volume, resolution-independent)
    r_centers, fe_bin = radial_profile(
        m_fe, r_signed, r_max, cross_section_area, n_bins=n_bins
    )

    # Exact self-normalised references: superpose each seed's own exact response (linearity), summed over periodic images so the fronts wrap correctly at the seam.
    r_sol = np.linspace(-r_max, r_max, 400)
    seed_var = epsilon**2 / 5 if ndim == 3 else epsilon**2 / 3
    sigma2 = 2 * kappa * t + seed_var

    fe_sol = (
        periodic_parabolic(r_sol, t, total_mass_1, sigma2, x_1, L)
        + periodic_parabolic(r_sol, t, total_mass_2, sigma2, x_2, L)
    ) / cross_section_area
    fe_sol_hyperbolic = (
        periodic_hyperbolic(r_sol, t, total_mass_1, epsilon, x_1, tau, kappa, L, source_shape)
        + periodic_hyperbolic(r_sol, t, total_mass_2, epsilon, x_2, tau, kappa, L, source_shape)
    ) / cross_section_area

    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.clear()

    ax.plot(r_centers, fe_bin, label="Fe mass profile")
    ax.plot(r_sol, fe_sol, label="Exact parabolic (heat-kernel marginal)")
    ax.plot(r_sol, fe_sol_hyperbolic, label="Exact hyperbolic (front+interior marginal)")

    ax.axvline(
        x=x_1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, label=r"$x=L/4$"
    )
    ax.axvline(
        x=x_2, color="gray", linestyle=":", linewidth=0.8, alpha=0.7, label=r"$x=3L/4$"
    )

    ax.set_xlabel("$x$ [kpc]")
    ax.set_ylabel(r"$Fe$ density [M$_\odot$ kpc$^{-3}$]")
    ax.set_xlim(-r_max, r_max)
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name + ".png", format="png", bbox_inches="tight", dpi=300)
    plt.close()
