#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:17:37 2024

@author: darwinr
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
from tqdm import tqdm
import swiftsimio as sw

# %%


def get_fe_metal_mass(data):
    """Get the Fe metal mass (in code units) from swiftsimio data"""
    # Get the metal mass
    if hasattr(data.gas.metal_mass_fractions, "fe"):
        m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:  # This case happens when we run the simulation without --feedback
        m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses
    return m_fe


def x_profile(value, x_coord, x_min=None, x_max=None, n_bins=50):
    """
    Computes a 1D profile of 'value' along the x-axis (x_coord) using linear binning.

    The number of bins is increased (default 50) as we are binning in 1D space linearly.
    """
    # Handle default x_min and x_max
    if x_min is None:
        x_min = x_coord.min()
    if x_max is None:
        x_max = x_coord.max()

    # Create LINEAR bin edges
    x_bins = np.linspace(x_min, x_max, n_bins + 1)

    # Get bin indices for each x-coordinate
    # np.digitize returns indices such that bins[i-1] <= x < bins[i]
    bin_indices = np.digitize(x_coord, bins=x_bins) - 1

    # Handle edge cases where x_coord might fall outside the range
    # np.clip handles the case where x_coord == x_max (placed in the last bin)
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)

    # Calculate the centers of the bins for plotting
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    # Initialize profile array
    values = np.zeros(n_bins)

    # Calculate the average value (e.g., metal mass) in each bin
    for i in range(n_bins):
        in_bin = bin_indices == i
        if np.any(in_bin):  # Ensure there are particles in the bin
            # Instead of np.mean(value[in_bin]), we are calculating
            # the total mass/count in the bin, and dividing by the bin width if needed.
            # For a simple 'profile' calculation like your original radial_profile,
            # we will just take the mean metal mass concentration in the x-slice.
            # If you want DENSITY (mass / volume), you would need to divide by V_bin.
            # Assuming you want the average metal mass in the x-slice:
            values[i] = np.mean(value[in_bin])

    return x_centers, values


def gaussian_1d(x, t, q_0, x_0, kappa, epsilon):
    """1D Gaussian solution for Parabolic Diffusion (or its functional form).
    q_0 is the normalization constant, x_0 is the peak center."""
    sigma_sq = epsilon**2 + 2 * kappa * t

    # Use the 1D solution form:
    return q_0 / np.sqrt(sigma_sq) * np.exp(-0.5 * ((x - x_0) ** 2) / sigma_sq)


# %%


def parse_option():
    description = """"
Plot the Fe 1D density profile
    """
    epilog = """
Examples:
--------
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --x_min 1e-1 --x_max 1.1
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --x_min 1e-1 --x_max 1.1 --log
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")

    parser.add_argument(
        "--epsilon",
        action="store",
        type=float,
        default=0.00,
        help="Size of the initial homogeneous sphere seeded with metals",
    )

    parser.add_argument(
        "--n_bins", action="store", type=int, default=40, help="Number bins"
    )

    parser.add_argument(
        "--x_min", action="store", type=float, default=1e-1, help="Minimal x."
    )

    parser.add_argument(
        "--x_max", action="store", type=float, default=1.6, help="Maximal x."
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
x_min = args.x_min
x_max = args.x_max
n_bins = args.n_bins
epsilon = args.epsilon
log = args.log

# Define the figure size
figsize = (6.4, 4.8)

# Open the data in the first snapshot to grab some information
data_init = sw.load(files[0])
boxsize = data_init.metadata.boxsize.value[0]

# Read kappa from the parameter file
try:
    kappa = float(data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
except KeyError:
    kappa = args.kappa

print(f"Using kappa = {kappa:.3e}")

for filename in tqdm(files):
    snapshot_number = int(filename.split("_")[1].split(".")[0])
    output_name = "metal_profile" + str(snapshot_number)
    data = sw.load(filename)

    if log:
        output_name = "log_" + output_name

    # Get data
    m_fe = get_fe_metal_mass(data)

    # 1. Get coordinates and center the box (coordinates in SWIFT are 0 to L)
    # The new center of the box is at 0.0
    centered_coords = data.gas.coordinates.value - 0.5 * boxsize

    # 2. Apply periodic boundary conditions (smallest image convention)
    # The centered box ranges from -L/2 to L/2.
    # We only care about the X-axis component.
    x_coord_unwrapped = centered_coords[:, 0]

    # In a periodic box, every particle's position 'x' has images at x + n*L.
    # The smallest image convention means shifting the particles so their
    # position is in the range [-L/2, L/2). Since the coordinates were already
    # centered to [-L/2, L/2), the smallest image convention is already
    # applied, but we can do a double-check wrap:
    x_coord = np.fmod(x_coord_unwrapped + 0.5 * boxsize, boxsize) - 0.5 * boxsize

    # 3. Define the x-axis range for fitting/plotting, e.g., from -L/2 to L/2
    # The range should be -0.5 * boxsize to 0.5 * boxsize
    L_half = 0.5 * boxsize

    # Use the full range of the box for the profile calculation
    x_min_profile = -L_half
    x_max_profile = L_half

    t = data.metadata.time.value

    # Compute the 1D profile along the x-axis
    x_centers, fe_bin = x_profile(
        m_fe, x_coord, x_min=x_min_profile, x_max=x_max_profile, n_bins=n_bins
    )

    # The delta function peaks are at L/4 and 3L/4 in the [0, L] box.
    # In the centered [-L/2, L/2] box, they are at:
    # x_1 = L/4 - L/2 = -L/4
    # x_2 = 3L/4 - L/2 = L/4
    x_1 = -L_half / 2.0  # corresponds to L/4 in the [0,L] box
    x_2 = L_half / 2.0  # corresponds to 3L/4 in the [0,L] box

    x_sol = np.linspace(x_min_profile, x_max_profile, 100)

    # Perform the fit on the BINNED data (x_centers and fe_bin)
    def fit_q(x, q_1, q_2):
        # Use the 1D Gaussian solution
        return gaussian_1d(x, t, q_1, x_1, kappa, epsilon) + gaussian_1d(
            x, t, q_2, x_2, kappa, epsilon
        )

    initial_guess_q_1 = 1e-1
    initial_guess_q_2 = 1e-1

    popt, pcov = curve_fit(
        fit_q, x_centers, fe_bin, p0=[initial_guess_q_1, initial_guess_q_2]
    )

    # Extract the fitted q_0 value
    q_1 = popt[0]
    q_2 = popt[1]

    # Compute the analytical solution (using the new 1D function)
    fe_1 = gaussian_1d(x_sol, t, q_1, x_1, kappa, epsilon)
    fe_2 = gaussian_1d(x_sol, t, q_2, x_2, kappa, epsilon)
    fe_sol = fe_1 + fe_2

    # Now plot
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1, figsize=figsize, layout="tight")
    ax.clear()

    # Plot the binned x-profile
    ax.plot(x_centers, fe_bin, label="Fe X-profile (Binned Mean Mass)")

    # Plot the analytical fit
    ax.plot(x_sol, fe_sol, label="Parabolic diffusion solution (on X-axis)")

    # Add markers for the analytical peak positions (L/4 and 3L/4)
    ax.axvline(
        x=x_1, color="gray", linestyle="--", linewidth=0.8, alpha=0.7, label=r"$x=L/4$"
    )
    ax.axvline(
        x=x_2, color="gray", linestyle=":", linewidth=0.8, alpha=0.7, label=r"$x=3L/4$"
    )

    ax.set_xlabel("$x$ [kpc]")  # Change label to x
    ax.set_ylabel(r"$Fe$ [M$_\odot$] (Mean Mass per Slice)")  # Update Y-label
    ax.set_xlim(x_min_profile, x_max_profile)  # Set X-limit to the full box

    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name + ".png", format="png", bbox_inches="tight", dpi=300)
    plt.close()
