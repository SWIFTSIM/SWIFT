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
    if hasattr(data.gas.metal_mass_fractions, 'fe'):
        m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else: # This case happens when we run the simulation without --feedback
        m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses
    return m_fe

def radial_profile(value, r, r_min=None, r_max=None, n_bins=30):
    # Handle default r_min and r_max
    if r_min is None:
        r_min = r.min()
    if r_max is None:
        r_max = r.max()

    # Create logarithmic bin edges
    r_bins = np.logspace(np.log10(r_min), np.log10(r_max), n_bins + 1)

    # Get bin indices for each radius
    bin_indices = np.digitize(r, bins=r_bins) - 1

    # Handle edge cases where radius might fall outside the range
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)

    # Calculate the centers of the bins for plotting
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

    # Initialize temperature array
    values = np.zeros(n_bins)

    # Loop through each bin and accumulate temperatures
    for i in range(n_bins):
        in_bin = (bin_indices == i)
        if np.any(in_bin):  # Ensure there are values in the bin
            values[i] = np.mean(value[in_bin])  # Average temperature in bin

    return r_centers, values


def gaussian(r, t, q_0, r_0, kappa, epsilon):
    return q_0*(2*np.pi)**(-1.5) / (epsilon**2 + 2*kappa*t)**1.5 * np.exp(-0.5*((r-r_0)**2)/(epsilon**2 + 2*kappa*t))

# %%


def parse_option():
    description = """"
Plot the Fe 1D density profile
    """
    epilog = """
Examples:
--------
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_min 1e-1 --r_max 1.1
python3 metal_profile.py snap/snapshot_*0.hdf5 --n_bins 30 --r_min 1e-1 --r_max 1.1 --log
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument("--epsilon",
                        action="store",
                        type=float,
                        default=0.00,
                        help="Size of the initial homogeneous sphere seeded with metals")

    parser.add_argument("--n_bins",
                        action="store",
                        type=int,
                        default=40,
                        help="Number bins")

    parser.add_argument("--r_min",
                        action="store",
                        type=float,
                        default=1e-1,
                        help="Minimal r.")

    parser.add_argument("--r_max",
                        action="store",
                        type=float,
                        default=1.6,
                        help="Maximal r.")

    parser.add_argument('--log', default=False, action="store_true",
                        help="Density plot in log.")

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
r_min = args.r_min
r_max = args.r_max
n_bins = args.n_bins
epsilon = args.epsilon
log = args.log

# Define the figure size
figsize = (6.4, 4.8)

# Open the data in the first snapshot to grab some information
data_init = sw.load(files[0])
boxsize = data_init.metadata.boxsize

# Read kappa from the parameter file
try:
    kappa = float(
        data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
except KeyError:
    kappa = args.kappa

print(f"Using kappa = {kappa:.3e}")

for filename in tqdm(files):
    snapshot_number = int(filename.split('_')[1].split('.')[0])
    output_name = "metal_profile" + str(snapshot_number)
    data = sw.load(filename)

    if log:
        output_name = "log_" + output_name

    # Get data
    m_fe = get_fe_metal_mass(data)
    r = np.linalg.norm(data.gas.coordinates, axis=1)
    t = data.metadata.time.value

    # Compute the radial profile
    r_centers, fe_bin = radial_profile(m_fe, r, r_min, r_max, n_bins=n_bins)

    # These values were set by hand to match the positions because this scripts
    #does not set the referential in the middle of the boxsize. Changing the
    # referential shows only one peak...
    r_1 = 0.54
    r_2 = 0.8
    r_sol = np.linspace(r_min, r_max, 100)

    # Perform the fit on all the data
    def fit_q(r, q_1, q_2):
        return gaussian(r, t, q_1, r_1, kappa, epsilon) + gaussian(r, t, q_2, r_2, kappa, epsilon)

    initial_guess_q_1 = 1e1
    initial_guess_q_2 = 1e1
    popt, pcov = curve_fit(fit_q, r, m_fe, p0=[initial_guess_q_1, initial_guess_q_2])

    # Extract the fitted q_0 value
    q_1 = popt[0]
    q_2 = popt[1]

    # Compute the analytical solution
    fe_1 = gaussian(r_sol, t, q_1, r_1, kappa, epsilon)
    fe_2 = gaussian(r_sol, t, q_2, r_2, kappa, epsilon)
    fe_sol = fe_1 + fe_2

    # Now plot
    fig, ax = plt.subplots(num=1, nrows=1, ncols=1,
                           figsize=figsize, layout="tight")
    ax.clear()
    ax.plot(r_centers, fe_bin, label='Fe mass profile')
    ax.plot(r_sol, fe_sol, label='Parabolic diffusion solution')
    ax.set_xlabel("$r$ [kpc]")
    ax.set_ylabel(r"$Fe$ [M$_\odot$]")
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name+".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()
