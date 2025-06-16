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
from scipy.special import erf
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

def x_profile(value, x, x_min=None, x_max=None, n_bins=30):
    if x_min is None:
        x_min = x.min()
    if x_max is None:
        x_max = x.max()

    x_bins = np.linspace(x_min, x_max, n_bins + 1)
    bin_indices = np.digitize(x, bins=x_bins) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])

    values = np.zeros(n_bins)
    for i in range(n_bins):
        in_bin = (bin_indices == i)
        if np.any(in_bin):
            values[i] = np.mean(value[in_bin])

    return x_centers, values

# Error function diffusion profile


def diffusion_erf_profile(x, t, qR, qL, kappa, x0):
    """Analytical solution to the discontinuity diffusion problem"""
    return 0.5 * (qR + qL) + 0.5 * (qR - qL) * erf((x - x0) / np.sqrt(4 * kappa * t))

# %%


def parse_option():
    parser = argparse.ArgumentParser(
        description="Plot the Fe 1D x-density profile")

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")
    parser.add_argument("--n_bins", type=int, default=40,
                        help="Number of bins")
    parser.add_argument("--x_min", type=float, default=-
                        1.0, help="Minimum x [kpc]")
    parser.add_argument("--x_max", type=float, default=1.0,
                        help="Maximum x [kpc]")
    parser.add_argument("--log", default=False,
                        action="store_true", help="Log scale")

    args = parser.parse_args()

    for f in args.files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found: {f}")

    return args, args.files


# %%
args, files = parse_option()
x_min = args.x_min
x_max = args.x_max
n_bins = args.n_bins
log = args.log
figsize = (6.4, 4.8)

#####################
print("Reading initial conditions from:", files[0])

# Open the data in the first snapshot to grab some information
data_init = sw.load(files[0])

# Read kappa from the parameter file
try:
    kappa = float(
        data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
except KeyError:
    kappa = args.kappa

# Get the data
m_fe_init = get_fe_metal_mass(data_init)
x_init = data_init.gas.coordinates[:, 0]

# Get the values. Here we assume that the the min is left and max is right.
qR = np.max(m_fe_init)
qL = np.min(m_fe_init)
print(f"Initial qL = {qL:.3e}, qR = {qR:.3e}")

# Get the middle of the domain
L = data_init.metadata.boxsize[0].value
x_0 = L/2

#####################

for filename in tqdm(files):
    snapshot_number = int(filename.split('_')[1].split('.')[0])
    output_name = f"metal_x_profile{snapshot_number}"
    data = sw.load(filename)

    if log:
        output_name = "log_" + output_name

    # Get data
    m_fe = get_fe_metal_mass(data)
    x = data.gas.coordinates[:, 0]
    t = data.metadata.time.value

    # Do the x-axis profile
    x_centers, fe_bin = x_profile(m_fe, x, x_min, x_max, n_bins=n_bins)

    # Compute the analytical solution
    # Note: We do not need to fit since we have all parameters from the ICs.
    x_sol = np.linspace(x_min, x_max, 1000)
    m_fe_sol = diffusion_erf_profile(x_sol, t, qR, qL, kappa, x_0)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, layout="tight")
    ax.plot(x_centers, fe_bin, label='Fe mass profile')
    ax.plot(x_sol, m_fe_sol, label='Diffusion fit (erf)', linestyle='--')
    ax.set_xlabel("$x$ [kpc]")
    ax.set_ylabel(r"$Fe$ [M$_\odot$]")
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name + ".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()
