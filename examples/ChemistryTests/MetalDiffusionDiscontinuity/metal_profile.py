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
from pNbody import *
from scipy.optimize import curve_fit
from scipy.special import erf
import argparse
from tqdm import tqdm
import swiftsimio as sw

# %%
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
    return 0.5 * (qR + qL) + 0.5 * (qR - qL) * erf((x - x0) / np.sqrt(4 * kappa * t))

# %%
def parse_option():
    parser = argparse.ArgumentParser(description="Plot the Fe 1D x-density profile")

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")
    parser.add_argument("--n_bins", type=int, default=40, help="Number of bins")
    parser.add_argument("--x_min", type=float, default=-1.0, help="Minimum x [kpc]")
    parser.add_argument("--x_max", type=float, default=1.0, help="Maximum x [kpc]")
    parser.add_argument("--log", default=False, action="store_true", help="Log scale")

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
# Calcule qR et qL une fois à partir du premier fichier
print("Reading initial conditions from:", files[0])
nb_init = Nbody(files[0], ptypes=[0])
fe_init = nb_init.metals[:, 0] * nb_init.Mass(units="Msol")
x_init = nb_init.pos[:, 0]

# Domaine total
x_domain = x_init.max() - x_init.min()
x_mid = 0.5 * (x_init.min() + x_init.max())

qR = np.max(fe_init[x_init >= x_mid])
qL = np.max(fe_init[x_init < x_mid])
print(f"Initial qL = {qL:.3e}, qR = {qR:.3e}")

# Chargement de kappa depuis le premier fichier également
data_init = sw.load(files[0])
try:
    kappa = float(data_init.metadata.parameters["GEARChemistry:diffusion_coefficient"])
except:
    kappa = 2.516846e-03  # Default fallback
print(f"Using kappa = {kappa:.3e}")

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

    nb = Nbody(filename, ptypes=[0])
    fe = nb.metals[:, 0] * nb.Mass(units="Msol")
    x = nb.pos[:, 0]  # Get x positions

    x_centers, fe_bin = x_profile(fe, x, x_min, x_max, n_bins=n_bins)

    t = nb.atime  # Time of the snapshot (scale factor)

    # Fit using the error function solution
    # def fit_func(x, qL, qR):
    #     return diffusion_erf_profile(x, t, qR, qL, kappa, x_0)

    # popt, pcov = curve_fit(fit_func, x_centers, fe_bin, p0=(qL, qR))
    # qR, qL = popt

    x_sol = np.linspace(x_min, x_max, 1000)
    fe_sol = diffusion_erf_profile(x_sol, t, qR, qL, kappa, x_0)

    # Plotting
    fig, ax = plt.subplots(figsize=figsize, layout="tight")
    ax.plot(x_centers, fe_bin, label='Fe mass profile')
    ax.plot(x_sol, fe_sol, label='Diffusion fit (erf)', linestyle='--')
    ax.set_xlabel("$x$ [kpc]")
    ax.set_ylabel("$Fe$ [M$_\odot$]")
    ax.legend()

    if log:
        ax.set_yscale("log")

    plt.savefig(output_name + ".png", format='png', bbox_inches='tight', dpi=300)
    plt.close()
