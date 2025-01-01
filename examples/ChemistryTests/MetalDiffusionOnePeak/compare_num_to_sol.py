#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:35:53 2024

@author: darwinr
"""
import argparse
import os
import cmasher as cmr
import numpy as np
import matplotlib.pyplot as plt
import swiftsimio as sw
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid
from scipy.optimize import curve_fit
from tqdm import tqdm


# %%

def gaussian(x, y, t, q_0, x_0, y_0, kappa):
    return q_0*(2*np.pi)**(-1.5) / (2*kappa*t) * np.exp(-0.5*((x-x_0)**2 + (y - y_0)**2)/(2*kappa*t))


# %%
def parse_option():
    description = """"
Plot the Fe mass diffusion and compare it with a Gaussian
    """
    epilog = """
Examples:
--------
python3 metal_projection.py snap/snapshot_*0.hdf5 --log
python3 metal_projection.py snap/snapshot_*0.hdf5
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument("--n_px",
                        action="store",
                        type=float,
                        default=2048,
                        help="Number of pixels for the surface density plots")

    parser.add_argument("--kappa",
                        action="store",
                        type=float,
                        default=2.516846e-03,
                        help="Diffusion coefficient")

    parser.add_argument('--log', default=False, action="store_true",
                        help="Density plot in log.")

    parser.add_argument("--x_min",
                        action="store",
                        type=float,
                        default=0.0,
                        help="Minimal x-axis value. In boxsize units")

    parser.add_argument("--x_max",
                        action="store",
                        type=float,
                        default=1.0,
                        help="Maximal x-axis value. In boxsize units")

    parser.add_argument("--y_min",
                        action="store",
                        type=float,
                        default=0.0,
                        help="Minimal y-axis value. In boxsize units")

    parser.add_argument("--y_max",
                        action="store",
                        type=float,
                        default=1.0,
                        help="Maximal y-axis value. In boxsize units")

    parser.parse_args()
    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError("You need to provide one file")

    return args, files

# %% Data


args, files = parse_option()
image_resolution = args.n_px
kappa = args.kappa  # Diffusion ocefficient set in the code

# Plot params
cmap = cmr.torch                   # CMasher
x_min = args.x_min
x_max = args.x_max
y_min = args.y_min
y_max = args.y_max
do_log = args.log

for filename in tqdm(files):
    snapshot_number = int(filename.split('_')[1].split('.')[0])
    output_name = "metal_surface_density_num_vs_sol_" + str(snapshot_number)

    # Get the data
    data = sw.load(filename)
    boxsize = data.metadata.boxsize

    if hasattr(data.gas.metal_mass_fractions, 'fe'):
        data.gas.m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:
        data.gas.m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses

    # Compute projected density
    metal_density = project_gas(
        data,
        resolution=image_resolution,
        project="m_fe",
        parallel=True,
        periodic=True,
    ).T

    metal_density.convert_to_cgs()
    metal_density = metal_density.value
    metal_density_log = np.log10(metal_density)

    figsize = (10, 10)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, num=1)
    vmin = None
    vmax = None

    if do_log:
        im = ax.imshow(metal_density_log,  zorder=1,
                       extent=[0, boxsize[0].value, 0, boxsize[1].value], aspect='equal',
                       cmap=cmap,
                       vmin=vmin, vmax=vmax)
    else:
        im = ax.imshow(metal_density,  zorder=1,
                       extent=[0, boxsize[0].value, 0, boxsize[1].value], aspect='equal',
                       cmap=cmap,
                       vmin=vmin, vmax=vmax)

    # Add a colorbar for all subplots
    # Adjust these values to position the colorbar
    cbar_ax = fig.add_axes([0.125, 0.9, 0.775, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')

    if do_log:
        cbar.set_label('$\log_{10}(\Sigma$ [g/cm$^2$]$)$', labelpad=1)
    else:
        cbar.set_label('$\Sigma$ [g/cm$^2$]', labelpad=1)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.xaxis.set_ticks_position('top')

    # Adjust the top to make space for the colorbar
    plt.subplots_adjust(top=0.9)

    # %Theoretical model

    # Create the 2D grid
    x = np.linspace(0, boxsize[0].value, image_resolution)
    y = np.linspace(0, boxsize[1].value, image_resolution)
    x, y = np.meshgrid(x, y)

    # Get the position of the initial particle-delta function
    index = np.argwhere(data.gas.particle_ids == 0).flatten()[0]
    pos = data.gas.coordinates[index]
    x_0 = boxsize[0].value/2
    y_0 = boxsize[1].value/2

    # Time
    t = data.metadata.time.value

    # Flatten the x, y, and metal_density arrays
    x_flat = x.flatten()
    y_flat = y.flatten()
    density_flat = metal_density.flatten()

    # Fit only q_0
    def fit_q_0(xy, q_0):
        x, y = xy
        return gaussian(x, y, t, q_0, x_0, y_0, kappa)

    # Initial guess for q_0
    initial_guess_q_0 = np.median(density_flat)

    # Perform the fit
    popt, pcov = curve_fit(fit_q_0, (x_flat, y_flat),
                           density_flat, p0=initial_guess_q_0)

    # Extract the fitted q_0 value
    q_0_fit = popt[0]

    metal_density_sol = gaussian(x, y, t, q_0_fit, x_0, y_0, kappa)
    metal_density_sol_log = np.log10(metal_density_sol)

    if do_log:
        ax.contour(metal_density_sol_log,  zorder=2,
                   extent=[0, boxsize[0].value, 0, boxsize[1].value],
                   cmap=cmap,
                   vmin=vmin, vmax=vmax)
    else:
        ax.contour(metal_density_sol,  zorder=2,
                   extent=[0, boxsize[0].value, 0, boxsize[1].value],
                   cmap=cmap,
                   vmin=vmin, vmax=vmax)

    ax.set_xlim([x_min*boxsize[0].value, x_max*boxsize[0].value])
    ax.set_ylim([y_min*boxsize[1].value, y_max*boxsize[1].value])

    ax.set_xlabel("$x$ [ckpc]")
    ax.set_ylabel("$y$ [ckpc]")

    plt.savefig(output_name+".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()
