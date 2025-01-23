#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import swiftsimio as sw
import numpy as np
import matplotlib.pyplot as plt
import unyt
from tqdm import tqdm
from swiftsimio.visualisation.projection import project_gas, project_pixel_grid

# %%


def metal_surface_density(data):

    if hasattr(data.gas.metal_mass_fractions, 'fe'):
        data.gas.m_fe = data.gas.metal_mass_fractions.fe * data.gas.masses
    else:
        data.gas.m_fe = data.gas.metal_mass_fractions[:, 0] * data.gas.masses

    # Compute projected density
    projected_density = project_gas(
        data,
        resolution=image_resolution,
        project="m_fe",
        parallel=True,
        periodic=True,
    ).T

    projected_density.convert_to_cgs()
    return projected_density


def smoothed_metal_surface_density(data):

    data.gas.m_fe = data.gas.smoothed_metal_mass_fractions.fe * data.gas.masses

    # Compute projected density
    projected_density = project_gas(
        data,
        resolution=image_resolution,
        project="m_fe",
        parallel=True,
        periodic=True,
    ).T

    projected_density.convert_to_cgs()
    return projected_density


def make_plot(boxsize, metal_density, output_name, log=False):
    figsize = (10, 10)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, num=1)

    vmin = None
    vmax = None

    if log:
        metal_density = np.log10(metal_density.value)

    im = ax.imshow(metal_density, origin='lower',  zorder=1,  # interpolation='None',
                   extent=[0, boxsize[0].value, 0, boxsize[1].value], aspect='equal',
                   cmap='inferno',
                   vmin=vmin, vmax=vmax)

    ax.set_xlim([0, boxsize[0]])
    ax.set_ylim([0, boxsize[1]])

    ax.set_xlabel("$x$ [ckpc]")
    ax.set_ylabel("$y$ [ckpc]")

    # Add a colorbar for all subplots
    # Adjust these values to position the colorbar
    cbar_ax = fig.add_axes([0.125, 0.9, 0.775, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label('$\log_{10}(\Sigma$ [g/cm$^2$]$)$', labelpad=1)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.xaxis.set_ticks_position('top')

    # Adjust the top to make space for the colorbar
    plt.subplots_adjust(top=0.9)
    plt.savefig(output_name + ".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()


# %%
def parse_option():
    description = """"
Plot the Fe mass
    """
    epilog = """
Examples:
--------
python3 metal_projection.py snap/snapshot_0000.hdf5 --n_px 2048
python3 metal_projection.py snap/snapshot_*.hdf5 --log
"""
    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files",
                        nargs="+",
                        type=str,
                        help="File name(s).")

    parser.add_argument("--n_px",
                        action="store",
                        type=float,
                        default=1024,
                        help="Number of pixels for the surface density plots")

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


args, files = parse_option()
image_resolution = args.n_px
log = args.log


for filename in tqdm(files):
    # Load data
    data = sw.load(filename)
    boxsize = data.metadata.boxsize

    # Plot with the metal density
    snapshot_number = int(filename.split('_')[1].split('.')[0])
    output_name = "metal_surface_density_" + str(snapshot_number)

    if log:
        output_name = "log_" + output_name

    projected_density = metal_surface_density(data)
    make_plot(boxsize, projected_density, output_name, log)
