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


def make_plot(boxsize, metal_density, output_name, log=False, vmin=None, vmax=None):
    figsize = (10, 10)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, num=1)

    if log:
        metal_density = np.log10(metal_density.value)
    else:
        metal_density = metal_density.value  # Just keep as numpy array

    im = ax.imshow(metal_density, origin='lower',
                   extent=[0, boxsize[0].value, 0, boxsize[1].value],
                   aspect='equal',
                   cmap='inferno',
                   vmin=vmin, vmax=vmax)

    ax.set_xlim([0, boxsize[0]])
    ax.set_ylim([0, boxsize[1]])

    ax.set_xlabel("$x$ [ckpc]")
    ax.set_ylabel("$y$ [ckpc]")

    # Add a colorbar for all subplots
    cbar_ax = fig.add_axes([0.125, 0.9, 0.775, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    label = '$\log_{10}(\Sigma$ [g/cm$^2$]$)$' if log else '$\Sigma$ [g/cm$^2$]'
    cbar.set_label(label, labelpad=1)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.xaxis.set_ticks_position('top')

    plt.subplots_adjust(top=0.9)
    plt.savefig(output_name + ".png", format='png',
                bbox_inches='tight', dpi=300)
    plt.close()


def parse_option():
    description = """Plot the Fe mass surface density."""
    epilog = """Examples:
    python3 metal_projection.py snap/snapshot_0000.hdf5 --n_px 2048
    python3 metal_projection.py snap/snapshot_*.hdf5 --log --vmin -5 --vmax -2"""

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")
    parser.add_argument("--n_px", type=float, default=1024,
                        help="Number of pixels for the surface density plots")
    parser.add_argument("--log", default=False, action="store_true",
                        help="Density plot in log scale")
    parser.add_argument("--vmin", type=float, default=None,
                        help="Minimum value for colorbar scale")
    parser.add_argument("--vmax", type=float, default=None,
                        help="Maximum value for colorbar scale")

    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found: {f}")

    return args, files


# Main script
args, files = parse_option()
image_resolution = args.n_px
log = args.log
vmin = args.vmin
vmax = args.vmax

for filename in tqdm(files):
    data = sw.load(filename)
    boxsize = data.metadata.boxsize

    snapshot_number = int(filename.split('_')[1].split('.')[0])
    output_name = f"{'log_' if log else ''}metal_surface_density_{
        snapshot_number}"

    projected_density = metal_surface_density(data)
    make_plot(boxsize, projected_density, output_name,
              log=log, vmin=vmin, vmax=vmax)
