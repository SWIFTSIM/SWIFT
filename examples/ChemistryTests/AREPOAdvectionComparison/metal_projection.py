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
"""Project the surface density of a chosen GEAR metal (or total metallicity) from SWIFT snapshots, one PNG per snapshot."""

import argparse
import os
import swiftsimio as sw
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from swiftsimio.visualisation.projection import project_gas

# %%


def select_metal_mass_fraction(mass_fractions, metal):
    """Return the per-particle mass fraction column for ``metal``: 'total', a column index, or a named element."""
    if metal == "total":
        return mass_fractions[:, -1]

    try:
        index = int(metal)
        return mass_fractions[:, index]
    except ValueError:
        pass

    if not hasattr(mass_fractions, metal):
        raise ValueError(
            f"Unknown metal '{metal}': not a valid column index, 'total', "
            "or named element."
        )
    return getattr(mass_fractions, metal)


def compute_metal_surface_density_map(data, image_resolution, metal):
    """Project ``metal``'s mass and the gas mass onto a grid and divide to get the metallicity map."""
    metal_mass_fraction = select_metal_mass_fraction(
        data.gas.metal_mass_fractions, metal
    )
    data.gas.metal_masses = metal_mass_fraction * data.gas.masses

    gas_mass_map = project_gas(
        data,
        resolution=image_resolution,
        project="masses",
        parallel=True,
        periodic=True,
    ).T

    metal_mass_map = project_gas(
        data,
        resolution=image_resolution,
        project="metal_masses",
        parallel=True,
        periodic=True,
    ).T

    return metal_mass_map / gas_mass_map


def make_plot(boxsize, metal_map, output_name, log=False, vmin=None, vmax=None):
    """Plot the projected metallicity map and save it to ``output_name``.png."""
    figsize = (10, 10)
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize, num=1)

    if log:
        metal_map = np.log10(metal_map.value)
    else:
        metal_map = metal_map.value  # Just keep as numpy array

    im = ax.imshow(
        metal_map,
        origin="lower",
        extent=[0, boxsize[0].value, 0, boxsize[1].value],
        aspect="equal",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )

    ax.set_xlim([0, boxsize[0]])
    ax.set_ylim([0, boxsize[1]])

    ax.set_xlabel("$x$ [ckpc]")
    ax.set_ylabel("$y$ [ckpc]")

    # Add a colorbar for all subplots
    cbar_ax = fig.add_axes([0.125, 0.9, 0.775, 0.05])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
    label = r"$\log_{10}(Z)$" if log else r"$Z$"
    cbar.set_label(label, labelpad=1)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("top")

    plt.subplots_adjust(top=0.9)
    plt.savefig(output_name + ".png", format="png", bbox_inches="tight", dpi=300)
    plt.close()


def parse_option():
    """Parse and validate the command-line arguments."""
    description = """Plot the metal mass surface density."""
    epilog = """Examples:
    python3 metal_projection.py snap/snapshot_0000.hdf5 --n_px 2048
    python3 metal_projection.py snap/snapshot_*.hdf5 --log --vmin -5 --vmax -2
    python3 metal_projection.py snap/snapshot_*.hdf5 --metal total"""

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument("files", nargs="+", type=str, help="File name(s).")
    parser.add_argument(
        "--n_px",
        type=int,
        default=1024,
        help="Number of pixels for the surface density plots",
    )
    parser.add_argument(
        "--metal",
        type=str,
        default="0",
        help="Metal to project: a column index (0-based), "
        "'total' for the total metallicity (last "
        "column), or a named element (e.g. 'fe') if "
        "element names were read from the parameter "
        "file. Default: '0'.",
    )
    parser.add_argument(
        "--log", default=False, action="store_true", help="Density plot in log scale"
    )
    parser.add_argument(
        "--vmin", type=float, default=None, help="Minimum value for colorbar scale"
    )
    parser.add_argument(
        "--vmax", type=float, default=None, help="Maximum value for colorbar scale"
    )

    args = parser.parse_args()
    files = args.files

    for f in files:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File not found: {f}")

    return args, files


if __name__ == "__main__":
    args, files = parse_option()
    image_resolution = args.n_px
    metal = args.metal
    log = args.log
    vmin = args.vmin
    vmax = args.vmax

    for filename in tqdm(files):
        data = sw.load(filename)
        boxsize = data.metadata.boxsize

        snapshot_number = int(os.path.basename(filename).split("_")[1].split(".")[0])
        output_name = f"{'log_' if log else ''}metal_surface_density_{snapshot_number}"

        metal_map = compute_metal_surface_density_map(data, image_resolution, metal)

        # Auto-scale per snapshot when not explicitly requested: a fixed
        # vmin/vmax saturates badly since the peak density can vary by orders
        # of magnitude across snapshots/configurations (tau, epsilon, ...).
        positive = metal_map.value[metal_map.value > 0]
        if log:
            frame_vmax = np.log10(positive.max()) if positive.size else 0.0
            frame_vmin = frame_vmax - 6.0
        else:
            frame_vmax = positive.max() if positive.size else 1.0
            frame_vmin = 0.0
        this_vmin = vmin if vmin is not None else frame_vmin
        this_vmax = vmax if vmax is not None else frame_vmax

        make_plot(
            boxsize, metal_map, output_name, log=log, vmin=this_vmin, vmax=this_vmax
        )
