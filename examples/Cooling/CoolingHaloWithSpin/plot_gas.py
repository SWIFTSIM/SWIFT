#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Mladen Ivkovic (mivkov@protonmail.com)
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
##############################################################################


import gc
import os
import argparse

import matplotlib as mpl
import swiftsimio
import swiftsimio.visualisation.rotation
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# -----------------------------------------------------------------------
# Parameters users should/may tweak

snapshot_base = "CoolingHalo"  # snapshot basename

# parameters for imshow plots
imshow_kwargs = {"origin": "lower", "cmap": "viridis"}

# parameters for swiftsimio projection
projection_kwargs = {"resolution": 512, "parallel": True, "periodic": False}
# use periodic=False here to be able to play with rotations

# parameters for scatter plots
scatterkwargs = {"s": 1, "marker": ".", "color": "black"}

# -----------------------------------------------------------------------


mpl.rcParams["text.usetex"] = True


def get_snapshot_list(snapshot_basename=snapshot_base):
    """
    Find the snapshot(s) in this directory that are to be plotted and return
    their names as list.
    """

    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename + "_") and f.endswith("hdf5"):
            snaplist.append(f)

    snaplist = sorted(snaplist)

    return snaplist


def set_colorbar(ax, im):
    """
    Set a nice colorbar for an imshow plot.
    """

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_gas(filename, fullbox=False, extent=0.1, scatter=False):
    """
    Create the actual plot.

    filename:   file to work with

    fullbox:    If True, plot contents of full box. Otherwise, plot only small
                region around center.

    extent:     if fullbox=False, plot this fraction of the boxsize left and
                right of the center of the box instead.

    scatter:    if True, make a scatterplot of particles instead of a projection.
    """

    print("working on", filename)

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata

    # select region to plot
    region = [
        (0.5 - extent) * meta.boxsize[0],
        (0.5 + extent) * meta.boxsize[0],
        (0.5 - extent) * meta.boxsize[1],
        (0.5 + extent) * meta.boxsize[1],
        (0.5 - extent) * meta.boxsize[2],
        (0.5 + extent) * meta.boxsize[2],
    ]
    imshow_extent = [region[0].v, region[1].v, region[2].v, region[3].v]

    if fullbox:
        region = None
        imshow_extent = [0, meta.boxsize[0].v, 0, meta.boxsize[1].v]

    xlabel_units_str = meta.boxsize.units.latex_representation()
    time_units_str = meta.time.units.latex_representation()

    # update shared parameters
    global imshow_kwargs
    imshow_kwargs["extent"] = imshow_extent

    global projection_kwargs
    projection_kwargs["region"] = region

    # Make the actual plot
    fig = plt.figure(figsize=(10, 5), dpi=200)
    suffix = ".png"
    if scatter:
        suffix = "-Scatter.png"
    figname = filename[:-5] + suffix

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    if scatter:
        coords = data.gas.coordinates
        if fullbox:
            mask = np.ones(coords.shape[0], dtype=bool)
        else:
            mask = coords[:, 0] >= region[0]
            mask = np.logical_and(mask, coords[:, 0] <= region[1])
            mask = np.logical_and(mask, coords[:, 1] >= region[2])
            mask = np.logical_and(mask, coords[:, 1] <= region[3])
            mask = np.logical_and(mask, coords[:, 2] >= region[4])
            mask = np.logical_and(mask, coords[:, 2] <= region[5])

        ax1.scatter(coords[mask, 0], coords[mask, 1], **scatterkwargs)
        ax2.scatter(coords[mask, 0], coords[mask, 2], **scatterkwargs)

    else:
        rotation_center = 0.5 * meta.boxsize
        # get rotation matrix from a vector. [0, 0, 1] gets you the default projection.
        rotation_matrix = swiftsimio.visualisation.rotation.rotation_matrix_from_vector(
            swiftsimio.cosmo_array([0.0, 1.0, 0.0])
        )

        mass_map_XY = swiftsimio.visualisation.project_gas(
            data, project="masses", **projection_kwargs
        )
        mass_map_YZ = swiftsimio.visualisation.project_gas(
            data,
            project="masses",
            rotation_matrix=rotation_matrix,
            rotation_center=rotation_center,
            **projection_kwargs,
        )

        im1 = ax1.imshow(mass_map_XY.T, **imshow_kwargs, norm=LogNorm())
        set_colorbar(ax1, im1)
        im2 = ax2.imshow(mass_map_YZ.T, **imshow_kwargs, norm=LogNorm())
        set_colorbar(ax2, im2)

    # Cosmetics
    ax1.set_xlabel("x [" + xlabel_units_str + "]")
    ax1.set_ylabel("y [" + xlabel_units_str + "]")
    ax2.set_xlabel("x [" + xlabel_units_str + "]")
    ax2.set_ylabel("z [" + xlabel_units_str + "]")

    # Add title
    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    #  if meta.cosmology is not None:
    #      title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time.to("Gyr"))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    print("Saved", figname)
    gc.collect()

    return


if __name__ == "__main__":

    # Setup arguments and parse them.
    parser = argparse.ArgumentParser(
        description="""
    Plots halo disk face-on and edge-on.
    Provide snapshot number as cmdline arg to plot a single snapshot.
    Alternatively, this script will plot all snapshots available in this dir.
    """
    )

    parser.add_argument(
        "snap", nargs="?", help="snapshot file name or integer index", default=None
    )
    parser.add_argument(
        "-s", "--scatter", action="store_true", help="make a scatter plot instead"
    )
    parser.add_argument(
        "-b",
        "--basename",
        nargs="?",
        help="base name of outputs",
        default=snapshot_base,
    )
    parser.add_argument(
        "-f",
        "--full",
        action="store_true",
        help="Plot full box instead of just center +- extent",
    )
    parser.add_argument(
        "-e",
        "--extent",
        nargs=1,
        help="Extent to plot around the center of the box (in units of boxsize)",
        default=0.01,
        type=float,
    )

    args = parser.parse_args()
    snapshot_base = args.basename
    try:
        # if arg is provided, argparse will return it as a list.
        extent = float(args.extent[0])
    except TypeError:
        # if arg isn't provided, argparse will return the default.
        extent = args.extent

    plot_all = True
    snapshot = ""

    # Are we plotting one snapshot, or all available?
    if args.snap is not None:
        try:
            snapint = int(args.snap)
            snapshot = snapshot_base + f"_{snapint:04d}.hdf5"
            if not os.path.exists(snapshot):
                raise FileNotFoundError(f"Snapshot '{snapshot}' not found.")
        except ValueError:
            if not os.path.exists(args.snap):
                raise ValueError(
                    "Unknown provided argument for snapshot: Must either be snapshot index integer or file name. It's neither."
                )
            snapshot = args.snap
        plot_all = False

    # get all files to be plotted into a list
    if plot_all:
        snaplist = get_snapshot_list(snapshot_base)
    else:
        snaplist = [snapshot]

    # now plot them
    for f in snaplist:
        plot_gas(f, fullbox=args.full, extent=extent, scatter=args.scatter)
