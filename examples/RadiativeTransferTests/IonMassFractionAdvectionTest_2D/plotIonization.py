#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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


# ----------------------------------------------------
# Plot the ionizing species
# ----------------------------------------------------

import os
import sys

import matplotlib as mpl
import swiftsimio
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

# Parameters users should/may tweak

# plot all groups and all photon quantities
plot_all_data = True
# fancy up the plots a bit?
fancy = True

# parameters for imshow plots
imshow_kwargs = {"origin": "lower", "cmap": "brg"}

# parameters for swiftsimio projections
projection_kwargs = {"resolution": 512, "parallel": True}

# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True
# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--snapshot-number", help="Number of snapshot to plot", type=int
    )
    parser.add_argument("-b", "--basename", help="Snapshot basename", default="output")

    args = parser.parse_args()
    return args


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)
        if len(snaplist) == 0:
            print(f"No snapshots with base '{snapshot_basename}' found")
            sys.exit(1)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_ionization(filename):
    """
    Create and save the plot
    """
    print("working on", filename)

    data = swiftsimio.load(filename)
    meta = data.metadata

    global imshow_kwargs
    imshow_kwargs["extent"] = [
        0.01 * meta.boxsize[0].v,
        0.99 * meta.boxsize[0].v,
        0.01 * meta.boxsize[1].v,
        0.99 * meta.boxsize[1].v,
    ]

    mass_map = swiftsimio.visualisation.projection.project_gas(
        data, project="masses", **projection_kwargs
    )

    data.gas.mXHI = data.gas.ion_mass_fractions.HI * data.gas.masses
    data.gas.mXHII = data.gas.ion_mass_fractions.HII * data.gas.masses
    data.gas.mXHeI = data.gas.ion_mass_fractions.HeI * data.gas.masses
    data.gas.mXHeII = data.gas.ion_mass_fractions.HeII * data.gas.masses
    data.gas.mXHeIII = data.gas.ion_mass_fractions.HeIII * data.gas.masses

    mass_weighted_XHImap = swiftsimio.visualisation.projection.project_gas(
        data, project="mXHI", **projection_kwargs
    )
    mass_weighted_XHIImap = swiftsimio.visualisation.projection.project_gas(
        data, project="mXHII", **projection_kwargs
    )
    mass_weighted_XHeImap = swiftsimio.visualisation.projection.project_gas(
        data, project="mXHeI", **projection_kwargs
    )
    mass_weighted_XHeIImap = swiftsimio.visualisation.projection.project_gas(
        data, project="mXHeII", **projection_kwargs
    )
    mass_weighted_XHeIIImap = swiftsimio.visualisation.projection.project_gas(
        data, project="mXHeIII", **projection_kwargs
    )

    XHImap = mass_weighted_XHImap / mass_map
    XHIImap = mass_weighted_XHIImap / mass_map
    XHeImap = mass_weighted_XHeImap / mass_map
    XHeIImap = mass_weighted_XHeIImap / mass_map
    XHeIIImap = mass_weighted_XHeIIImap / mass_map

    fig = plt.figure(figsize=(20, 8), dpi=200)
    figname = filename[:-5] + "-mass-fractions.png"

    ax1 = fig.add_subplot(251)
    ax2 = fig.add_subplot(252)
    ax3 = fig.add_subplot(253)
    ax4 = fig.add_subplot(254)
    ax5 = fig.add_subplot(255)
    ax6 = fig.add_subplot(256)
    ax7 = fig.add_subplot(257)
    ax8 = fig.add_subplot(258)
    ax9 = fig.add_subplot(259)
    ax10 = fig.add_subplot(2, 5, 10)

    #  im0 = ax0.imshow(mass_map.T, **imshow_kwargs)
    #  set_colorbar(ax0, im0)
    #  ax0.set_title("Mass")

    imshow_kwargs_ions = imshow_kwargs.copy()
    imshow_kwargs_ions["norm"] = SymLogNorm(vmin=0.0, vmax=1.0, linthresh=1e-3, base=10)

    im1 = ax1.imshow(XHImap.T, **imshow_kwargs_ions)
    set_colorbar(ax1, im1)
    ax1.set_title("HI Mass Fraction")

    im2 = ax2.imshow(XHIImap.T, **imshow_kwargs_ions)
    set_colorbar(ax2, im2)
    ax2.set_title("HII Mass Fraction")

    im3 = ax3.imshow(XHeImap.T, **imshow_kwargs_ions)
    set_colorbar(ax3, im3)
    ax3.set_title("HeI Mass Fraction")

    im4 = ax4.imshow(XHeIImap.T, **imshow_kwargs_ions)
    set_colorbar(ax4, im4)
    ax4.set_title("HeII Mass Fraction")

    im5 = ax5.imshow(XHeIIImap.T, **imshow_kwargs_ions)
    set_colorbar(ax5, im5)
    ax5.set_title("HeIII Mass Fraction")

    im6 = ax6.imshow(mass_weighted_XHImap.T, **imshow_kwargs)
    set_colorbar(ax6, im6)
    ax6.set_title("HI Mass")

    im7 = ax7.imshow(mass_weighted_XHIImap.T, **imshow_kwargs)
    set_colorbar(ax7, im7)
    ax7.set_title("HII Mass")

    im8 = ax8.imshow(mass_weighted_XHeImap.T, **imshow_kwargs)
    set_colorbar(ax8, im8)
    ax8.set_title("HeI Mass")

    im9 = ax9.imshow(mass_weighted_XHeIImap.T, **imshow_kwargs)
    set_colorbar(ax9, im9)
    ax9.set_title("HeII Mass")

    im10 = ax10.imshow(mass_weighted_XHeIIImap.T, **imshow_kwargs)
    set_colorbar(ax10, im10)
    ax10.set_title("HeIII Mass")

    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]:
        ax.set_xlabel("[kpc]")
        ax.set_ylabel("[kpc]")

    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time.to("Myr"))
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    return


if __name__ == "__main__":
    # Get command line arguments
    args = parse_args()

    if args.snapshot_number:
        plot_all = False
        snapnr = int(args.snapshot_numbeR)
    else:
        plot_all = True

    snapshot_base = args.basename
    snaplist = get_snapshot_list(snapshot_base)

    for f in snaplist:
        plot_ionization(f)
