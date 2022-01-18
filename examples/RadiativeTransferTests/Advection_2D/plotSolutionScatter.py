#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
# plot 2D photon data, using 1D scatterplots
# give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir
# ----------------------------------------------------

import sys
import os
import swiftsimio
import numpy as np
import gc
from matplotlib import pyplot as plt
import matplotlib as mpl

# Parameters users should/may tweak
plot_all_data = True  # plot all groups and all photon quantities
snapshot_base = "output"  # snapshot basename
fancy = True  # fancy up the plots a bit?

# parameters for imshow plots

scatterplot_kwargs = {
    "alpha": 0.6,
    "s": 4,
    "marker": ".",
    "linewidth": 0.0,
    "facecolor": "blue",
}
# -----------------------------------------------------------------------


# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True


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

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def set_colorbar(ax, im):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


def plot_photons(filename, energy_boundaries=None, flux_boundaries=None):
    """
    Create the actual plot.

    filename: file to work with
    energy_boundaries:  list of [E_min, E_max] for each photon group. 
                        If none, limits are set automatically.
    flux_boundaries:    list of [F_min, F_max] for each photon group. 
                        If none, limits are set automatically.
    """

    print("working on", filename)

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
    xlabel_units_str = meta.boxsize.units.latex_representation()
    x_coordinates = data.gas.coordinates[:, 0]

    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        # add mass weights to remove surface density dependence in images
        new_attribute_str = "radiation_energy" + str(g + 1)
        en = getattr(data.gas.photon_energies, "group" + str(g + 1))
        setattr(data.gas, new_attribute_str, en)

        if plot_all_data:
            # prepare also the fluxes
            #  for direction in ["X", "Y", "Z"]:
            for direction in ["X", "Y"]:
                new_attribute_str = "radiation_flux" + str(g + 1) + direction
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                setattr(data.gas, new_attribute_str, f)

    if plot_all_data:
        fig = plt.figure(figsize=(5 * 3, 5.05 * ngroups), dpi=200)
        figname = filename[:-5] + "-scatter-all-quantities.png"

        for g in range(ngroups):

            # get energy projection
            new_attribute_str = "radiation_energy" + str(g + 1)
            energies = getattr(data.gas, new_attribute_str)

            ax = fig.add_subplot(ngroups, 3, g * 3 + 1)
            ax.scatter(x_coordinates, energies, **scatterplot_kwargs)
            ax.set_ylabel("Group {0:2d}".format(g + 1))
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Energies")

            # get flux X projection
            new_attribute_str = "radiation_flux" + str(g + 1) + "X"
            fluxX = getattr(data.gas, new_attribute_str)
            flux_units_str = fluxX.units.latex_representation()

            ax = fig.add_subplot(ngroups, 3, g * 3 + 2)
            ax.scatter(x_coordinates, fluxX, **scatterplot_kwargs)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("flux x [$" + flux_units_str + "$]")
            if g == 0:
                ax.set_title("Flux X")

            # get flux Y projection
            new_attribute_str = "radiation_flux" + str(g + 1) + "Y"
            fluxX = getattr(data.gas, new_attribute_str)
            flux_units_str = fluxX.units.latex_representation()

            ax = fig.add_subplot(ngroups, 3, g * 3 + 3)
            ax.scatter(x_coordinates, fluxX, **scatterplot_kwargs)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("flux y [$" + flux_units_str + "$]")
            if g == 0:
                ax.set_title("Flux Y")

    else:  # plot just energies

        fig = plt.figure(figsize=(5 * ngroups, 5), dpi=200)
        figname = filename[:-5] + "-scatter.png"

        for g in range(ngroups):

            # get energy projection
            new_attribute_str = "radiation_energy" + str(g + 1)
            energies = getattr(data.gas, new_attribute_str)
            energy_units_str = energies.units.latex_representation()

            ax = fig.add_subplot(ngroups, 3, g * 3 + 1)
            ax.scatter(x_coordinates, energies, **scatterplot_kwargs)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_title("Group {0:2d}".format(g + 1))
            ax.set_ylabel("Energies [$" + energy_units_str + "$]")

    # Add title
    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time)
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()

    return


def get_minmax_vals(snaplist):
    """
    Find minimal and maximal values for energy and flux
    so you can fix axes limits over all snapshots

    snaplist: list of snapshot filenames

    returns:

    energy_boundaries: list of [E_min, E_max] for each photon group
    flux_boundaries: list of [Fx_min, Fy_max] for each photon group
    """

    emins = []
    emaxs = []
    fmins = []
    fmaxs = []

    for filename in snaplist:

        data = swiftsimio.load(filename)
        meta = data.metadata

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
        emin_group = []
        emax_group = []
        fluxmin_group = []
        fluxmax_group = []

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            emin_group.append(en.min())
            emax_group.append(en.max())

            dirmin = []
            dirmax = []
            for direction in ["X", "Y"]:
                new_attribute_str = "radiation_flux" + str(g + 1) + direction
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                dirmin.append(f.min())
                dirmax.append(f.max())
            fluxmin_group.append(min(dirmin))
            fluxmax_group.append(max(dirmax))

        emins.append(emin_group)
        emaxs.append(emax_group)
        fmins.append(fluxmin_group)
        fmaxs.append(fluxmax_group)

    energy_boundaries = []
    flux_boundaries = []
    for g in range(ngroups):
        emin = min([emins[f][g] for f in range(len(snaplist))])
        emax = max([emaxs[f][g] for f in range(len(snaplist))])
        energy_boundaries.append([emin, emax])
        fmin = min([fmins[f][g] for f in range(len(snaplist))])
        fmax = max([fmaxs[f][g] for f in range(len(snaplist))])
        flux_boundaries.append([fmin, fmax])

    return energy_boundaries, flux_boundaries


if __name__ == "__main__":

    snaplist = get_snapshot_list(snapshot_base)
    if fancy:
        energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)
    else:
        energy_boundaries = None
        flux_boundaries = None

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
