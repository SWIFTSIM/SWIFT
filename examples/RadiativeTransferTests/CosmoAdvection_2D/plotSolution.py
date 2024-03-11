#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2024 Stan Verhoeve (s06verhoeve@gmail.com)
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
# plot photon data for 2D problems
# give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir
# ----------------------------------------------------

import gc
import os
import sys
import argparse

import unyt
import numpy as np
import matplotlib as mpl
import swiftsimio
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Parameters users should/may tweak
plot_all_data = True  # plot all groups and all photon quantities
snapshot_base = "output"  # snapshot basename
fancy = True  # fancy up the plots a bit?
plot_physical_quantities = True

# parameters for imshow plots
imshow_kwargs = {"origin": "lower", "cmap": "viridis"}


projection_kwargs = {"resolution": 1024, "parallel": True}
# -----------------------------------------------------------------------


plot_all = False
mpl.rcParams["text.usetex"] = True


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--snapshot-number", help="Number of snapshot to plot", type=int
    )
    parser.add_argument(
        "-z",
        "--redshift",
        help="Redshift domain to plot advection for",
        default="high_redshift",
    )

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
            print(f"No snapshots with base {snapshot_basename} found!")
            sys.exit(1)

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

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])
    xlabel_units_str = meta.boxsize.units.latex_representation()

    global imshow_kwargs
    imshow_kwargs["extent"] = [0, meta.boxsize[0].v, 0, meta.boxsize[1].v]

    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        # add mass weights to remove surface density dependence in images
        new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
        en = getattr(data.gas.photon_energies, "group" + str(g + 1))
        en = en * data.gas.masses
        setattr(data.gas, new_attribute_str, en)

        if plot_all_data:
            # prepare also the fluxes
            #  for direction in ["X", "Y", "Z"]:
            for direction in ["X", "Y"]:
                new_attribute_str = (
                    "mass_weighted_radiation_flux" + str(g + 1) + direction
                )
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                # f *= data.gas.masses
                f = f * data.gas.masses
                setattr(data.gas, new_attribute_str, f)

    # get mass surface density projection that we'll use to remove density dependence in image
    mass_map = swiftsimio.visualisation.projection.project_gas(
        data, project="masses", **projection_kwargs
    )

    if plot_all_data:
        fig = plt.figure(figsize=(5 * 3, 5.05 * ngroups), dpi=200)
        figname = filename[:-5] + f"-all-quantities.png"

        for g in range(ngroups):

            # get energy projection
            new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
            photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            photon_map = photon_map / mass_map

            ax = fig.add_subplot(ngroups, 3, g * 3 + 1)
            if energy_boundaries is not None:
                imshow_kwargs["vmin"] = energy_boundaries[g][0]
                imshow_kwargs["vmax"] = energy_boundaries[g][1]
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_ylabel("Group {0:2d}".format(g + 1))
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Energies")

            # get flux X projection
            new_attribute_str = "mass_weighted_radiation_flux" + str(g + 1) + "X"
            photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            photon_map = photon_map / mass_map

            ax = fig.add_subplot(ngroups, 3, g * 3 + 2)
            if flux_boundaries is not None:
                imshow_kwargs["vmin"] = flux_boundaries[g][0]
                imshow_kwargs["vmax"] = flux_boundaries[g][1]
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("y [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Flux X")

            # get flux Y projection
            new_attribute_str = "mass_weighted_radiation_flux" + str(g + 1) + "Y"
            photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            photon_map = photon_map / mass_map

            ax = fig.add_subplot(ngroups, 3, g * 3 + 3)
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_xlabel("x [$" + xlabel_units_str + "$]")
            ax.set_ylabel("y [$" + xlabel_units_str + "$]")
            if g == 0:
                ax.set_title("Flux Y")

    else:  # plot just energies

        fig = plt.figure(figsize=(5 * ngroups, 5), dpi=200)
        figname = filename[:-5] + ".png"

        for g in range(ngroups):

            # get projection
            new_attribute_str = "mass_weighted_radiation_energy" + str(g + 1)
            photon_map = swiftsimio.visualisation.projection.project_gas(
                data, project=new_attribute_str, **projection_kwargs
            )
            photon_map = photon_map / mass_map

            ax = fig.add_subplot(1, ngroups, g + 1)
            if energy_boundaries is not None:
                imshow_kwargs["vmin"] = energy_boundaries[g][0]
                imshow_kwargs["vmax"] = energy_boundaries[g][1]
            im = ax.imshow(photon_map.T, **imshow_kwargs)
            set_colorbar(ax, im)
            ax.set_title("Group {0:2d}".format(g + 1))
            if g == 0:
                ax.set_ylabel("Energies")

    # Add title
    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(1 * meta.time)
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()
    gc.collect()

    return


def get_minmax_vals(snaplist):
    """
    Find minimal and maximal values for energy and flux,
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

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])
        emin_group = []
        emax_group = []
        fluxmin_group = []
        fluxmax_group = []

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            emin_group.append((1 * en.min()).value)
            emax_group.append((1 * en.max()).value)

            dirmin = []
            dirmax = []
            for direction in ["X", "Y"]:
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                dirmin.append((1 * f.min()).value)
                dirmax.append((1 * f.max()).value)
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
    # Get command line arguments
    args = parse_args()

    if args.snapshot_number:
        plot_all = False
        snapnr = int(args.snapshot_number)
    else:
        plot_all = True

    domain = args.redshift
    if domain in ("low", "l", "low_redshift", "low redshift", "low-redshift"):
        redshift_domain = "low_redshift"
    elif domain in (
        "medium",
        "m",
        "medium_redshift",
        "medium redshift",
        "medium-redshift",
    ):
        redshift_domain = "medium_redshift"
    elif domain in ("high", "h", "high_redshift", "high redshift", "high-redshift"):
        redshift_domain = "high_redshift"
    else:
        print("Redshift domain not recognised!")
        sys.exit(1)

    snaplist = get_snapshot_list(snapshot_base + f"_{domain}")
    if fancy:
        energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)
    else:
        energy_boundaries = None
        flux_boundaries = None

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
