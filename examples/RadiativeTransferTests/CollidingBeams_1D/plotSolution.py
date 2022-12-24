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


# ----------------------------------------------
# plot photon data assuming a 1D problem
# give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir
# ----------------------------------------------

import os
import sys

import matplotlib as mpl
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt

# Parameters users should/may tweak
snapshot_base = "output"  # snapshot basename

# properties for all scatterplots
scatterplot_kwargs = {
    "facecolor": "red",
    "s": 4,
    "alpha": 0.6,
    "linewidth": 0.0,
    "marker": ".",
}

# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False  # plot all snapshots
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True


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

    # Read in data firt
    data = swiftsimio.load(filename)
    meta = data.metadata
    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])

    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        new_attribute_str = "radiation_energy" + str(g + 1)
        en = getattr(data.gas.photon_energies, "group" + str(g + 1))
        setattr(data.gas, new_attribute_str, en)

        # prepare also the fluxes
        for direction in ["X"]:
            new_attribute_str = "radiation_flux" + str(g + 1) + direction
            f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
            setattr(data.gas, new_attribute_str, f)

    part_positions = data.gas.coordinates[:, 0].copy()

    # Plot plot plot!
    fig = plt.figure(figsize=(5.0 * ngroups, 5.4), dpi=200)
    figname = filename[:-5] + "-all-quantities.png"

    for g in range(ngroups):

        # plot energy
        new_attribute_str = "radiation_energy" + str(g + 1)
        photon_energy = getattr(data.gas, new_attribute_str)

        ax = fig.add_subplot(2, ngroups, g + 1)

        ax.scatter(
            part_positions, photon_energy, **scatterplot_kwargs, label="simulation"
        )
        ax.legend()

        ax.set_title("Group {0:2d}".format(g + 1))
        if g == 0:
            ax.set_ylabel(
                "Energies [$" + photon_energy.units.latex_representation() + "$]"
            )
        ax.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

        if energy_boundaries is not None:

            if abs(energy_boundaries[g][1]) > abs(energy_boundaries[g][0]):
                fixed_min = energy_boundaries[g][0] - 0.1 * abs(energy_boundaries[g][1])
                fixed_max = energy_boundaries[g][1] * 1.1
            else:
                fixed_min = energy_boundaries[g][0] * 1.1
                fixed_max = energy_boundaries[g][1] + 0.1 * abs(energy_boundaries[g][0])
            ax.set_ylim(fixed_min, fixed_max)

        # plot flux X
        new_attribute_str = "radiation_flux" + str(g + 1) + "X"
        photon_flux = getattr(data.gas, new_attribute_str)
        if scheme.startswith("GEAR M1closure"):
            photon_flux = photon_flux.to("erg/cm**2/s")
        elif scheme.startswith("SPH M1closure"):
            photon_flux = photon_flux.to("erg*cm/s")
        else:
            print("Error: Unknown RT scheme " + scheme)
            exit()

        ax = fig.add_subplot(2, ngroups, g + 1 + ngroups)
        ax.scatter(part_positions, photon_flux, **scatterplot_kwargs)

        if g == 0:
            ax.set_ylabel("Flux X [$" + photon_flux.units.latex_representation() + "$]")
        ax.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

        if flux_boundaries is not None:

            if abs(flux_boundaries[g][1]) > abs(flux_boundaries[g][0]):
                fixed_min = flux_boundaries[g][0] - 0.1 * abs(flux_boundaries[g][1])
                fixed_max = flux_boundaries[g][1] * 1.1
            else:
                fixed_min = flux_boundaries[g][0] * 1.1
                fixed_max = flux_boundaries[g][1] + 0.1 * abs(flux_boundaries[g][0])

            ax.set_ylim(fixed_min, fixed_max)

    # add title
    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time)
    fig.suptitle(title)

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

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

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
        emin_group = []
        emax_group = []
        fluxmin_group = []
        fluxmax_group = []

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            emin_group.append(en.min())
            emax_group.append(en.max())

            for direction in ["X"]:
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                fluxmin_group.append(f.min())
                fluxmax_group.append(f.max())

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
    energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
