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
fancy = True  # fancy up the plots a bit
plot_analytical_solutions = True  # overplot analytical solution

# properties for all scatterplots
scatterplot_kwargs = {
    "facecolor": "red",
    "s": 4,
    "alpha": 0.6,
    "linewidth": 0.0,
    "marker": ".",
}

# properties for all analytical solution plots
analytical_solution_kwargs = {"linewidth": 1.0, "ls": "--", "c": "k", "alpha": 0.5}

# -----------------------------------------------------------------------

if plot_analytical_solutions:
    from makeIC import initial_condition

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

    boxsize = meta.boxsize[0]
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

    # get analytical solutions
    if plot_analytical_solutions:

        time = meta.time
        speed = meta.reduced_lightspeed

        advected_positions = data.gas.coordinates[:].copy()
        advected_positions[:, 0] -= speed * time
        nparts = advected_positions.shape[0]
        # add periodicity corrections
        negatives = advected_positions < 0.0
        if negatives.any():
            while advected_positions.min() < 0.0:
                advected_positions[negatives] += boxsize
        overshooters = advected_positions > boxsize
        if overshooters.any():
            while advected_positions.max() > boxsize:
                advected_positions[overshooters] -= boxsize

        analytical_solutions = np.zeros((nparts, ngroups), dtype=np.float64)
        for p in range(part_positions.shape[0]):
            E, F = initial_condition(advected_positions[p])
            for g in range(ngroups):
                analytical_solutions[p, g] = E[g]

    fig = plt.figure(figsize=(5.05 * ngroups, 5.4), dpi=200)
    figname = filename[:-5] + ".png"

    for g in range(ngroups):

        # plot energy density
        new_attribute_str = "radiation_energy" + str(g + 1)
        photon_energy = getattr(data.gas, new_attribute_str)

        volumes = data.gas.masses / data.gas.densities
        photon_energy_density = photon_energy / volumes

        ax = fig.add_subplot(2, ngroups, g + 1)
        s = np.argsort(part_positions)
        if plot_analytical_solutions:
            ax.plot(
                part_positions[s],
                analytical_solutions[s, g],
                **analytical_solution_kwargs,
                label="analytical solution",
            )
        ax.scatter(
            part_positions,
            photon_energy_density,
            **scatterplot_kwargs,
            label="simulation",
        )
        ax.legend()

        ax.set_title("Group {0:2d}".format(g + 1))
        if g == 0:
            ax.set_ylabel(
                "Energy Density [$"
                + photon_energy_density.units.latex_representation()
                + "$]"
            )
        ax.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

        # plot flux X
        new_attribute_str = "radiation_flux" + str(g + 1) + "X"
        photon_flux = getattr(data.gas, new_attribute_str)

        if scheme.startswith("SPH M1closure"):
            photon_flux = photon_flux / volumes

        photon_flux = photon_flux.to("erg/cm**2/s")

        ax2 = fig.add_subplot(2, ngroups, g + 1 + ngroups)
        ax2.scatter(part_positions, photon_flux, **scatterplot_kwargs)

        if g == 0:
            ax2.set_ylabel(
                "Flux X [$" + photon_flux.units.latex_representation() + "$]"
            )
        ax2.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

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


if __name__ == "__main__":

    snaplist = get_snapshot_list(snapshot_base)

    for f in snaplist:
        plot_photons(f)
