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


# ----------------------------------------------
# plot photon data assuming a 1D problem
# give snapshot number as cmdline arg to plot
# single snapshot, otherwise this script plots
# all snapshots available in the workdir
# ----------------------------------------------

import os
import sys
import argparse

import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt

# Parameters users should/may tweak
plot_all_data = True  # plot all groups and all photon quantities
snapshot_base = "output"  # snapshot basename
fancy = True  # fancy up the plots a bit
plot_analytical_solutions = True  # overplot analytical solution
plot_physical_quantities = True  # Plot physiical or comoving energy densities

time_first = 0

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
plot_all = False  # plot all snapshots


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-n", "--snapshot-number", help="Number of snaphot to plot", type=int
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


def plot_photons(filename, energy_boundaries=None, flux_boundaries=None):
    """
    Create the actual plot.

    filename: file to work with
    energy_boundaries:  list of [E_min, E_max] for each photon group.
                        If none, limits are set automatically.
    flux_boundaries:    list of [F_min, F_max] for each photon group.
                        If none, limits are set automatically.
    """
    global time_first
    print("working on", filename)

    # Read in data firt
    data = swiftsimio.load(filename)
    meta = data.metadata
    cosmology = meta.cosmology_raw

    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    boxsize = meta.boxsize[0]

    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])
    # Currently, SPHM1RT only works for frequency group = 4 in the code
    # However, we only plot 3 frequency groups here, so
    # we set ngroups = 3:
    if scheme.startswith("SPH M1closure"):
        ngroups = 3

    for g in range(ngroups):
        # workaround to access named columns data with swiftsimio visualisaiton
        new_attribute_str = "radiation_energy" + str(g + 1)
        en = getattr(data.gas.photon_energies, "group" + str(g + 1))
        setattr(data.gas, new_attribute_str, en)

        if plot_all_data:
            # prepare also the fluxes
            for direction in ["X"]:
                new_attribute_str = "radiation_flux" + str(g + 1) + direction
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                setattr(data.gas, new_attribute_str, f)

    part_positions = data.gas.coordinates[:, 0].copy()

    # get analytical solutions
    if plot_analytical_solutions:
        # Grab unit system used in IC
        IC_units = swiftsimio.units.cosmo_units

        # Grab unit system used in snapshot
        snapshot_units = meta.units

        snapshot_unit_energy = (
            snapshot_units.mass * snapshot_units.length ** 2 / snapshot_units.time ** 2
        )

        time = meta.time.copy()

        # Take care of simulation time not starting at 0
        if time_first == 0:
            time_first = time
            time = 0 * snapshot_units.time
        else:
            time -= time_first

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
            E, F = initial_condition(advected_positions[p], IC_units)
            for g in range(ngroups):
                analytical_solutions[p, g] = E[g] / snapshot_unit_energy

    # Plot plot plot!
    if plot_all_data:

        fig = plt.figure(figsize=(5.05 * ngroups, 5.4), dpi=200)
        figname = filename[:-5] + f"-all-quantities.png"

        for g in range(ngroups):

            # plot energy
            new_attribute_str = "radiation_energy" + str(g + 1)
            photon_energy = getattr(data.gas, new_attribute_str)

            # physical_photon_energy = photon_energy.to_physical()
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
                    fixed_min = energy_boundaries[g][0] - 0.1 * abs(
                        energy_boundaries[g][1]
                    )
                    fixed_max = energy_boundaries[g][1] * 1.1
                else:
                    fixed_min = energy_boundaries[g][0] * 1.1
                    fixed_max = energy_boundaries[g][1] + 0.1 * abs(
                        energy_boundaries[g][0]
                    )
                ax.set_ylim(fixed_min, fixed_max)

            # plot flux X
            new_attribute_str = "radiation_flux" + str(g + 1) + "X"
            photon_flux = getattr(data.gas, new_attribute_str)
            physical_photon_flux = photon_flux.to_physical()

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
                ax.set_ylabel(
                    "Flux X [$" + photon_flux.units.latex_representation() + "$]"
                )
            ax.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

            if flux_boundaries is not None:

                if abs(flux_boundaries[g][1]) > abs(flux_boundaries[g][0]):
                    fixed_min = flux_boundaries[g][0] - 0.1 * abs(flux_boundaries[g][1])
                    fixed_max = flux_boundaries[g][1] * 1.1
                else:
                    fixed_min = flux_boundaries[g][0] * 1.1
                    fixed_max = flux_boundaries[g][1] + 0.1 * abs(flux_boundaries[g][0])

                ax.set_ylim(fixed_min, fixed_max)

    else:  # plot just energies

        fig = plt.figure(figsize=(5 * ngroups, 5), dpi=200)
        figname = filename[:-5] + ".png"

        for g in range(ngroups):

            ax = fig.add_subplot(1, ngroups, g + 1)

            new_attribute_str = "radiation_energy" + str(g + 1)
            photon_energy = getattr(data.gas, new_attribute_str).to_physical()

            s = np.argsort(part_positions)
            if plot_analytical_solutions:
                ax.plot(
                    part_positions[s],
                    analytical_solutions[s, g],
                    **analytical_solution_kwargs,
                    label="analytical solution",
                )
            ax.scatter(
                part_positions, photon_energy, **scatterplot_kwargs, label="simulation"
            )
            ax.set_title("Group {0:2d}".format(g + 1))
            if g == 0:
                ax.set_ylabel(
                    "Energies [$" + photon_energy.units.latex_representation() + "$]"
                )
            ax.set_xlabel("x [$" + part_positions.units.latex_representation() + "$]")

            if energy_boundaries is not None:

                if abs(energy_boundaries[g][1]) > abs(energy_boundaries[g][0]):
                    fixed_min = energy_boundaries[g][0] - 0.1 * abs(
                        energy_boundaries[g][1]
                    )
                    fixed_max = energy_boundaries[g][1] * 1.1
                else:
                    fixed_min = energy_boundaries[g][0] * 1.1
                    fixed_max = energy_boundaries[g][1] + 0.1 * abs(
                        energy_boundaries[g][0]
                    )
                ax.set_ylim(fixed_min, fixed_max)

    # add title
    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.3e}".format(1 * meta.time)
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

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])
        emin_group = []
        emax_group = []
        fluxmin_group = []
        fluxmax_group = []

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            emin_group.append((1 * en.min()).value)
            emax_group.append((1 * en.max()).value)

            for direction in ["X"]:
                #  for direction in ["X", "Y", "Z"]:
                f = getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + direction)
                fluxmin_group.append(
                    (1 * f.min()).to(unyt.erg / unyt.cm ** 2 / unyt.s).value
                )
                fluxmax_group.append(
                    (1 * f.max()).to(unyt.erg / unyt.cm ** 2 / unyt.s).value
                )

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
    snaplist = get_snapshot_list(snapshot_base + f"_{domain}")

    if fancy:
        energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)
    else:
        energy_boundaries, flux_boundaries = (None, None)

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
