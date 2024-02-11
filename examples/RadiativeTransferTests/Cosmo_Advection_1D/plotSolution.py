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
import unyt
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as cf

# Parameters users should/may tweak
plot_all_data = True  # plot all groups and all photon quantities
snapshot_base = "output"  # snapshot basename
fancy = True  # fancy up the plots a bit
plot_analytical_solutions = True  # overplot analytical solution
plot_physical_quantities = True  # Plot physiical or comoving energy densities

time_first = 0
time_prev = 0
H_first = 0

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


def plot_param_over_time(snapshot_list, param="energy density"):
    print(f"Now plotting {param} over time")
    # Grab number of photon groups
    data = swiftsimio.load(snapshot_list[0])
    meta = data.metadata
    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"][0])

    # Number of rows and columns
    nrows = 1 + int(plot_physical_quantities)
    ncols = ngroups
    # Create figure and axes
    fig, axs = plt.subplots(nrows, ncols, figsize=(5.04 * ncols, 5.4 * nrows), dpi=200)

    # Iterate over all photon groups
    for n in range(ngroups):
        # Arrays to keep track of plot_param and scale factor
        plot_param = [[], []]
        scale_factor = []
        analytic_exponent = [0, 0]

        # Functions to convert between scale factor and redshift
        a2z = lambda a: 1 / a - 1
        z2a = lambda z: 1 / (z + 1)

        for file in snapshot_list:
            data = swiftsimio.load(file)
            meta = data.metadata

            # Read comoving energy density and energy
            energy_density = data.gas.photon_energy_densities[:, n]
            energy = getattr(data.gas.photon_energies, f"group{n+1}")

            if plot_physical_quantities:
                # The SWIFT cosmology module assumes 3-dimensional lengths and volumes,
                # so multiply by a**2 to get the correct relations
                physical_energy_density = (
                    energy_density.to_physical() * meta.scale_factor ** 2
                )
                physical_energy = energy.to_physical()

                if param == "energy density":
                    plot_param[1].append(
                        1
                        * physical_energy_density.sum()
                        / physical_energy_density.shape[0]
                    )
                    analytic_exponent[1] = -1.0
                elif param == "total energy":
                    plot_param[1].append(1 * physical_energy.sum())
                    analytic_exponent[1] = 0.0

            if param == "energy density":
                plot_param[0].append(1 * energy_density.sum() / energy_density.shape[0])
                analytic_exponent[0] = 0.0
            elif param == "total energy":
                plot_param[0].append(1 * energy.sum())
                analytic_exponent[0] = 0.0

            scale_factor.append(meta.scale_factor)

        if param == "energy density":
            titles = [
                "Comoving energy density",
                "Physical energy density $\\times a^2$",
            ]
            ylabel = "Average energy density"
            figname = "output_energy_density_over_time.png"
        elif param == "total energy":
            titles = ["Comoving total energy", "Physical total energy"]
            ylabel = "Total energy"
            figname = "output_total_energy_over_time.png"

        # Analytic scale factor
        analytic_scale_factor = np.linspace(min(scale_factor), max(scale_factor), 1000)

        for i in range(nrows):
            ax = axs[i, n]
            ax.scatter(scale_factor, plot_param[i], label="Simulation")

            # Analytic scale factor relation
            analytic = analytic_scale_factor ** analytic_exponent[i]

            # Scale solution to correct offset
            analytic = analytic / analytic[0] * plot_param[i][0]
            ax.plot(
                analytic_scale_factor,
                analytic,
                c="r",
                label=f"Analytic solution $\propto a^{{{analytic_exponent[i]}}}$",
            )

            ax.legend()
            ax.set_title(titles[i] + f" group {n+1}")

            ax.set_xlabel("Scale factor")
            secax = ax.secondary_xaxis("top", functions=(a2z, z2a))
            secax.set_xlabel("Redshift")

            ax.yaxis.get_offset_text().set_position((-0.05, 1))

            if analytic_exponent[i] == 0.0:
                ax.set_ylim(plot_param[i][0] * 0.95, plot_param[i][0] * 1.05)
            if n == 0:
                units = plot_param[i][0].units.latex_representation()
                ax.set_ylabel(f"{ylabel} [${units}$]")
    plt.tight_layout()
    plt.savefig(figname)
    plt.close()


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
        # Grab unit system used in snapshot
        snapshot_units = meta.units
        snapshot_unit_energy = (
            snapshot_units.mass * snapshot_units.length ** 2 / snapshot_units.time ** 2
        )

        # Grab unit system used in IC
        IC_units = swiftsimio.units.cosmo_units
        plot_unit_energy = (
            IC_units["mass"] * IC_units["length"] ** 2 / IC_units["time"] ** 2
        )

        # Factor to convert between the unit systems
        conversion_factor = plot_unit_energy / snapshot_unit_energy

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
                analytical_solutions[p, g] = E[g] * conversion_factor

    # Plot plot plot!
    if plot_all_data:
        fig = plt.figure(figsize=(5.05 * ngroups, 5.4), dpi=200)
        figname = filename[:-5] + "-all-quantities.png"

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

        ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])
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

    snaplist = get_snapshot_list(snapshot_base)
    if fancy:
        energy_boundaries, flux_boundaries = get_minmax_vals(snaplist)
    else:
        energy_boundaries, flux_boundaries = (None, None)

    for f in snaplist:
        plot_photons(
            f, energy_boundaries=energy_boundaries, flux_boundaries=flux_boundaries
        )
    # Only plot over time if more than 2 snapshots
    if len(snaplist) > 2:
        for param in ["energy density", "total energy"]:
            plot_param_over_time(snaplist, param)
