#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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

# ----------------------------------------------------------------------
# plots
#   - radiation energies of particles as function of radius
#   - magnitude of radiation fluxes of particles as function of radius
#   - total energy in radial bins
#   - total vectorial sum of fluxes in radial bins
# and compare with expected propagation speed solution.
# Usage:
#   give snapshot number as cmdline arg to plot
#   single snapshot, otherwise this script plots
#   all snapshots available in the workdir.
#   Make sure to select the photon group to plot that
#   doesn't interact with gas to check the *propagation*
#   correctly.
# ----------------------------------------------------------------------

import sys
import os
import gc
import swiftsimio
import numpy as np
import unyt
from matplotlib import pyplot as plt
import matplotlib as mpl
from scipy import stats

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "propagation_test"

# additional anisotropy estimate plot?
plot_anisotropy_estimate = False

# which photon group to use.
# NOTE: array index, not group number (which starts at 1 for GEAR)
group_index = 0

scatterplot_kwargs = {
    "alpha": 0.3,
    "s": 2,
    "marker": ".",
    "linewidth": 0.0,
    "facecolor": "blue",
}

lineplot_kwargs = {"linewidth": 2}

# -----------------------------------------------------------------------

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True

mpl.rcParams["text.usetex"] = True


def analytical_integrated_energy_solution(L, time, r, rmax):
    """
    Compute analytical solution for the sum of the energy
    in bins for given injection rate <L> at time <time> 
    at bin edges <r> and maximal radius <rmax>
    """

    r_center = 0.5 * (r[:-1] + r[1:])
    r0 = r[0]
    Etot = L * time

    if rmax == 0:
        return r_center, np.zeros(r.shape[0] - 1) * Etot.units

    E = np.zeros(r.shape[0] - 1) * Etot.units
    mask = r_center <= rmax
    E[mask] = Etot / (rmax - r0) * (r[1:] - r[:-1])[mask]

    return r_center, E


def analytical_energy_solution(L, time, r, rmax):
    """
    Compute analytical solution for the energy distribution
    for given injection rate <L> at time <time> at radii <r>
    """

    r_center = 0.5 * (r[:-1] + r[1:])
    r0 = r[0]
    Etot = L * time

    if rmax == 0:
        return r_center, np.zeros(r.shape[0] - 1) * Etot.units

    E_fraction_bin = np.zeros(r.shape[0] - 1) * Etot.units
    mask = r_center <= rmax
    dr = r[1:] - r[:-1]
    E_fraction_bin[mask] = 1.0 / (rmax - r0) * dr[mask]
    bin_surface = dr
    total_weight = Etot / np.sum(E_fraction_bin / bin_surface)
    E = E_fraction_bin / bin_surface * total_weight

    return r_center, E


def analytical_flux_magnitude_solution(L, time, r, rmax, scheme):
    """
    For radiation that doesn't interact with the gas, the
    flux should correspond to the free streaming (optically
    thin) limit. So compute and return that.
    """
    r, E = analytical_energy_solution(L, time, r, rmax)
    if scheme.startswith("GEAR M1closure"):
        F = unyt.c.to(r.units / time.units) * E / r.units ** 3
    elif scheme.startswith("SPH M1closure"):
        F = unyt.c.to(r.units / time.units) * E
    else:
        print("Error: Unknown RT scheme " + scheme)
        exit()

    return r, F


def line(x, a, b):
    return a * x + b


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


def plot_photons(filename, emin, emax, fmin, fmax):
    """
    Create the actual plot.

    filename: file to work with
    emin: list of minimal nonzero energy of all snapshots
    emax: list of maximal energy of all snapshots
    fmin: list of minimal flux magnitude of all snapshots
    fmax: list of maximal flux magnitude of all snapshots
    """

    print("working on", filename)

    # Read in data first
    data = swiftsimio.load(filename)
    meta = data.metadata
    boxsize = meta.boxsize
    edgelen = min(boxsize[0], boxsize[1])

    scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))

    xstar = data.stars.coordinates
    xpart = data.gas.coordinates
    dxp = xpart - xstar
    r = np.sqrt(np.sum(dxp ** 2, axis=1))

    time = meta.time
    r_expect = meta.time * meta.reduced_lightspeed

    use_const_emission_rates = False
    if scheme.startswith("GEAR M1closure"):
        use_const_emission_rates = bool(
            meta.parameters["GEARRT:use_const_emission_rates"]
        )
    elif scheme.startswith("SPH M1closure"):
        use_const_emission_rates = bool(
            meta.parameters["SPHM1RT:use_const_emission_rates"]
        )
    else:
        print("RT scheme not identified. Exit.")
        exit()
    L = None

    if use_const_emission_rates:
        # read emission rate parameter as string
        if scheme.startswith("GEAR M1closure"):
            emissionstr = meta.parameters["GEARRT:star_emission_rates_LSol"].decode(
                "utf-8"
            )
            # clean string up
            if emissionstr.startswith("["):
                emissionstr = emissionstr[1:]
            if emissionstr.endswith("]"):
                emissionstr = emissionstr[:-1]

            # transform string values to floats with unyts
            emissions = emissionstr.split(",")
            emlist = []
            for er in emissions:
                emlist.append(float(er))
            const_emission_rates = unyt.unyt_array(emlist, unyt.L_Sun)
            L = const_emission_rates[group_index]
        elif scheme.startswith("SPH M1closure"):
            unit_l_in_cgs = meta.units.length
            unit_v_in_cgs = meta.units.length / meta.units.time
            unit_m_in_cgs = meta.units.mass
            emissionstr = meta.parameters["SPHM1RT:star_emission_rates"].decode("utf-8")
            # clean string up
            if emissionstr.startswith("["):
                emissionstr = emissionstr[1:]
            if emissionstr.endswith("]"):
                emissionstr = emissionstr[:-1]

            # transform string values to floats with unyts
            emissions = emissionstr.split(",")
            emlist = []
            for er in emissions:
                emlist.append(
                    float(er) * unit_m_in_cgs * unit_v_in_cgs ** 3 / unit_l_in_cgs
                )
            const_emission_rates = unyt.unyt_array(emlist, "erg/s")
            L = const_emission_rates[group_index]
        else:
            print("Error: Unknown RT scheme " + scheme)
            exit()

    if plot_anisotropy_estimate:
        ncols = 4
    else:
        ncols = 3
    fig = plt.figure(figsize=(5 * ncols, 5.5), dpi=200)

    nbins = 100
    r_bin_edges = np.linspace(0.5 * edgelen * 1e-2, 0.507 * edgelen, nbins + 1)
    r_bin_centres = 0.5 * (r_bin_edges[1:] + r_bin_edges[:-1])
    r_analytical_bin_edges = np.linspace(
        0.5 * edgelen * 1e-2, 0.507 * edgelen, nbins + 1
    )

    # --------------------------
    # Read in and process data
    # --------------------------

    energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
    Fx = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "X")
    Fy = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "Y")

    fmag = np.sqrt(Fx ** 2 + Fy ** 2)
    particle_count, _ = np.histogram(
        r,
        bins=r_analytical_bin_edges,
        range=(r_analytical_bin_edges[0], r_analytical_bin_edges[-1]),
    )
    L = L.to(energies.units / time.units)

    xlabel_units_str = boxsize.units.latex_representation()
    energy_units_str = energies.units.latex_representation()
    flux_units_str = Fx.units.latex_representation()

    # ------------------------
    # Plot photon energies
    # ------------------------
    ax1 = fig.add_subplot(1, ncols, 1)
    ax1.set_title("Particle Radiation Energies")
    ax1.set_ylabel("Photon Energy [$" + energy_units_str + "$]")

    # don't expect more than float precision
    emin_to_use = max(emin, 1e-5 * emax)

    if use_const_emission_rates:
        # plot entire expected solution
        rA, EA = analytical_energy_solution(L, time, r_analytical_bin_edges, r_expect)

        mask = particle_count > 0
        if mask.any():
            EA = EA[mask].to(energies.units)
            rA = rA[mask]
            pcount = particle_count[mask]

            # the particle bin counts will introduce noise.
            # So use a linear fit for the plot. I assume here
            # that the particle number per bin increases
            # proprtional to r, which should roughly be the
            # case for the underlying glass particle distribution.
            lin_res = stats.linregress(rA, pcount)

            ax1.plot(
                rA,
                EA / line(rA.v, lin_res.slope, lin_res.intercept),
                **lineplot_kwargs,
                linestyle="--",
                c="red",
                label="Analytical Solution",
            )

    else:
        # just plot where photon front should be
        ax1.plot(
            [r_expect, r_expect],
            [emin_to_use, emax * 1.1],
            label="expected photon front",
            color="red",
        )

    ax1.scatter(r, energies, **scatterplot_kwargs)
    energies_binned, _, _ = stats.binned_statistic(
        r,
        energies,
        statistic="mean",
        bins=r_bin_edges,
        range=(r_bin_edges[0], r_bin_edges[-1]),
    )
    ax1.plot(
        r_bin_centres, energies_binned, **lineplot_kwargs, label="Mean Radiation Energy"
    )
    ax1.set_ylim(emin_to_use, emax * 1.1)

    # ------------------------------
    # Plot binned photon energies
    # ------------------------------
    ax2 = fig.add_subplot(1, ncols, 2)
    ax2.set_title("Total Radiation Energy in radial bins")
    ax2.set_ylabel("Total Photon Energy [$" + energy_units_str + "$]")

    energies_summed_bin, _, _ = stats.binned_statistic(
        r,
        energies,
        statistic="sum",
        bins=r_bin_edges,
        range=(r_bin_edges[0], r_bin_edges[-1]),
    )
    ax2.plot(
        r_bin_centres,
        energies_summed_bin,
        **lineplot_kwargs,
        label="Total Energy in Bin",
    )
    current_ylims = ax2.get_ylim()
    ax2.set_ylim(emin_to_use, current_ylims[1])

    if use_const_emission_rates:
        # plot entire expected solution
        # Note: you need to use the same bins as for the actual results
        rA, EA = analytical_integrated_energy_solution(L, time, r_bin_edges, r_expect)

        ax2.plot(
            rA,
            EA.to(energies.units),
            **lineplot_kwargs,
            linestyle="--",
            c="red",
            label="Analytical Solution",
        )
    else:
        # just plot where photon front should be
        ax2.plot(
            [r_expect, r_expect],
            ax2.get_ylim(r),
            label="Expected Photon Front",
            color="red",
        )

    # ------------------------------
    # Plot photon fluxes
    # ------------------------------
    ax3 = fig.add_subplot(1, ncols, 3)
    ax3.set_title("Particle Radiation Flux Magnitudes")
    ax3.set_ylabel("Photon Flux Magnitude [$" + flux_units_str + "$]")

    fmin_to_use = max(fmin, 1e-5 * fmax)
    ax3.set_ylim(fmin_to_use, fmax * 1.1)

    ax3.scatter(r, fmag, **scatterplot_kwargs)

    fmag_mean_bin, _, _ = stats.binned_statistic(
        r,
        fmag,
        statistic="mean",
        bins=r_bin_edges,
        range=(r_bin_edges[0], r_bin_edges[-1]),
    )
    ax3.plot(
        r_bin_centres,
        fmag_mean_bin,
        **lineplot_kwargs,
        label="Mean Radiation Flux of particles",
    )

    if use_const_emission_rates:
        # plot entire expected solution
        rA, FA = analytical_flux_magnitude_solution(
            L, time, r_analytical_bin_edges, r_expect, scheme
        )

        mask = particle_count > 0
        if mask.any():
            FA = FA[mask].to(Fx.units)
            rA = rA[mask]
            pcount = particle_count[mask]

            # the particle bin counts will introduce noise.
            # So use a linear fit for the plot. I assume here
            # that the particle number per bin increases
            # proprtional to r, which should roughly be the
            # case for the underlying glass particle distribution.
            lin_res = stats.linregress(rA, pcount)

            ax3.plot(
                rA,
                FA / line(rA.v, lin_res.slope, lin_res.intercept),
                **lineplot_kwargs,
                linestyle="--",
                c="red",
                label="analytical solution",
            )

    else:
        # just plot where photon front should be
        ax1.plot(
            [r_expect, r_expect],
            [emin_to_use, emax * 1.1],
            label="expected photon front",
            color="red",
        )

    # ------------------------------
    # Plot photon flux sum
    # ------------------------------

    if plot_anisotropy_estimate:

        ax4 = fig.add_subplot(1, ncols, 4)
        ax4.set_title("Vectorial Sum of Radiation Flux in radial bins")
        ax4.set_ylabel("[1]")

        fmag_sum_bin, _, _ = stats.binned_statistic(
            r,
            fmag,
            statistic="sum",
            bins=r_bin_edges,
            range=(r_bin_edges[0], r_bin_edges[-1]),
        )
        mask_sum = fmag_sum_bin > 0
        fmag_max_bin, _, _ = stats.binned_statistic(
            r,
            fmag,
            statistic="max",
            bins=r_bin_edges,
            range=(r_bin_edges[0], r_bin_edges[-1]),
        )
        mask_max = fmag_max_bin > 0
        Fx_sum_bin, _, _ = stats.binned_statistic(
            r,
            Fx,
            statistic="sum",
            bins=r_bin_edges,
            range=(r_bin_edges[0], r_bin_edges[-1]),
        )
        Fy_sum_bin, _, _ = stats.binned_statistic(
            r,
            Fy,
            statistic="sum",
            bins=r_bin_edges,
            range=(r_bin_edges[0], r_bin_edges[-1]),
        )
        F_sum_bin = np.sqrt(Fx_sum_bin ** 2 + Fy_sum_bin ** 2)

        ax4.plot(
            r_bin_centres[mask_sum],
            F_sum_bin[mask_sum] / fmag_sum_bin[mask_sum],
            **lineplot_kwargs,
            label="$\left| \sum_{i \in \mathrm{particles \ in \ bin}} \mathbf{F}_i \\right| $ / $\sum_{i \in \mathrm{particles \ in \ bin}} \left| \mathbf{F}_{i} \\right| $",
        )
        ax4.plot(
            r_bin_centres[mask_max],
            F_sum_bin[mask_max] / fmag_max_bin[mask_max],
            **lineplot_kwargs,
            linestyle="--",
            label="$\left| \sum_{i \in \mathrm{particles \ in \ bin}} \mathbf{F}_i \\right| $ / $\max_{i \in \mathrm{particles \ in \ bin}} \left| \mathbf{F}_{i} \\right| $",
        )

    # -------------------------------------------
    # Cosmetics that all axes have in common
    # -------------------------------------------
    for ax in fig.axes:
        ax.set_xlabel("r [$" + xlabel_units_str + "$]")
        ax.set_yscale("log")
        ax.set_xlim(0.0, 0.501 * edgelen)
        ax.legend(fontsize="x-small")

    # Add title
    title = filename.replace("_", "\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time)
    fig.suptitle(title)

    plt.tight_layout()
    figname = filename[:-5]
    figname += "-PhotonPropagation.png"
    plt.savefig(figname)
    plt.close()
    gc.collect()

    return


def get_plot_boundaries(filenames):
    """
    Get minimal and maximal nonzero photon energy values
    """

    data = swiftsimio.load(filenames[0])
    energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
    emaxguess = energies.max()

    emin = emaxguess
    emax = 0.0
    fmagmin = 1e30
    fmagmax = -10.0

    for f in filenames:
        data = swiftsimio.load(f)

        energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
        mask = energies > 0.0

        if mask.any():

            nonzero_energies = energies[mask]
            this_emin = nonzero_energies.min()
            emin = min(this_emin, emin)

            this_emax = energies.max()
            emax = max(emax, this_emax)

        fx = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "X")
        fy = getattr(data.gas.photon_fluxes, "Group" + str(group_index + 1) + "Y")
        fmag = np.sqrt(fx ** 2 + fy ** 2)

        fmagmin = min(fmagmin, fmag.min())
        fmagmax = max(fmagmax, fmag.max())

    return emin, emax, fmagmin, fmagmax


if __name__ == "__main__":

    print(
        "REMINDER: Make sure you selected the correct photon group",
        "to plot, which is hardcoded in this script.",
    )
    snaplist = get_snapshot_list(snapshot_base)
    emin, emax, fmagmin, fmagmax = get_plot_boundaries(snaplist)
    for f in snaplist:
        plot_photons(f, emin, emax, fmagmin, fmagmax)
