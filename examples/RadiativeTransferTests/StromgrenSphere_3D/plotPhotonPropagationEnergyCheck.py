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

import gc
import os
import sys

import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

import stromgren_plotting_tools as spt

# Parameters users should/may tweak

# snapshot basename
snapshot_base = "propagation_test"

# additional anisotropy estimate plot?
plot_anisotropy_estimate = False

# which photon group to use.
# NOTE: array index, not group number (which starts at 1 for GEAR)
group_index = 0

scatterplot_kwargs = {
    "alpha": 0.1,
    "s": 1,
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


def plot_photons():
    """
    Create the actual plot.
    """
    ncols = 2
    fig = plt.figure(figsize=(5 * ncols, 5.5), dpi=200)

    snaplist = get_snapshot_list(snapshot_base)

    timeplot = unyt.unyt_array([])
    Etot_analytic = unyt.unyt_array([])
    Etot_simulation = unyt.unyt_array([])

    for filename in snaplist:
        print("working on", filename)

        # Read in data first
        data = swiftsimio.load(filename)
        meta = data.metadata
        scheme = str(meta.subgrid_scheme["RT Scheme"].decode("utf-8"))
        boxsize = meta.boxsize


        time = meta.time

        L = None

        use_const_emission_rates = False
        if scheme.startswith("GEAR M1closure"):
            luminosity_model = meta.parameters["GEARRT:stellar_luminosity_model"]
            use_const_emission_rates = luminosity_model.decode("utf-8") == "const"
        elif scheme.startswith("SPH M1closure"):
            use_const_emission_rates = bool(
                meta.parameters["SPHM1RT:use_const_emission_rates"]
            )
        else:
            print("Error: Unknown RT scheme " + scheme)
            exit()

        if use_const_emission_rates:
            # read emission rate parameter as string
            if scheme.startswith("GEAR M1closure"):
                const_emission_rates = (
                    spt.trim_paramstr(
                        meta.parameters["GEARRT:const_stellar_luminosities_LSol"].decode(
                            "utf-8"
                        )
                    )
                    * unyt.L_Sun
                )
                L = const_emission_rates[group_index]
            elif scheme.startswith("SPH M1closure"):
                units = data.units
                unit_l_in_cgs = units.length.in_cgs()
                unit_v_in_cgs = (units.length / units.time).in_cgs()
                unit_m_in_cgs = units.mass.in_cgs()
                const_emission_rates = (
                    spt.trim_paramstr(
                        meta.parameters["SPHM1RT:star_emission_rates"].decode("utf-8")
                    )
                    * unit_m_in_cgs
                    * unit_v_in_cgs ** 3
                    / unit_l_in_cgs
                )
                L = const_emission_rates[group_index]
            else:
                print("Error: Unknown RT scheme " + scheme)
                exit()

        # --------------------------
        # Read in and process data
        # --------------------------

        energies = getattr(data.gas.photon_energies, "group" + str(group_index + 1))
        L = L.to(energies.units / time.units)

        xlabel_units_str = time.units.latex_representation()
        energy_units_str = energies.units.latex_representation()
        Etot = L * time

        timeplot = np.append(timeplot,time)  
        Etot_analytic = np.append(Etot_analytic,Etot) 
        Etot_simulation = np.append(Etot_simulation,np.sum(energies)) 

    ax1 = fig.add_subplot(1, ncols, 1)
    ax1.set_title(scheme)
    ax1.set_ylabel("Total Photon Energy [$" + energy_units_str + "$]")

    ax1.plot(timeplot,
        Etot_analytic,
        label="Analytic",
        color="red",
    )
    ax1.plot(timeplot,
        Etot_simulation,
        label="simulations",
        color="blue",
    )

    ax2 = fig.add_subplot(1, ncols, 2)
    ax2.set_title(scheme)
    ax2.set_ylabel("Relative difference")

    ax2.plot(timeplot,
        (Etot_simulation-Etot_analytic)/Etot_analytic,
        color="blue",
    )

    # -------------------------------------------
    # Cosmetics that all axes have in common
    # -------------------------------------------
    for ax in fig.axes:
        ax.set_xlabel("time [$" + xlabel_units_str + "$]")
        ax.legend(fontsize="x-small")

    # Add title
    title = filename.replace("_", r"\_")  # exception handle underscore for latex
    if meta.cosmology is not None:
        title += ", $z$ = {0:.2e}".format(meta.z)
    title += ", $t$ = {0:.2e}".format(meta.time)
    fig.suptitle(title)

    plt.tight_layout()
    figname = filename[:-5]
    figname += "-PhotonPropagationEnergyCheck.png"
    plt.savefig(figname)
    plt.close()
    gc.collect()
    return



if __name__ == "__main__":

    print(
        "REMINDER: Make sure you selected the correct photon group",
        "to plot, which is hardcoded in this script.",
    )
    plot_photons()
