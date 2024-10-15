#!/usr/bin/env python3

###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Stan Verhoeve (s06verhoeve@gmail.com)
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

import os
import sys
import argparse

import matplotlib as mpl
import numpy as np
import swiftsimio
import unyt
from matplotlib import pyplot as plt

# Parameters users should/may tweak
snapshot_base = "output"  # Snapshot basename
plot_physical_quantities = True  # Plot physical or comoving quantities

time_first = 0  # Time of first snapshot


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-z", "--redshift", help="Redshift domain to plot advection for", default="high"
    )

    args = parser.parse_args()
    return args


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted
    and return their names as list
    """

    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            snaplist.append(f)

    snaplist = sorted(snaplist)
    if len(snaplist) == 0:
        print(f"No snapshots with base {snapshot_basename} found!")
        sys.exit(1)
    return snaplist


def plot_param_over_time(snapshot_list, param="energy density", redshift_domain="high"):
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

            # Read comoving variables
            energy = getattr(data.gas.photon_energies, f"group{n+1}")
            mass = data.gas.masses
            rho = data.gas.densities
            volume = mass / rho

            energy_density = energy / volume

            if plot_physical_quantities:
                # The SWIFT cosmology module assumes 3-dimensional lengths and volumes,
                # so multiply by a to get the correct relations
                physical_energy_density = (
                    energy_density.to_physical() * meta.scale_factor
                )
                physical_energy = energy.to_physical()

                if param == "energy density":
                    plot_param[1].append(
                        1
                        * physical_energy_density.sum()
                        / physical_energy_density.shape[0]
                    )
                    analytic_exponent[1] = -2.0
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
            titles = ["Comoving energy density", "Physical energy density $\\times a$"]
            ylabel = "Average energy density"
            figname = "output_energy_density_over_time"
        elif param == "total energy":
            titles = ["Comoving total energy", "Physical total energy"]
            ylabel = "Total energy"
            figname = "output_total_energy_over_time"

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
                ylabel_scale = ""
            else:
                ylabel_scale = "$\\times a$"
            if n == 0:
                units = plot_param[i][0].units.latex_representation()
                ax.set_ylabel(f"{ylabel} [${units}$] {ylabel_scale}")
    plt.tight_layout()
    plt.savefig(f"{figname}-{redshift_domain}.png")
    plt.close()


if __name__ in ("__main__"):
    # Get command line args
    args = parse_args()
    domain = args.redshift.lower()
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

    snaplist = get_snapshot_list(snapshot_base + f"_{redshift_domain}")

    for param in ["energy density", "total energy"]:
        plot_param_over_time(snaplist, param, redshift_domain)
