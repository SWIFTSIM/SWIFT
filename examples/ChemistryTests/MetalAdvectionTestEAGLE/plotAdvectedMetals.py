#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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

import math

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import unyt

import swiftsimio
from swiftsimio.visualisation import project_gas
from pathlib import Path

from makeIC import VELOCITY


def plot_single(ax_mass, ax_fraction, name, title, data, mass_map, kwargs_inner):
    mass_weighted_map = project_gas(data, project=name, **kwargs_inner["projection"])
    ax_mass.imshow(mass_weighted_map.in_cgs().value.T, **kwargs_inner["imshow_mass"])
    ax_mass.set_title(title)
    ax_fraction.imshow(
        (mass_weighted_map / mass_map).T, **kwargs_inner["imshow_fraction"]
    )
    ax_fraction.set_title(title)
    ax_mass.axis("off")
    ax_fraction.axis("off")


def plot_all(fname, savename):
    data = swiftsimio.load(fname)

    # Shift coordinates to starting position
    velocity = 0.5 * math.sqrt(2) * VELOCITY * np.array([1, 1, 0]) * unyt.cm / unyt.s
    data.gas.coordinates -= data.metadata.time * velocity

    # Add mass weighted element mass fractions to gas dataset
    masses = data.gas.masses
    element_mass_fractions = data.gas.element_mass_fractions
    data.gas.element_mass_H = element_mass_fractions.hydrogen * masses
    data.gas.element_mass_He = element_mass_fractions.helium * masses
    data.gas.element_mass_C = element_mass_fractions.carbon * masses
    data.gas.element_mass_N = element_mass_fractions.nitrogen * masses
    data.gas.element_mass_O = element_mass_fractions.oxygen * masses
    data.gas.element_mass_Ne = element_mass_fractions.neon * masses
    data.gas.element_mass_Mg = element_mass_fractions.magnesium * masses
    data.gas.element_mass_Si = element_mass_fractions.silicon * masses
    data.gas.element_mass_Fe = element_mass_fractions.iron * masses
    data.gas.total_metal_mass = data.gas.metal_mass_fractions * masses

    # Create necessary figures and axes
    fig = plt.figure(layout="constrained", figsize=(16, 8))
    fig.suptitle(
        f"Profiles shifted to starting position after t={data.metadata.time:.2f}",
        fontsize=14,
    )
    fig_ratios, fig_masses = fig.subfigures(2, 1)
    fig_ratios.suptitle("Mass ratio of elements")
    fig_masses.suptitle("Surface density in elements")
    axes_ratios = fig_ratios.subplots(2, 5, sharex=True, sharey=True)
    axes_masses = fig_masses.subplots(2, 5, sharex=True, sharey=True)

    # parameters for swiftsimio projections
    projection_kwargs = {
        "region": np.array([0, 2, 0, 1, 0, 1]) * unyt.cm,
        "resolution": 500,
        "parallel": True,
    }
    # Parameters for imshow
    norm_ratios = Normalize(0, 0.95)
    norm_masses = Normalize(0, 0.9)
    imshow_fraction_kwargs = dict(norm=norm_ratios, cmap="rainbow")
    imshow_mass_kwargs = dict(norm=norm_masses, cmap="turbo")

    # Plot the quantities
    mass_map = project_gas(data, project="masses", **projection_kwargs)

    # Parameters for plotting:
    plotting_kwargs = dict(
        data=data,
        mass_map=mass_map,
        kwargs_inner=dict(
            projection=projection_kwargs,
            imshow_mass=imshow_mass_kwargs,
            imshow_fraction=imshow_fraction_kwargs,
        ),
    )

    plot_single(
        axes_masses[0][0],
        axes_ratios[0][0],
        "total_metal_mass",
        "All metals",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[0][1],
        axes_ratios[0][1],
        "element_mass_H",
        "Hydrogen",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[0][2],
        axes_ratios[0][2],
        "element_mass_He",
        "Helium",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[0][3],
        axes_ratios[0][3],
        "element_mass_C",
        "Carbon",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[0][4],
        axes_ratios[0][4],
        "element_mass_N",
        "Nitrogen",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[1][0],
        axes_ratios[1][0],
        "element_mass_O",
        "Oxygen",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[1][1],
        axes_ratios[1][1],
        "element_mass_Ne",
        "Neon",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[1][2],
        axes_ratios[1][2],
        "element_mass_Mg",
        "Magnesium",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[1][3],
        axes_ratios[1][3],
        "element_mass_Si",
        "Silicon",
        **plotting_kwargs,
    )
    plot_single(
        axes_masses[1][4],
        axes_ratios[1][4],
        "element_mass_Fe",
        "Iron",
        **plotting_kwargs,
    )

    # Add Colorbars
    cb_masses = fig_masses.colorbar(
        ScalarMappable(**imshow_mass_kwargs), shrink=0.75, pad=0.01, ax=axes_masses
    )
    cb_ratios = fig_ratios.colorbar(
        ScalarMappable(**imshow_fraction_kwargs), shrink=0.75, pad=0.01, ax=axes_ratios
    )
    cb_masses.ax.set_ylabel("Surface density (g/cm^2)", rotation=270, labelpad=15)
    cb_ratios.ax.set_ylabel("Mass ratio", rotation=270, labelpad=15)

    # Save output
    fig.savefig(savename, dpi=300)


if __name__ == "__main__":
    cwd = Path(__file__).parent
    print("Plotting metals for output_0000.hdf5...")
    plot_all(cwd / "output_0000.hdf5", savename=cwd / "metal_advection_0000.png")
    print("Plotting metals for output_0001.hdf5...")
    plot_all(cwd / "output_0001.hdf5", savename=cwd / "metal_advection_0001.png")
    print("Done!")
