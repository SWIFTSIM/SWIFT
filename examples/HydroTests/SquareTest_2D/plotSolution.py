###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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
"""
Plots the solution of the square test in a smoothed way using SWIFTsimIO's 
smoothed plotting.
"""

import sys
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid

snap = int(sys.argv[1])

sim = load(f"square_{snap:04d}.hdf5")
resolution = 512

# First create a grid that gets the particle density so we can divide it out later
unweighted_grid = project_gas_pixel_grid(sim, 512, None)

# Set up plotting stuff
try:
    plt.style.use("mnras_durham")
except:
    rcParams = {
        "font.serif": ["STIX", "Times New Roman", "Times"],
        "font.family": ["serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
    }
    plt.rcParams.update(rcParams)


def get_data_dump(metadata):
    """
    Gets a big data dump from the SWIFT metadata
    """

    try:
        viscosity = metadata.viscosity_info
    except:
        viscosity = "No info"

    try:
        diffusion = metadata.diffusion_info
    except:
        diffusion = "No info"

    output = (
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    return output


# Now we can do the plotting.
fig, ax = plt.subplots(2, 3, figsize=(6.974, 6.974 * (2.0 / 3.0)))
ax = ax.flatten()

# These are stored in priority order
plot = dict(
    internal_energies="Internal Energy ($u$)",
    densities=r"Density ($\rho$)",
    pressures="Pressure ($P$)",
    entropies="Entropy ($A$)",
    conduction_parameters=r"Diffusion ($\alpha$)",
)

current_axis = 0

for key, label in plot.items():
    if current_axis > 4:
        break
    else:
        axis = ax[current_axis]

    try:
        # Raw data
        try:
            grid = (
                project_gas_pixel_grid(sim, resolution=resolution, project=key)
                / unweighted_grid
            )
            axis.imshow(grid, origin="lower", extent=[0, 1, 0, 1], cmap="cividis")
        except:
            continue

        # Exact solution, a square!
        axis.plot(
            [0.25, 0.75, 0.75, 0.25, 0.25],
            [0.25, 0.25, 0.75, 0.75, 0.25],
            linestyle="dashed",
            color="white",
        )

        circle = plt.Circle(
            (0.8, 0.8),
            radius=np.sqrt(1.0 / sim.metadata.n_gas) * 1.25,
            linestyle="solid",
            color="white",
            fill=False,
        )

        axis.add_artist(circle)

        axis.tick_params(
            axis="both",
            which="both",
            labelbottom=False,
            labelleft=False,
            bottom=False,
            left=False,
        )

        axis.set_title(label)

        current_axis += 1
    except KeyError:
        # Mustn't have that data!
        continue


info_axis = ax[-1]

info = get_data_dump(sim.metadata)

info_axis.text(
    0.5, 0.45, info, ha="center", va="center", fontsize=5, transform=info_axis.transAxes
)

info_axis.axis("off")


fig.tight_layout(pad=1.0)
fig.savefig("SquareTest.pdf", dpi=300)
