###############################################################################
# This file is part of the ANARCHY paper.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#               2019 Josh Borrow (joshua.boorrow@durham.ac.uk)
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
################################################################################

# Compares the swift result for the 2D spherical Sod shock with a high
# resolution 2D reference result

import sys
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from swiftsimio import load
from analyticSolution import analytic

snap = int(sys.argv[1])

sim = load(f"sodshock_{snap:04d}.hdf5")

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


# See analyticSolution for params.


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


# Read the simulation data
time = sim.metadata.t.value

data = dict(
    x=sim.gas.coordinates.value[:, 0],
    v=sim.gas.velocities.value[:, 0],
    u=sim.gas.internal_energy.value,
    S=sim.gas.entropy.value,
    P=sim.gas.pressure.value,
    rho=sim.gas.density.value,
    y=sim.gas.coordinates.value[:, 1],
    z=sim.gas.coordinates.value[:, 2],
)

# Try to add on the viscosity and diffusion.
try:
    data["visc"] = sim.gas.viscosity.value
except:
    pass

try:
    data["diff"] = 100.0 * sim.gas.diffusion.value
except:
    pass

# Read in the "solution" data and calculate those that don't exist.

ref = analytic(time=time)

# We only want to plot this for the region that we actually have data for, hence the masking.
mask = np.logical_and(ref["x"] < np.max(data["x"]), ref["x"] > np.min(data["x"]))
ref = {k: v[mask] for k, v in ref.items()}

# Now we can do the plotting.
fig, ax = plt.subplots(2, 3, figsize=(6.974, 6.974 * (2.0 / 3.0)))
ax = ax.flatten()

# These are stored in priority order
plot = dict(
    v="Velocity ($v_x$)",
    u="Internal Energy ($u$)",
    rho=r"Density ($\rho$)",
    P="Pressure ($P$)",
    diff=r"100$\times$ Diffusion Coefficient ($\alpha_D$)",
    visc=r"Viscosity Coefficient ($\alpha_V$)",
    S="Entropy ($A$)",
)

log = dict(v=False, u=False, S=False, P=False, rho=False, visc=False, diff=False)
ylim = dict(v=(-0.05, 1.0), diff=(0.0, None), visc=(0.0, None))

current_axis = 0

for key, label in plot.items():
    if current_axis > 4:
        break
    else:
        axis = ax[current_axis]

    try:
        if log[key]:
            axis.semilogy()

        # Raw data
        axis.plot(
            data["x"],
            data[key],
            ".",
            color="C1",
            markersize=0.5,
            alpha=0.5,
            rasterized=True,
            markeredgecolor="none",
            zorder=-1,
        )

        mask_noraster = np.logical_and.reduce([
            data["y"] < 0.52,
            data["y"] > 0.48,
            data["z"] < 0.52,
            data["z"] > 0.48
        ])

        axis.plot(
            data["x"][mask_noraster],
            data[key][mask_noraster],
            ".",
            color="C3",
            rasterized=False,
            markeredgecolor="none",
            markersize=3,
            zorder=0,
        )

        # Exact solution
        try:
            axis.plot(ref["x"], ref[key], c="C0", ls="dashed", zorder=1, lw=1)
        except KeyError:
            # No solution :(
            pass

        axis.set_xlabel("Position ($x$)", labelpad=0)
        axis.set_ylabel(label, labelpad=0)

        axis.set_xlim(0.6, 1.5)

        try:
            axis.set_ylim(*ylim[key])
        except KeyError:
            # No worries pal
            pass

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


fig.tight_layout(pad=0.5)
fig.savefig("SodShock.pdf", dpi=300)
