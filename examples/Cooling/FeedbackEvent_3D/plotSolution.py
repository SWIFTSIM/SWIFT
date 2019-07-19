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

import sys
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from unyt import cm, s, km, kpc, Pa, msun, K, keV, mh

kPa = 1000 * Pa

from swiftsimio import load

snap = int(sys.argv[1])

sim = load(f"feedback_{snap:04d}.hdf5")

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


# Read the simulation data
boxSize = sim.metadata.boxsize[0]

x = sim.gas.coordinates[:, 0] - boxSize / 2
y = sim.gas.coordinates[:, 1] - boxSize / 2
z = sim.gas.coordinates[:, 2] - boxSize / 2
r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
vel = sim.gas.velocities
v_r = (x * vel[:, 0] + y * vel[:, 1] + z * vel[:, 2]) / r

# Remove unyt information
v_r = v_r.to(km / s).value
r = r.to(kpc).value


data = dict(
    x=r,
    v_r=v_r,
    u=sim.gas.temperature.to(K).value,
    S=sim.gas.entropy.to(keV / K).value,
    P=sim.gas.pressure.to(kPa).value,
    rho=sim.gas.density.to(mh / (cm ** 3)).value,
)

# Try to add on the viscosity and diffusion.
try:
    data["visc"] = sim.gas.viscosity.value
except:
    pass

try:
    data["diff"] = sim.gas.diffusion.value
except:
    pass

# Bin the data
x_bin_edge = np.linspace(0.0, boxSize.to(kpc).value)
x_bin = 0.5 * (x_bin_edge[1:] + x_bin_edge[:-1])
binned = {
    k: stats.binned_statistic(data["x"], v, statistic="mean", bins=x_bin_edge)[0]
    for k, v in data.items()
}
square_binned = {
    k: stats.binned_statistic(data["x"], v ** 2, statistic="mean", bins=x_bin_edge)[0]
    for k, v in data.items()
}
sigma = {
    k: np.sqrt(v2 - v ** 2)
    for k, v2, v in zip(binned.keys(), square_binned.values(), binned.values())
}

# Now we can do the plotting.
fig, ax = plt.subplots(2, 3, figsize=(6.974, 6.974 * (2.0 / 3.0)))
ax = ax.flatten()

# These are stored in priority order
plot = dict(
    v_r="Radial Velocity ($v_r$) [km s$^{-1}$]",
    u="Temperature ($T$) [K]",
    rho=r"Density ($\rho$) [cm$^{-1}$]",
    visc=r"Viscosity Coefficient ($\alpha_V$)",
    diff=r"Diffusion Coefficient ($\alpha_D$)",
    P="Pressure ($P$) [kPa]",
    S="Entropy ($A$) [keV K$^{-1}$]",
)

log = dict(
    v_r=False, v_phi=False, u=False, S=False, P=False, rho=False, visc=False, diff=False
)
ylim = dict(
    v_r=[-8, 5], u=[3500, 5500], rho=[0.02, 0.15], visc=[0, 2.0], diff=[0, 0.25],
    P=[3e-18, 10e-18], S=[-0.5e60, 4e60] 
)

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
            ms=0.5,
            alpha=0.5,
            markeredgecolor="none",
            rasterized=True,
            zorder=0,
        )
        # Binned data
        axis.errorbar(
            x_bin,
            binned[key],
            yerr=sigma[key],
            fmt=".",
            ms=3.0,
            color="C3",
            lw=0.5,
            zorder=2,
        )

        axis.set_xlabel("Radius ($r$) [kpc]", labelpad=0)
        axis.set_ylabel(label, labelpad=0)

        axis.set_xlim(0.0, 0.7 * boxSize.to(kpc).value)

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
fig.savefig("FeedbackEvent.pdf", dpi=300)
