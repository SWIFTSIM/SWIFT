import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# Load snapshot
filename = args.input
data = load(filename)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]

# Retrieve particle attributes of interest
v = data.gas.velocities
normv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
B = data.gas.magnetic_flux_densities
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
divB = data.gas.magnetic_divergences
h = data.gas.smoothing_lengths

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities

data.gas.mass_weighted_pressures = data.gas.masses * data.gas.pressures

data.gas.mass_weighted_normv = data.gas.masses * normv

data.gas.mass_weighted_normB = data.gas.masses * normB

data.gas.mass_weighted_error = data.gas.masses * np.log10(
    np.maximum(h * abs(divB) / (normB + 0.01 * np.max(normB)), 1e-6)
)

common_arguments = dict(data=data, resolution=512, parallel=True)

mass_map = project_gas(**common_arguments, project="masses")

mass_weighted_density_map = project_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_pressure_map = project_gas(
    **common_arguments, project="mass_weighted_pressures"
)

mass_weighted_normv_map = project_gas(**common_arguments, project="mass_weighted_normv")

mass_weighted_normB_map = project_gas(**common_arguments, project="mass_weighted_normB")

mass_weighted_error_map = project_gas(**common_arguments, project="mass_weighted_error")

# Take out mass dependence
density_map = mass_weighted_density_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map
normv_map = mass_weighted_normv_map / mass_map
normB_map = mass_weighted_normB_map / mass_map
error_map = mass_weighted_error_map / mass_map

# Plot maps
plt.rcParams.update({"font.size": 16})
fig, ax = plt.subplots(3, 2, figsize=(12.25, 17))

a00 = ax[0, 0].contourf(
    density_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 12.0, 100)
)
a01 = ax[0, 1].contourf(
    pressure_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 3.0, 100)
)
a10 = ax[1, 0].contourf(
    normv_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 1.5, 100)
)
a11 = ax[1, 1].contourf(
    normB_map.value.T, cmap="gist_heat", levels=np.linspace(0.0, 2.5, 100)
)
a20 = ax[2, 0].contourf(error_map.value.T, cmap="jet", levels=np.linspace(-5.0, 0.0, 6))

# Add panel with infromation about the run
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=ax[2, 1].transAxes
)

ax[2, 1].text(
    0.5,
    0.8,
    "Orszag Tang Vortex at time $t=%.2f$" % data.metadata.time,
    **text_common_args,
)
ax[2, 1].text(0.5, 0.7, "SWIFT %s" % git.decode("utf-8"), **text_common_args)
ax[2, 1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
ax[2, 1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
ax[2, 1].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
ax[2, 1].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
ax[2, 1].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
ax[2, 1].text(
    0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
)
ax[2, 1].axis("off")

for axi in ax:
    for axii in axi:
        axii.set_xticks([])
        axii.set_yticks([])
        axii.set_aspect("equal")

# Set appropriate colourbars
fig.colorbar(
    a00,
    ax=ax[0, 0],
    label=r"$\rho$",
    fraction=0.042,
    pad=0.04,
    location="left",
    ticks=np.linspace(0.0, 12.0, 7),
)

fig.colorbar(
    a01,
    ax=ax[0, 1],
    label=r"$P$",
    fraction=0.042,
    pad=0.04,
    ticks=np.linspace(0.0, 3.0, 7),
)

fig.colorbar(
    a10,
    ax=ax[1, 0],
    label=r"$|\mathbf{v}|$",
    fraction=0.042,
    pad=0.04,
    location="left",
    ticks=np.linspace(0.0, 1.5, 7),
)

fig.colorbar(
    a11,
    ax=ax[1, 1],
    label=r"$|\mathbf{B}|$",
    fraction=0.042,
    pad=0.04,
    ticks=np.linspace(0.0, 2.5, 6),
)

fig.colorbar(
    a20,
    ax=ax[2, 0],
    label=r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$",
    fraction=0.042,
    pad=0.04,
    location="left",
)

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig(args.output)
