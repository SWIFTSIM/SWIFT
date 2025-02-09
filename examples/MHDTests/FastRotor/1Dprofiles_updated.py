import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# Where we want to slice the slice
y0 = 0.5 #0.5

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
rho = data.gas.densities

P = data.gas.pressures

B = data.gas.magnetic_flux_densities

Bx, By = B[:, 0], B[:, 1]

normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)

divB = data.gas.magnetic_divergences

h = data.gas.smoothing_lengths

errB = np.log10(h * abs(divB) / normB)

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_densities = data.gas.masses * rho

data.gas.mass_weighted_pressures = data.gas.masses * P

data.gas.mass_weighted_Bx = data.gas.masses * Bx

data.gas.mass_weighted_By = data.gas.masses * By

data.gas.mass_weighted_errB = data.gas.masses * errB

common_arguments = dict(
    data=data, z_slice=0.5 * data.metadata.boxsize[2], resolution=512, parallel=True
)

mass_map = slice_gas(**common_arguments, project="masses")

mass_weighted_density_map = slice_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_pressure_map = slice_gas(
    **common_arguments, project="mass_weighted_pressures"
)

mass_weighted_Bx_map = slice_gas(**common_arguments, project="mass_weighted_Bx")

mass_weighted_By_map = slice_gas(**common_arguments, project="mass_weighted_By")

mass_weighted_errB_map = slice_gas(**common_arguments, project="mass_weighted_errB")

# Take out mass dependence
density_map = mass_weighted_density_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map
Bx_map = mass_weighted_Bx_map / mass_map
By_map = mass_weighted_By_map / mass_map
errB_map = mass_weighted_errB_map / mass_map

map_pixel_length = len(mass_map)

x = np.linspace(0.0, 1.0, map_pixel_length)
slice_ind = int(np.floor(y0 * map_pixel_length))

# Plot maps
plt.rcParams.update({"font.size": 16})

fig, axs = plt.subplots(5, 1, figsize=((2*4, 4*4)), sharex=True)
fig.subplots_adjust(hspace=0.1)


#axs[0].plot(x, density_map[:, slice_ind], "k-", lw=0.5)
#axs[0].set_yticks(np.arange(2.0, 14.0, 2.0))
#axs[0].set_ylabel(r"$\rho$")
#axs[0].set_ylim(0.0, 13.0)
#axs[0].plot(x, pressure_map[:, slice_ind], "k-", lw=0.5)
#axs[0].set_yticks(np.arange(0.0, 0.35, 0.05))
#axs[0].set_ylabel(r"$P$")
#axs[0].set_ylim(0.0, 0.3)

#axs[1].plot(x, Bx_map[:, slice_ind], "k-", lw=0.5)
#axs[1].set_yticks(np.arange(-0.5, 0.55, 0.25))
#axs[1].set_ylabel(r"$B_x$")
#axs[1].set_ylim(-0.5, 0.5)

#axs[2].plot(x, By_map[:, slice_ind], "k-", lw=0.5)
#axs[2].set_yticks(np.arange(-0.5, 0.55, 0.25))
#axs[2].set_ylabel(r"$B_y$")
#axs[2].set_ylim(-0.5, 0.5)

#axs[3].plot(x, errB_map[:, slice_ind], "k-", lw=0.5)
#axs[3].set_yticks(np.arange(-6.0, 1.0, 1.0))
#axs[3].set_xlabel(r"$x$")
#axs[3].set_ylabel(r"$\mathrm{log}_{10} \left( h \quad \nabla \cdot B / |B| \right)$")
#axs[3].set_ylim(-6.5, 0.5)
axs[0].plot(x, density_map[:, slice_ind], "k-", lw=0.5)
axs[0].set_yticks(np.arange(2.0, 14.0, 2.0))
axs[0].set_ylabel(r"$\rho$")
axs[0].set_ylim(0.0, 13.0)


axs[1].plot(x, Bx_map[:, slice_ind], "k-", lw=0.5)
axs[1].set_yticks(np.arange(0.0, 2.0, 0.5))
axs[1].set_ylabel(r"$B_x$")
axs[1].set_ylim(-0.4, 1.6)


axs[2].plot(x, By_map[:, slice_ind], "k-", lw=0.5)
axs[2].set_yticks(np.arange(-1.0, 1.5, 0.5))
axs[2].set_ylabel(r"$B_y$")
axs[2].set_ylim(-1.2, 1.4)

axs[3].plot(x, errB_map[:, slice_ind], "k-", lw=0.5)
axs[3].set_yticks(np.arange(-6.0, 1.0, 1.0))
axs[3].set_xlabel(r"$x$")
axs[3].set_ylabel(r"$\mathrm{log}_{10} \left( h \quad \nabla \cdot B / |B| \right)$")
axs[3].set_ylim(-6.5, 0.5)


# Add panel with infromation about the run
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=axs[4].transAxes
)

axs[4].text(
    0.5,
    0.8,
    "Orszag Tang Vortex at $t=%.2f$, $y_0 = %.4f$" % (data.metadata.time, y0),
    **text_common_args,
)
axs[4].text(0.5, 0.7, "SWIFT %s" % git.decode("utf-8"), **text_common_args)
axs[4].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
axs[4].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
axs[4].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
axs[4].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
axs[4].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
axs[4].text(
    0.5, 0.1, r"Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
)
axs[4].axis("off")



for ax in axs:
    ax.minorticks_on()
    ax.grid()

fig.tight_layout()
plt.savefig(args.output)
