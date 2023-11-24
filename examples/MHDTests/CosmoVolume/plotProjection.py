#!/usr/bin/env python3

# ------------------------------------------------
# Plots a projection of the DM and baryon mass
# usage:
# $ python3 plotProjection.py <snapnr>
# where <snapnr> is number of snapshot to plot
# ------------------------------------------------

import numpy as np
import sys
import os
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from swiftsimio import load
from swiftsimio.visualisation.projection import project_pixel_grid, project_gas
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)

from unyt import msun, kpc


def set_colorbar(ax, im):
    """
    Adapt the colorbar a bit for axis object <ax> and
    imshow instance <im>
    """
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    return


# Grab snapshot

snapshot_basename = "snapshots/snap"

try:
    snapnr = sys.argv[1]
except IndexError:
    print("You need to provide the index of the snapshot to plot.")
    print("E.g. python3 plotProjection.py 3")
    quit()

try:
    snapnr_int = int(snapnr)
except ValueError:
    print("<snapnr> must be an integer.")
    print("You provided :'" + snapnr + "'")

file = snapshot_basename + "_" + str(snapnr_int).zfill(4) + ".hdf5"

if not os.path.isfile(file):
    print("Didn't find snapshot", file)
    quit()


# Load data
data = load(file)
meta = data.metadata
boxsize = meta.boxsize
extent = [0, boxsize[0].v, 0, boxsize[1].v]

# Generate smoothing lengths for the dark matter
data.dark_matter.smoothing_length = generate_smoothing_lengths(
    data.dark_matter.coordinates,
    data.metadata.boxsize,
    kernel_gamma=1.8,
    neighbours=57,
    speedup_fac=2,
    dimension=3,
)

# Project the dark matter mass
dm_mass = project_pixel_grid(
    # Note here that we pass in the dark matter dataset not the whole
    # data object, to specify what particle type we wish to visualise
    data=data.dark_matter,
    boxsize=data.metadata.boxsize,
    resolution=1024,
    project="masses",
    parallel=True,
    region=None,
)


# Project the gas mass
mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)
mass_map.convert_to_units(msun / kpc ** 2)

######## Magnetic STuff

B = data.gas.magnetic_flux_densities
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
divB = data.gas.magnetic_divergences
h = data.gas.smoothing_lengths

data.gas.B_Mag = normB
data.gas.DivB_error = np.maximum(h * abs(divB) / normB, 1e-6)


divb_map = project_gas(data, resolution=1024, project="DivB_error", parallel=True)
# divb_map.convert_to_units(msun / kpc ** 2)

bfld_map = project_gas(data, resolution=1024, project="B_Mag", parallel=True)
# bfld_map.convert_to_units(msun / kpc ** 2)

# Make figure an plot
# fig = plt.figure(figsize=(12, 5), dpi=200)
fig = plt.figure(figsize=(12, 12), dpi=200)

ax1 = fig.add_subplot(221)
im1 = ax1.imshow(
    dm_mass.T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm()
)
ax1.set_title("Dark Matter Mass")
set_colorbar(ax1, im1)

ax2 = fig.add_subplot(222)
im2 = ax2.imshow(
    mass_map.T,
    origin="lower",
    extent=extent,
    cmap="magma",
    norm=LogNorm(vmax=5e7, vmin=1e5),
)
ax2.set_title("Baryon Mass")
set_colorbar(ax2, im2)

ax3 = fig.add_subplot(223)
im3 = ax3.imshow(
    divb_map.T,
    origin="lower",
    extent=extent,
    cmap="cividis",
    norm=LogNorm(vmax=1e7, vmin=1e5),
)
ax3.set_title("divB")
set_colorbar(ax3, im3)

ax4 = fig.add_subplot(224)
im4 = ax4.imshow(
    bfld_map.T,
    origin="lower",
    extent=extent,
    cmap="magma",
    norm=LogNorm(vmax=1e1, vmin=1e-8),
)
ax4.set_title("Magnetic Field")
set_colorbar(ax4, im4)

# Add xlabels
xunits = data.dark_matter.coordinates.units.latex_representation()
for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlabel("x [" + xunits + "]")
    ax.set_ylabel("y [" + xunits + "]")


# Add title
title = file.replace("_", r"\_")  # exception handle underscore for latex
if meta.cosmology is not None:
    title += ", $z$ = {0:.3f}".format(meta.z)
title += ", $t$ = {0:.2e}".format(meta.time.to("Gyr"))
fig.suptitle(title)

# Save figure
plt.tight_layout()
figname = file[:-5] + ".png"
plt.savefig(figname)
plt.close()

print("saved", figname)
