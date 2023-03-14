#!/usr/bin/env python3

# ------------------------------------------------
# Plots a projection of the DM and baryon mass
# usage:
# $ python3 plotProjection.py <snapnr>
# where <snapnr> is number of snapshot to plot
# ------------------------------------------------

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


# Make figure an plot
fig = plt.figure(figsize=(12, 5), dpi=200)

ax1 = fig.add_subplot(121)
im1 = ax1.imshow(
    dm_mass.T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm()
)
ax1.set_title("Dark Matter Mass", usetex=True)
set_colorbar(ax1, im1)

ax2 = fig.add_subplot(122)
im2 = ax2.imshow(
    mass_map.T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm()
)
ax2.set_title("Baryon Mass", usetex=True)
set_colorbar(ax2, im2)


# Add xlabels
xunits = data.dark_matter.coordinates.units.latex_representation()
for ax in [ax1, ax2]:
    ax.set_xlabel("x [" + xunits + "]", usetex=True)
    ax.set_ylabel("y [" + xunits + "]", usetex=True)


# Add title
title = file.replace("_", r"\_")  # exception handle underscore for latex
if meta.cosmology is not None:
    title += ", $z$ = {0:.3f}".format(meta.z)
title += ", $t$ = {0:.2e}".format(meta.time.to("Gyr"))
fig.suptitle(title, usetex=True)

# Save figure
plt.tight_layout()
figname = file[:-5] + ".png"
plt.savefig(figname)
plt.close()

print("saved", figname)
