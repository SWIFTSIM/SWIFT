#!/usr/bin/env python3

# ------------------------------------------------
# Plots a projection of the DM mass.
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

from swiftsimio import load
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)


# Grab snapshot

snapshot_basename = "snap"

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


# Make figure an plot
fig = plt.figure(figsize=(6, 5), dpi=200)
ax = fig.add_subplot(111)
im = ax.imshow(dm_mass.T, origin="lower", extent=extent, cmap="inferno", norm=LogNorm())
plt.colorbar(im)

# Add xlabels
xunits = data.dark_matter.coordinates.units.latex_representation()
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
