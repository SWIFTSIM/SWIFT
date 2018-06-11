#!/usr/bin/env python3

import yt

# TODO
# Gas:
# Surface density vs R
# Average vertical height vs R
# Rotational velocity vs R
# Velocity dispersion vs R
# Temperature vs density
# mass distribution as function of density
# mass distribution difference between SN+SFR?
# Newly stars formed vs R
# Clumps count vs R

# Do the same for stars
# Star formation rate
# and more...

filename = "./agora_disk_0000.hdf5"

snap = yt.load(filename, over_refine_factor=2)

density = snap.all_data()[("gas", "density")]
x = snap.all_data()[("gas", "x")]

projection_plot = yt.ProjectionPlot(
    snap,
    "z",
    ("gas", "cell_mass"),
    width=5.5
).show()
