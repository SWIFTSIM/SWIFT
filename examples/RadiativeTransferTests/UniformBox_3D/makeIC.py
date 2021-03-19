#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Create a uniform grid of hydro particles and a smaller box of
# star particles in the center of the box. The main idea here is that
# there are some hydro particles that have no star neighbours.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np

# Box is 1 Mpc
boxsize = 1 * unyt.Mpc

# number of hydro particles in each dimension
n_p = 40

# number of star particles in each dimension
n_s = 20

xp = []
xs = []

dx = boxsize / n_p
ds = 0.2 * boxsize / n_s

# Generate hydro coordinates
for i in range(n_p):
    x = (i + 0.5) * dx
    for j in range(n_p):
        y = (j + 0.5) * dx
        for k in range(n_p):
            z = (k + 0.5) * dx
            xp.append(np.array([x, y, z], dtype=np.float))

# Generate star coordinates
for i in range(n_s):
    # factor 0.52 below: shift particles a bit so they don't overlap with hydro
    x = 0.4 * boxsize + (i + 0.52) * ds
    for j in range(n_s):
        y = 0.4 * boxsize + (j + 0.52) * ds
        for k in range(n_s):
            z = 0.4 * boxsize + (k + 0.52) * ds
            xs.append(np.array([x, y, z], dtype=np.float))

xp = unyt.unyt_array(xp, boxsize.units)
xs = unyt.unyt_array(xs, boxsize.units)


# Generate object. cosmo_units corresponds to default Gadget-oid units
# of 10^10 Msun, Mpc, and km/s
w = Writer(cosmo_units, boxsize)

w.gas.coordinates = xp
w.stars.coordinates = xs


# Random velocities from 0 to 1 km/s
w.gas.velocities = np.zeros(xp.shape) * (unyt.km / unyt.s)
w.stars.velocities = np.zeros(xs.shape) * (unyt.km / unyt.s)

# Generate uniform masses as 10^6 solar masses for each particle
w.gas.masses = np.ones(xp.shape[0], dtype=float) * (1e6 * unyt.msun)
w.stars.masses = np.ones(xs.shape[0], dtype=float) * (1e6 * unyt.msun)

# Generate internal energy corresponding to 10^4 K
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=float) * (1e4 * unyt.kb * unyt.K) / (1e6 * unyt.msun)
)

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
w.stars.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write("uniformBox-rt.hdf5")
