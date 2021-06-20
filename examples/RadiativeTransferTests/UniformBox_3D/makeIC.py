#!/usr/bin/env python3

# ---------------------------------------------------------------------
# Create a uniform grid of hydro particles and a smaller box of
# star particles in the center of the box. The main idea here is that
# there are some hydro particles that have no star neighbours.
# ---------------------------------------------------------------------

from swiftsimio import Writer

import unyt
import numpy as np
import h5py

# Box is 1 Mpc
boxsize = 100 * unyt.m

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


w = Writer(unyt.unit_systems.cgs_unit_system, boxsize)

w.gas.coordinates = xp
w.stars.coordinates = xs
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
w.stars.velocities = np.zeros(xs.shape) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(xp.shape[0], dtype=np.float) * 1000 * unyt.g
w.stars.masses = np.ones(xs.shape[0], dtype=np.float) * 1000 * unyt.g
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=np.float) * (300.0 * unyt.kb * unyt.K) / (unyt.g)
)

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)
w.stars.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
w.write("uniformBox-rt.hdf5")


# Writing Photon Initial Conditions
# --------------------------------------

# This part is specific for the GEAR RT scheme, assuming we're
# running with nPhotonGroups photon groups.

nPhotonGroups = 4

# Now open file back up again and add photon groups
# you can make them unitless, the units have already been
# written down in the writer. In this case, it's in cgs.

F = h5py.File("uniformBox-rt.hdf5", "r+")
header = F["Header"]
nparts = header.attrs["NumPart_ThisFile"][0]
nstars = header.attrs["NumPart_ThisFile"][4]
parts = F["/PartType0"]


for grp in range(nPhotonGroups):
    dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
    energydata = np.ones((nparts), dtype=np.float) * (grp + 1)
    parts.create_dataset(dsetname, data=energydata)

    dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
    #  if dsetname not in parts.keys():
    fluxdata = np.ones((nparts, 3), dtype=np.float) * (grp + 1)
    fluxdata[:, 1] *= 2.0
    fluxdata[:, 2] *= 3.0
    parts.create_dataset(dsetname, data=fluxdata)

F.close()
