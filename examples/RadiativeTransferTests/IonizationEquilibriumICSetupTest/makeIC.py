#!/usr/bin/env python3

# -----------------------------------------------------------
# Use 10k particles in 1D with a wide range of different
# internal energies (hence temperatures) to test the
# ionization equilibrium IC setup within GEAR-RT.
# -----------------------------------------------------------

from swiftsimio import Writer

import unyt
import numpy as np
import h5py

# define unit system to use
unitsystem = unyt.unit_systems.cgs_unit_system

# set minimal and maximal specific internal energy
umin = 1e9 * unyt.erg / unyt.g
umax = 1e20 * unyt.erg / unyt.g

# Box is 1 Mpc
boxsize = 1e15 * unitsystem["length"]

# number of photon groups
nPhotonGroups = 1

# number of particles in each dimension
n_p = 10000

# filename of ICs to be generated
outputfilename = "ionization_equilibrium_test.hdf5"

# particle positions
xp = unyt.unyt_array(np.zeros((n_p, 3), dtype=np.float64), boxsize.units)
dx = boxsize / n_p
for i in range(n_p):
    xp[i, 0] = (i + 0.5) * dx

w = Writer(unyt.unit_systems.cgs_unit_system, boxsize, dimension=1)

w.gas.coordinates = xp
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * 1000 * unyt.g
w.gas.internal_energy = (
    np.logspace(np.log10(umin.v), np.log10(umax.v), n_p) * unyt.erg / unyt.g
)

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=1)

# If IDs are not present, this automatically generates
w.write(outputfilename)

# Generate initial ionization mass fractions.
# These are deliberately bogus so you can make sure that
# they have been changed by SWIFT.

F = h5py.File(outputfilename, "r+")
header = F["Header"]
nparts = header.attrs["NumPart_ThisFile"][0]
parts = F["/PartType0"]

# this is bogus data for debugging checks.
HIdata = np.ones((nparts), dtype=np.float32) * 0.11
parts.create_dataset("MassFractionHI", data=HIdata)
HIIdata = np.ones((nparts), dtype=np.float32) * 0.22
parts.create_dataset("MassFractionHII", data=HIIdata)
HeIdata = np.ones((nparts), dtype=np.float32) * 0.123
parts.create_dataset("MassFractionHeI", data=HeIdata)
HeIIdata = np.ones((nparts), dtype=np.float32) * 0.234
parts.create_dataset("MassFractionHeII", data=HeIIdata)
HeIIIdata = np.ones((nparts), dtype=np.float32) * 0.345
parts.create_dataset("MassFractionHeIII", data=HeIIIdata)

F.close()
