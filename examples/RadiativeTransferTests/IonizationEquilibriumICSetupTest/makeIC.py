#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
##############################################################################


# -----------------------------------------------------------
# Use particles in with a wide range of different
# internal energies (hence temperatures) to test the
# ionization equilibrium IC setup within GEAR-RT.
# -----------------------------------------------------------

import h5py
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio.metadata.writer.unit_systems import cosmo_units

# define unit system to use
#  unitsystem = unyt.unit_systems.cgs_unit_system
unitsystem = cosmo_units

# set minimal and maximal specific internal energy
umin = 1e9 * unyt.erg / unyt.g
umin = umin.to(unitsystem["energy"] / unitsystem["mass"])
umax = 1e20 * unyt.erg / unyt.g
umax = umax.to(unitsystem["energy"] / unitsystem["mass"])

# Box is 1 Mpc
boxsize = 1 * unyt.kpc
boxsize = boxsize.to(unitsystem["length"])

# number of particles in each dimension
n_p = 15
nparts = n_p**3

# filename of ICs to be generated
outputfilename = "ionization_equilibrium_test.hdf5"

# particle positions
xp = unyt.unyt_array(np.zeros((nparts, 3), dtype=np.float64), boxsize.units)
dx = boxsize / n_p

pmass = 1.0 * unyt.M_Sun
pmass = pmass.to(unitsystem["mass"])

ind = 0
for i in range(n_p):
    x = (i + 0.5) * dx
    for j in range(n_p):
        y = (j + 0.5) * dx
        for k in range(n_p):
            z = (k + 0.5) * dx

            xp[ind] = (x, y, z)
            ind += 1


boxsize_cosmo = sw.cosmo_array(
    [boxsize.value, boxsize.value, boxsize.value],
    boxsize.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=1,
)
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system, boxsize=boxsize_cosmo, dimension=3
)

w.gas.coordinates = sw.cosmo_array(
    xp.value, xp.units, comoving=True, scale_factor=1.0, scale_exponent=1
)
w.gas.velocities = sw.cosmo_array(
    np.zeros(xp.shape),
    unyt.cm / unyt.s,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
w.gas.masses = sw.cosmo_array(
    np.ones(xp.shape[0], dtype=np.float64) * pmass.value,
    pmass.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
ie_vals = np.logspace(np.log10(umin.v), np.log10(umax.v), nparts)
w.gas.internal_energy = sw.cosmo_array(
    ie_vals, umax.units, comoving=True, scale_factor=1.0, scale_exponent=-2
)

# Generate initial guess for smoothing lengths based on MIPS
w.gas.generate_smoothing_lengths()

# If IDs are not present, this automatically generates
w.write(outputfilename)

# Generate initial ionization mass fractions.
# These are deliberately bogus, so you can make sure that
# they have been changed by SWIFT.

F = h5py.File(outputfilename, "r+")
header = F["Header"]
nparts = header.attrs["NumPart_ThisFile"][0]
parts = F["/PartType0"]

# this is bogus data for debugging checks.
HIdata = np.ones(nparts, dtype=np.float32) * 0.11
parts.create_dataset("MassFractionHI", data=HIdata)
HIIdata = np.ones(nparts, dtype=np.float32) * 0.22
parts.create_dataset("MassFractionHII", data=HIIdata)
HeIdata = np.ones(nparts, dtype=np.float32) * 0.123
parts.create_dataset("MassFractionHeI", data=HeIdata)
HeIIdata = np.ones(nparts, dtype=np.float32) * 0.234
parts.create_dataset("MassFractionHeII", data=HeIIdata)
HeIIIdata = np.ones(nparts, dtype=np.float32) * 0.345
parts.create_dataset("MassFractionHeIII", data=HeIIIdata)

F.close()
