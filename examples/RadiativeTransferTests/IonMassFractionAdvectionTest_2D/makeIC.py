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


# ---------------------------------------------------------------------
# Test the diffusion/advection of ions by creating regions with/without
# initial ions. Run without radiation.
# ---------------------------------------------------------------------

from swiftsimio import Writer
from swiftsimio.units import cosmo_units
import unyt
import numpy as np
import h5py

gamma = 5.0 / 3.0
outputfilename = "advect_ions.hdf5"

if __name__ == "__main__":

    glass = h5py.File("glassPlane_64.hdf5", "r")
    parts = glass["PartType0"]
    xp = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()

    # Set up metadata
    unitL = unyt.Mpc
    edgelen = 1 * unyt.Mpc
    edgelen = edgelen.to(unitL)
    boxsize = np.array([1.0, 1.0, 0.0]) * edgelen

    xp *= edgelen
    h *= edgelen

    w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=2)

    # write particle positions and smoothing lengths
    w.gas.coordinates = xp
    w.gas.smoothing_length = h

    # get gas masses
    mpart = 1.6e5 * unyt.Msun
    masses = np.ones(xp.shape[0], dtype=np.float64) * mpart
    # change some gas masses
    mask = xp[:, 0] > 0.5 * edgelen
    masses[mask] *= 3
    w.gas.masses = masses

    w.gas.internal_energy = (
        np.ones(xp.shape[0], dtype=np.float64) * 1.25e6 * unyt.m ** 2 / unyt.s ** 2
    )

    # get velocities
    vels = np.zeros((xp.shape[0], 3))
    vels[:, 0] = -1.0
    vels[:, 1] = +1.0
    w.gas.velocities = vels * 1000 * cosmo_units["length"] / cosmo_units["time"]

    w.write(outputfilename)

    # Now open file back up again and add RT data.
    F = h5py.File(outputfilename, "r+")
    header = F["Header"]
    nparts = header.attrs["NumPart_ThisFile"][0]
    parts = F["/PartType0"]
    pos = parts["Coordinates"]

    # Create initial ionization species mass fractions.
    HIdata = np.ones(nparts, dtype=np.float32) * 0.76
    HIIdata = np.zeros(nparts, dtype=np.float32)
    HeIdata = np.ones(nparts, dtype=np.float32) * 0.24
    HeIIdata = np.zeros(nparts, dtype=np.float32)
    HeIIIdata = np.zeros(nparts, dtype=np.float32)

    mask1 = pos[:, 0] > edgelen / 3
    mask1 = np.logical_and(mask1, pos[:, 0] < 2 * edgelen / 3)
    HIdata[mask1] = 0.26
    HIIdata[mask1] = 0.5

    mask2 = pos[:, 1] > edgelen / 4
    mask2 = np.logical_and(mask2, pos[:, 1] < edgelen / 2)
    HeIdata[mask2] = 0.1
    HeIIdata[mask2] = 0.14

    mask3 = pos[:, 1] > edgelen / 2
    mask3 = np.logical_and(mask3, pos[:, 1] < 3 * edgelen / 4)
    HeIdata[mask3] = 0.05
    HeIIdata[mask3] = 0.09
    HeIIIdata[mask3] = 0.1

    parts.create_dataset("MassFractionHI", data=HIdata)
    parts.create_dataset("MassFractionHII", data=HIIdata)
    parts.create_dataset("MassFractionHeI", data=HeIdata)
    parts.create_dataset("MassFractionHeII", data=HeIIdata)
    parts.create_dataset("MassFractionHeIII", data=HeIIIdata)

    # close up, and we're done!
    F.close()
