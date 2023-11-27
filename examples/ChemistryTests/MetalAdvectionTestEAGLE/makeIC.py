#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttehove@ugent.be)
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
# Test the diffusion/advection of metals by creating regions with/without
# initial metallicities. Run with EAGLE chemistry.
# ---------------------------------------------------------------------

import numpy as np
import h5py

GAMMA = 5 / 3
RHO = 1
P = 1
VELOCITY = 1
BOX_SIZE = 1

outputfilename = "advect_metals.hdf5"

if __name__ == "__main__":

    glass = h5py.File("glassPlane_128.hdf5", "r")
    parts = glass["PartType0"]
    pos = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()

    # Set up metadata
    boxsize = np.array([1.0, 1.0, 0.0])
    n_part = len(h)
    ids = np.arange(n_part) + 1

    # Setup other particle quantities
    masses = RHO * BOX_SIZE ** 3 / n_part * np.ones(n_part)
    velocities = np.zeros((n_part, 3))
    velocities[:, 0] = VELOCITY
    internal_energy = P / (RHO * (GAMMA - 1)) * np.ones(n_part)

    # Setup metallicities
    metallicities = np.zeros(n_part)
    # Mask for middle square
    mask = ((1 / 3 * BOX_SIZE < pos[:, 0]) & (pos[:, 0] < 2 / 3 * BOX_SIZE) &
            (1 / 3 * BOX_SIZE < pos[:, 1]) & (pos[:, 1] < 2 / 3 * BOX_SIZE))
    metallicities[mask] = 1

    # Now open the file and write the data.
    file = h5py.File(outputfilename, "w")

    # Header
    grp = file.create_group("/Header")
    grp.attrs["BoxSize"] = [1.0, 1.0, 1.0]
    grp.attrs["NumPart_Total"] = [n_part, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [n_part, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = 0
    grp.attrs["Dimension"] = 2

    # Units
    grp = file.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = file.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=velocities, dtype="f")
    grp.create_dataset("Masses", data=masses, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=internal_energy, dtype="f")
    grp.create_dataset("ParticleIDs", data=ids, dtype="L")
    grp.create_dataset("Metallicity", data=metallicities, dtype="L")

    # TODO: add ElementAbundance arrays and IronMassFracFromSNIa (see EAGLE/chemistry_io.h)

    file.close()

    # close up, and we're done!
    file.close()
