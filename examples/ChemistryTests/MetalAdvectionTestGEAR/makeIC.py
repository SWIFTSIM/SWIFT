#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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
import math

# ---------------------------------------------------------------------
# Test the diffusion/advection of metals by creating regions with/without
# initial metallicities. Run with EAGLE chemistry.
# ---------------------------------------------------------------------

import numpy as np
import h5py

GAMMA = 5 / 3
RHO = 1
P = 1
VELOCITY = 2.5
# Set ELEMENT_COUNT to 2 for AGORA runs, or for GEAR compile swift with `--with-chemistry=GEAR_X` with X equal to
# ELEMENT_COUNT.
ELEMENT_COUNT = 4

outputfilename = "advect_metals.hdf5"


def get_mask(boxsize, pos, element_idx):
    x = boxsize[0]
    right_half_mask = pos[:, 0] > x / 3

    if element_idx == ELEMENT_COUNT - 1:
        return right_half_mask
    else:
        periods = 2 * (element_idx + 2)
        periods_mask = (pos[:, 0] * periods // x) % 2 == 0
        return right_half_mask & periods_mask


def get_abundance(element_idx):
    if element_idx == ELEMENT_COUNT - 1:
        return 0.2
    else:
        return 0.1 * 0.5 ** element_idx


def get_element_abundances_metallicity(pos, boxsize):
    element_abundances = []
    for i in range(ELEMENT_COUNT):
        element_abundances.append(
            np.where(get_mask(boxsize, pos, i), get_abundance(i), 0)
        )

    return np.stack(element_abundances, axis=1)


if __name__ == "__main__":
    glass = h5py.File("glassPlane_64.hdf5", "r")
    parts = glass["PartType0"]
    pos = parts["Coordinates"][:]
    pos = np.concatenate([pos, pos + np.array([1, 0, 0])])
    h = parts["SmoothingLength"][:]
    h = np.concatenate([h, h])
    glass.close()

    # Set up metadata
    boxsize = np.array([2.0, 1.0, 1.0])
    n_part = len(h)
    ids = np.arange(n_part) + 1

    # Setup other particle quantities
    rho = RHO * np.ones_like(h)
    rho[pos[:, 1] < 0.5 * boxsize[1]] *= 0.5
    masses = rho * np.prod(boxsize) / n_part
    velocities = np.zeros((n_part, 3))
    velocities[:, :] = 0.5 * math.sqrt(2) * VELOCITY * np.array([1.0, 1.0, 0.0])
    internal_energy = P / (rho * (GAMMA - 1))

    # Setup metallicities
    metallicities = get_element_abundances_metallicity(pos, boxsize)

    # Now open the file and write the data.
    file = h5py.File(outputfilename, "w")

    # Header
    grp = file.create_group("/Header")
    grp.attrs["BoxSize"] = boxsize
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
    grp.create_dataset("MetalMassFraction", data=metallicities, dtype="f")

    file.close()

    # close up, and we're done!
    file.close()

    if ELEMENT_COUNT == 2:
        print(
            f"Make sure swift was compiled with `--with-chemistry=GEAR_2` or `--with-chemistry=AGORA`"
        )
    else:
        print(
            f"Make sure swift was compiled with `--with-chemistry=GEAR_{ELEMENT_COUNT}`"
        )
