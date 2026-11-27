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
# initial metallicities. Run with GEAR or AGORA chemistry.
# ---------------------------------------------------------------------

import numpy as np
import h5py
import unyt
import swiftsimio as sw

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
        return 0.1 * 0.5**element_idx


def get_element_abundances_metallicity(pos, boxsize):
    element_abundances = []
    for i in range(ELEMENT_COUNT):
        element_abundances.append(
            np.where(get_mask(boxsize, pos, i), get_abundance(i), 0)
        )

    return np.stack(element_abundances, axis=1)


if __name__ == "__main__":
    glass = sw.load("glassPlane_64.hdf5")
    parts = glass.gas
    pos = parts.coordinates.value
    pos = np.concatenate([pos, pos + np.array([1, 0, 0])])
    h = parts.smoothing_lengths.value
    h = np.concatenate([h, h])

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
    w = sw.Writer(
        unit_system=unyt.unit_systems.cgs_unit_system,
        boxsize=sw.cosmo_array(
            boxsize[:2], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
        ),
        dimension=2,
    )

    w.gas.coordinates = pos * unyt.cm * w.a.comoving
    w.gas.velocities = velocities * unyt.cm / unyt.s * w.a.comoving**0
    w.gas.masses = masses * unyt.g * w.a.comoving**0
    w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
    w.gas.internal_energy = internal_energy * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
    w.gas.particle_ids = ids.astype(np.uint64)

    w.write(outputfilename)

    # Add the metallicities
    with h5py.File(outputfilename, "r+") as f:
        grp = f["PartType0"]
        grp.create_dataset("MetalMassFraction", data=metallicities, dtype="f")

    if ELEMENT_COUNT == 2:
        print(
            "Make sure swift was compiled with `--with-chemistry=GEAR_2` or `--with-chemistry=AGORA`"
        )
    else:
        print(
            f"Make sure swift was compiled with `--with-chemistry=GEAR_{ELEMENT_COUNT}`"
        )
