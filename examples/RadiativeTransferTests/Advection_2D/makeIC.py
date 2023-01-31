#!/usr/bin/env python3

###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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


# -------------------------------------------------------------
# Add initial conditions for photon energies and fluxes
# for 2D advection of photons.
# First photon group: Top hat function with zero as the
#       baseline, advects in x direction
# Second photon group: Top hat function with nonzero value
#       as the baseline, advcts in y direction.
# Third photon group: Gaussian advecting diagonally
# Fourth photon group: Circle moving radially from the center
# -------------------------------------------------------------

import h5py
import numpy as np
import unyt
from swiftsimio import Writer

# define unit system to use
unitsystem = unyt.unit_systems.cgs_unit_system

# define box size
boxsize = 1e10 * unitsystem["length"]

# number of photon groups
nPhotonGroups = 4

# filename of ICs to be generated
outputfilename = "advection_2D.hdf5"


def initial_condition(x):
    """
    The initial conditions that will be advected

    x: particle position. 3D unyt array

    returns: 
    E: photon energy for each photon group. List of scalars with size of nPhotonGroups
    F: photon flux for each photon group. List with size of nPhotonGroups of numpy arrays of shape (3,)
    """

    # you can make the photon quantities unitless, the units will
    # already have been written down in the writer.

    E_list = []
    F_list = []
    c = unyt.c.to(unitsystem["length"] / unitsystem["time"])

    # Group 1 Photons:
    # -------------------

    in_x = 0.33 * boxsize < x[0] < 0.66 * boxsize
    in_y = 0.33 * boxsize < x[1] < 0.66 * boxsize
    if in_x and in_y:
        E = 1.0
    else:
        E = 0.0

    # Assuming all photons flow in only one direction
    # (optically thin regime, "free streaming limit"),
    #  we have that |F| = c * E
    F = np.zeros(3, dtype=np.float64)
    F[0] = c * E

    E_list.append(E)
    F_list.append(F)

    # Group 2 Photons:
    # -------------------

    in_x = 0.33 * boxsize < x[0] < 0.66 * boxsize
    in_y = 0.33 * boxsize < x[1] < 0.66 * boxsize
    if in_x and in_y:
        E = 2.0
    else:
        E = 1.0

    F = np.zeros(3, dtype=np.float64)
    F[1] = c * E

    E_list.append(E)
    F_list.append(F)

    # Group 3 Photons:
    # -------------------
    sigma = 0.1 * boxsize
    mean = 0.5 * boxsize
    amplitude = 2.0
    baseline = 1.0

    E = (
        amplitude
        * np.exp(-((x[0] - mean) ** 2 + (x[1] - mean) ** 2) / (2 * sigma ** 2))
        + baseline
    )
    F = np.zeros(3, dtype=np.float64)
    F[0] = c * E / 1.414213562  # sqrt(2)
    F[1] = c * E / 1.414213562  # sqrt(2)

    E_list.append(E)
    F_list.append(F)

    # Group 4 Photons:
    # -------------------

    circle_radius = 0.15 * boxsize
    center = 0.5 * boxsize
    dx = x[0] - center
    dy = x[1] - center
    r = np.sqrt(dx ** 2 + dy ** 2)
    if r <= circle_radius:
        unit_vector = (dx / r, dy / r)

        E = 1.0
        F = np.zeros(3, dtype=np.float64)
        F[0] = unit_vector[0] * c * E
        F[1] = unit_vector[1] * c * E

    else:
        E = 0.0
        F = np.zeros(3, dtype=np.float64)

    E_list.append(E)
    F_list.append(F)

    return E_list, F_list


if __name__ == "__main__":
    glass = h5py.File("glassPlane_128.hdf5", "r")

    # Read particle positions and h from the glass
    pos = glass["/PartType0/Coordinates"][:, :]
    h = glass["/PartType0/SmoothingLength"][:]
    glass.close()

    pos *= boxsize
    h *= boxsize

    numPart = np.size(h)

    w = Writer(unyt.unit_systems.cgs_unit_system, boxsize, dimension=2)

    w.gas.coordinates = pos
    w.gas.velocities = np.zeros((numPart, 3)) * (unyt.cm / unyt.s)
    w.gas.masses = np.ones(numPart, dtype=np.float64) * 1000 * unyt.g
    w.gas.internal_energy = (
        np.ones(numPart, dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / unyt.g
    )

    # Generate initial guess for smoothing lengths based on MIPS
    w.gas.smoothing_length = h

    # If IDs are not present, this automatically generates
    w.write(outputfilename)

    # Now open file back up again and add photon groups
    # you can make them unitless, the units have already been
    # written down in the writer. In this case, it's in cgs.

    F = h5py.File(outputfilename, "r+")
    header = F["Header"]
    nparts = header.attrs["NumPart_ThisFile"][0]
    parts = F["/PartType0"]

    for grp in range(nPhotonGroups):
        dsetname = "PhotonEnergiesGroup{0:d}".format(grp + 1)
        energydata = np.zeros(nparts, dtype=np.float32)
        parts.create_dataset(dsetname, data=energydata)

        dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
        #  if dsetname not in parts.keys():
        fluxdata = np.zeros((nparts, 3), dtype=np.float32)
        parts.create_dataset(dsetname, data=fluxdata)

    for p in range(nparts):
        E, Flux = initial_condition(pos[p])
        for g in range(nPhotonGroups):
            Esetname = "PhotonEnergiesGroup{0:d}".format(g + 1)
            parts[Esetname][p] = E[g]
            Fsetname = "PhotonFluxesGroup{0:d}".format(g + 1)
            parts[Fsetname][p] = Flux[g]

    F.close()
