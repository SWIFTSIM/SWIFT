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


# -----------------------------------------------------------
# Add initial conditions for photon energies and fluxes
# for 1D advection of photons.
# First photon group: Top hat function with zero as the
#       baseline.
# Second photon group: Top hat function with nonzero value
#       as the baseline.
# Third photon group: Gaussian.
# -----------------------------------------------------------

import h5py
import numpy as np
import unyt
from swiftsimio import Writer

# define unit system to use
unitsystem = unyt.unit_systems.cgs_unit_system

# Box is 1 Mpc
boxsize = 1e10 * unitsystem["length"]

# number of photon groups
nPhotonGroups = 3

# number of particles in each dimension
n_p = 1000

# filename of ICs to be generated
outputfilename = "collision_1D.hdf5"


def initial_condition(x):
    """
    The initial conditions that will be advected

    x: particle position. 3D unyt array

    returns: 
    E: photon energy density for each photon group. List of scalars with size of nPhotonGroups
    F: photon flux for each photon group. List with size of nPhotonGroups of numpy arrays of shape (3,)
    """

    # you can make the photon quantities unitless, the units will
    # already have been written down in the writer.

    E_list = []
    F_list = []

    # Group 1 Photons:
    # -------------------

    F = np.zeros(3, dtype=np.float64)

    if x[0] < 0.33 * boxsize:
        E = 1.0
        F[0] = unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E
    elif x[0] < 0.66 * boxsize:
        E = 0.0
        F[0] = unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E
    else:
        E = 1.0
        F[0] = -unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E

    E_list.append(E)
    F_list.append(F)

    # Group 2 Photons:
    # -------------------

    F = np.zeros(3, dtype=np.float64)
    if x[0] < 0.33 * boxsize:
        E = 3.0
        F[0] = unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E
    elif x[0] < 0.66 * boxsize:
        E = 1.0
        if x[0] < 0.5 * boxsize:
            F[0] = unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E
        else:
            F[0] = -unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E
    else:
        E = 3.0
        F[0] = -unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E

    E_list.append(E)
    F_list.append(F)

    # Group 3 Photons:
    # -------------------
    sigma = 0.05 * boxsize
    if x[0] < 0.5 * boxsize:
        mean = 0.33 / 2 * boxsize
        sign = 1
    else:
        mean = (5.0 / 6.0) * boxsize
        sign = -1
    amplitude = 2.0

    E = amplitude * np.exp(-(x[0] - mean) ** 2 / (2 * sigma ** 2))
    F = np.zeros(3, dtype=np.float64)
    F[0] = sign * unyt.c.to(unitsystem["length"] / unitsystem["time"]) * E

    E_list.append(E)
    F_list.append(F)

    return E_list, F_list


if __name__ == "__main__":

    xp = unyt.unyt_array(np.zeros((n_p, 3), dtype=np.float64), boxsize.units)

    dx = boxsize / n_p

    for i in range(n_p):
        xp[i, 0] = (i + 0.5) * dx

    w = Writer(unyt.unit_systems.cgs_unit_system, boxsize, dimension=1)

    w.gas.coordinates = xp
    w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
    w.gas.masses = np.ones(xp.shape[0], dtype=np.float64) * 1000 * unyt.g
    w.gas.internal_energy = (
        np.ones(xp.shape[0], dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / unyt.g
    )

    # Generate initial guess for smoothing lengths based on MIPS
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=1)

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
        energydata = np.zeros((nparts), dtype=np.float32)
        parts.create_dataset(dsetname, data=energydata)

        dsetname = "PhotonFluxesGroup{0:d}".format(grp + 1)
        #  if dsetname not in parts.keys():
        fluxdata = np.zeros((nparts, 3), dtype=np.float32)
        parts.create_dataset(dsetname, data=fluxdata)

    for p in range(nparts):
        E, Flux = initial_condition(xp[p])
        for g in range(nPhotonGroups):
            Esetname = "PhotonEnergiesGroup{0:d}".format(g + 1)
            parts[Esetname][p] = E[g]
            Fsetname = "PhotonFluxesGroup{0:d}".format(g + 1)
            parts[Fsetname][p] = Flux[g]

    #  from matplotlib import pyplot as plt
    #  plt.figure()
    #  for g in range(nPhotonGroups):
    #      Esetname = "PhotonEnergiesGroup{0:d}".format(g+1)
    #      plt.plot(xp[:,0], parts[Esetname], label="E "+str(g+1))
    #      #  Fsetname = "PhotonFluxesGroup{0:d}".format(g+1)
    #      #  plt.plot(xp[:,0], parts[Fsetname][:,0], label="F "+str(g+1))
    #  plt.legend()
    #  plt.show()

    F.close()
