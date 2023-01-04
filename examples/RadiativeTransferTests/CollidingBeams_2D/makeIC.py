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
outputfilename = "collision_2D.hdf5"


def get_radiation_IC(pos, numPart, ngroups):

    Elist = []
    Flist = []

    x = pos[:, 0]
    y = pos[:, 1]
    c = unyt.c.to(unitsystem["length"] / unitsystem["time"])

    for g in range(ngroups):

        E = np.zeros(numPart)
        F = np.zeros((numPart, 3))

        # prepare particle array masks
        flux_sign = np.ones(numPart, dtype=float)

        if g == 0:
            # horizontal beams
            # -----------------------

            left = np.logical_and(x > 0.0, x < 1.0 / 3)
            right = np.logical_and(x > 2.0 / 3, x < 1.0)
            xcriterion = np.logical_or(left, right)
            ycriterion = np.logical_and(y > 0.4, y < 0.6)
            is_max = np.logical_and(xcriterion, ycriterion)

            flux_sign[right] = -1.0

            Emax = 1.0
            E[is_max] = Emax
            F[is_max, 0] = Emax * c * flux_sign[is_max]

        if g == 1:
            # vertical beams, nonzero energy everywhere
            # -------------------------------------------

            top = np.logical_and(y > 2.0 / 3.0, y < 1.0)
            bottom = np.logical_and(y > 0.0, y < 1.0 / 3.0)

            xcriterion = np.logical_and(x > 0.4, x < 0.6)
            ycriterion = np.logical_or(top, bottom)
            is_max = np.logical_and(xcriterion, ycriterion)

            flux_sign[y > 0.5] = -1.0

            Emax = 2.0
            Emin = 1.0
            E[:] = Emin
            E[is_max] = Emax
            F[:, 1] = E * c * flux_sign

        if g == 2:
            # diagonal beams
            # -------------------------------------------

            width = 0.1
            l = np.sqrt(2) / 3.0
            sin = np.sin(np.pi / 4)
            xprime = l * sin
            x0 = xprime - 0.5 * width * sin
            x1 = xprime + 0.5 * width * sin
            x2 = 1 - xprime - 0.5 * width * sin
            x3 = 1 - xprime + 0.5 * width * sin

            def upper_line(x):
                return x + 0.5 * width

            def lower_line(x):
                return x - 0.5 * width

            def descending_left(x, xprime):
                return -(x - xprime) + xprime

            def descending_right(x, xprime):
                return -(x - 1.0 + xprime) + 1.0 - xprime

            first = x < x0
            lower_first = np.logical_and(y < upper_line(x), y > lower_line(x))
            firstcond = np.logical_and(first, lower_first)

            second = np.logical_and(x > x0, x < x1)
            lower_second = np.logical_and(
                y > lower_line(x), y < descending_left(x, xprime)
            )
            secondcond = np.logical_and(second, lower_second)

            third = np.logical_and(x > x2, x < x3)
            upper_third = np.logical_and(
                y < upper_line(x), y > descending_right(x, xprime)
            )
            thirdcond = np.logical_and(third, upper_third)

            fourth = np.logical_and(x > x3, x < 1.0)
            upper_fourth = np.logical_and(y < upper_line(x), y > lower_line(x))
            fourthcond = np.logical_and(fourth, upper_fourth)

            Emax = 1.0
            E[firstcond] = Emax
            E[secondcond] = Emax
            E[thirdcond] = Emax
            E[fourthcond] = Emax

            flux_sign[thirdcond] = -1.0
            flux_sign[fourthcond] = -1.0

            F[:, 0] = E * c * flux_sign / np.sqrt(2)
            F[:, 1] = E * c * flux_sign / np.sqrt(2)

        if g == 3:
            # diagonal beams that meed in the middle
            # -------------------------------------------

            width = 0.1
            l = np.sqrt(2) / 3.0
            sin = np.sin(np.pi / 4)
            xprime = l * sin
            x0 = xprime - 0.5 * width * sin
            x1 = xprime + 0.5 * width * sin
            x2 = 1 - xprime - 0.5 * width * sin
            x3 = 1 - xprime + 0.5 * width * sin

            def upper_line(x):
                return x + 0.5 * width

            def lower_line(x):
                return x - 0.5 * width

            def descending_left(x, xprime):
                return -(x - xprime) + xprime

            def descending_lower(x, xprime):
                return -(x - 1 + xprime) + xprime - 0.5 * width

            def descending_upper(x, xprime):
                return -(x - 1.0 + xprime) + xprime + 0.5 * width

            def ascending_right(x, xprime):
                return (x - 1 + xprime) + xprime

            first = x < x0
            lower_first = np.logical_and(y < upper_line(x), y > lower_line(x))
            firstcond = np.logical_and(first, lower_first)

            second = np.logical_and(x > x0, x < x1)
            lower_second = np.logical_and(
                y > lower_line(x), y < descending_left(x, xprime)
            )
            secondcond = np.logical_and(second, lower_second)

            third = np.logical_and(x > x2, x < x3)
            upper_third = np.logical_and(
                y < ascending_right(x, xprime), y > descending_lower(x, xprime)
            )
            thirdcond = np.logical_and(third, upper_third)

            fourth = np.logical_and(x > x3, x < 1.0)
            upper_fourth = np.logical_and(
                y < descending_upper(x, xprime), y > descending_lower(x, xprime)
            )
            fourthcond = np.logical_and(fourth, upper_fourth)

            Emax = 1.0
            E[firstcond] = Emax
            E[secondcond] = Emax
            E[thirdcond] = Emax
            E[fourthcond] = Emax

            flux_sign[thirdcond] = -1.0
            flux_sign[fourthcond] = -1.0

            F[:, 0] = E * c * flux_sign / np.sqrt(2)
            F[:, 1] = E * c / np.sqrt(2)

        #  histE, _, _ = np.histogram2d(pos[:,0], pos[:,1], weights=E, bins=50)
        #  histFx, _, _ = np.histogram2d(pos[:,0], pos[:,1], weights=F[:,0], bins=50)
        #  histFy, _, _ = np.histogram2d(pos[:,0], pos[:,1], weights=F[:,1], bins=50)
        #  from matplotlib import pyplot as plt
        #
        #  fig = plt.figure()
        #  ax1 = fig.add_subplot(1, 3, 1)
        #  ax1.imshow(histE.T, origin="lower")
        #  ax2 = fig.add_subplot(1, 3, 2)
        #  ax2.imshow(histFx.T, origin="lower")
        #  ax3 = fig.add_subplot(1, 3, 3)
        #  ax3.imshow(histFy.T, origin="lower")
        #  plt.show()

        Elist.append(E)
        Flist.append(F)

    return Elist, Flist


if __name__ == "__main__":
    glass = h5py.File("glassPlane_128.hdf5", "r")

    # Read particle positions and h from the glass
    pos = glass["/PartType0/Coordinates"][:, :]
    h = glass["/PartType0/SmoothingLength"][:]
    glass.close()

    numPart = np.size(h)

    # get radiation IC while particle coordinates are still [0, 1)
    Elist, Flist = get_radiation_IC(pos, numPart, nPhotonGroups)

    pos *= boxsize
    h *= boxsize

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

    for g in range(nPhotonGroups):
        Esetname = "PhotonEnergiesGroup{0:d}".format(g + 1)
        parts[Esetname][:] = Elist[g][:]
        Fsetname = "PhotonFluxesGroup{0:d}".format(g + 1)
        parts[Fsetname][:, :] = Flist[g][:, :]

    F.close()
