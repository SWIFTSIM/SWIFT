#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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


# --------------------------------------------------
# 'Module' containing RT I/O routines for the RT
# debugging scheme
# --------------------------------------------------


import os
import h5py
import numpy as np


class RTGasData(object):
    """
    Object to store RT gas particle data of a snapshot
    """

    def __init__(self):

        self.IDs = None
        self.coords = None
        self.h = None

        self.RTCallsIactGradientInteraction = None
        self.RTCallsIactTransportInteraction = None

        self.InjectionDone = None
        self.ThermochemistryDone = None
        self.TransportDone = None
        self.GradientsDone = None

        self.RadiationAbsorbedTot = None

        self.nsubcycles = None

        return


class RTStarData(object):
    """
    Object to store RT star particle data of a snapshot
    """

    def __init__(self):

        self.IDs = None
        self.coords = None
        self.h = None

        self.EmissionRateSet = None
        self.InjectionInteractions = None
        self.RadiationEmittedTot = None

        return


class RTSnapData(object):
    """
    Object to store RT snapshot data
    """

    def __init__(self):
        self.snapnr = None
        self.ncells = None
        self.boxsize = None
        self.has_stars = True
        self.stars = RTStarData()
        self.gas = RTGasData()
        return


class Rundata(object):
    """
    Object to store some global swift data
    """

    def __init__(self):
        self.has_stars = False  # assume we don't have stars, check while reading in
        self.with_mpi = False

        return


def get_snap_data(prefix="output", skip_snap_zero=False, skip_last_snap=False):
    """
    Finds all prefix_XXXX.hdf5 files and reads
    the RT data in.

    Parameters
    ----------

    prefix: str
        file name prefix of snapshots. Don't include the
        underscore!

    skip_snap_zero: bool
        whether to skip the snapshot prefix_0000.hdf5

    skip_last_snap: bool
        whether to skip the last available snapshot


    Returns
    -------

    snapdata: list
        list of RTSnapData objects filled out with actual
        snapshot data

    rundata: rundata object
        object containing data of the entire run/simulation/compilation
    """

    snapdata = []

    ls = os.listdir()
    hdf5files = []
    for f in ls:
        if f.startswith(prefix + "_") and f.endswith(".hdf5"):
            hdf5files.append(f)

    if len(hdf5files) == 0:
        raise IOError("No " + prefix + "_XXXX.hdf5 files found in this directory")

    hdf5files.sort()

    # --------------------------------------
    # Get global data from first snapshot
    # --------------------------------------

    rundata = Rundata()
    firstfile = hdf5files[0]
    F = h5py.File(firstfile, "r")

    try:
        scheme = str(F["SubgridScheme"].attrs["RT Scheme"])
    except KeyError:
        raise ValueError(
            "These tests only work for the debug RT scheme.",
            "Compile swift --with-rt=debug",
        )

    if "debug" not in scheme and "GEAR" not in scheme:
        raise ValueError(
            "These tests only work for the debug RT scheme.",
            "Compile swift --with-rt=debug",
        )

    with_mpi = False
    mpistr = F["Code"].attrs["MPI library"]
    if mpistr != b"Non-MPI version of SWIFT":
        with_mpi = True
    rundata.with_mpi = with_mpi

    F.close()

    for f in hdf5files:
        snapnrstr = f[len(prefix) + 1 : len(prefix) + 5]
        snapnr = int(snapnrstr)

        if skip_snap_zero and snapnr == 0:
            continue
        if skip_last_snap and f == hdf5files[-1]:
            continue

        newsnap = RTSnapData()
        newsnap.snapnr = snapnr

        F = h5py.File(f, "r")
        newsnap.boxsize = F["Header"].attrs["BoxSize"]
        newsnap.ncells = F["Cells"]
        Gas = F["PartType0"]
        ids = Gas["ParticleIDs"][:]
        inds = np.argsort(ids)

        newsnap.gas.IDs = ids[inds]
        newsnap.gas.coords = Gas["Coordinates"][:][inds]
        newsnap.gas.h = Gas["SmoothingLengths"][:][inds]

        newsnap.gas.RTCallsIactGradientInteraction = Gas[
            "RTDebugCallsIactGradientInteractions"
        ][:][inds]
        newsnap.gas.RTCallsIactTransportInteraction = Gas[
            "RTDebugCallsIactTransportInteractions"
        ][:][inds]
        newsnap.gas.InjectionDone = Gas["RTDebugInjectionDone"][:][inds]
        newsnap.gas.GradientsDone = Gas["RTDebugGradientsDone"][:][inds]
        newsnap.gas.TransportDone = Gas["RTDebugTransportDone"][:][inds]
        newsnap.gas.ThermochemistryDone = Gas["RTDebugThermochemistryDone"][:][inds]

        newsnap.gas.RadiationAbsorbedTot = Gas["RTDebugRadAbsorbedTot"][:][inds]
        newsnap.gas.nsubcycles = Gas["RTDebugSubcycles"][:][inds]

        try:
            Stars = F["PartType4"]
            ids = Stars["ParticleIDs"][:]
            has_stars = True

            inds = np.argsort(ids)

            newsnap.stars.IDs = ids[inds]
            newsnap.stars.coords = Stars["Coordinates"][:][inds]
            newsnap.stars.h = Stars["SmoothingLengths"][:][inds]

            newsnap.stars.EmissionRateSet = Stars["RTDebugEmissionRateSet"][:][inds]
            newsnap.stars.InjectionInteractions = Stars["RTDebugHydroIact"][:][inds]
            newsnap.stars.RadiationEmittedTot = Stars["RTDebugRadEmittedTot"][:][inds]

        except KeyError:
            has_stars = False

        newsnap.has_stars = has_stars
        snapdata.append(newsnap)

    for snap in snapdata:
        rundata.has_stars = rundata.has_stars or snap.has_stars

    if len(snapdata) == 0:
        print("Didn't read in snapshot data.")
        print("Do you only have 2 snapshots and are skipping the first and the last?")
        quit()

    return snapdata, rundata
