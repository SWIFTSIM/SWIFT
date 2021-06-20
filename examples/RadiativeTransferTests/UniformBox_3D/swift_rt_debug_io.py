#!/usr/bin/env python3

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

        self.RTStarIact = None
        self.RTCallsIactGradient = None
        self.RTCallsIactTransport = None
        self.RTCallsIactGradientInteraction = None
        self.RTCallsIactTransportInteraction = None

        self.InjectionDone = None
        self.ThermochemistryDone = None
        self.TransportDone = None
        self.GradientsDone = None

        self.RadiationAbsorbedTot = None

        return


class RTStarData(object):
    """
    Object to store RT star particle data of a snapshot
    """

    def __init__(self):

        self.IDs = None
        self.coords = None
        self.h = None

        self.RTHydroIact = None
        self.EmissionRateSet = None

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
        self.stars = RTStarData()
        self.gas = RTGasData()
        return


class Rundata(object):
    """
    Object to store some global swift data
    """

    def __init__(self):
        self.hydro_controlled_injection = False

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
        print(
            "These tests only work for the debug RT scheme.",
            "Compile swift --with-rt=debug",
        )
        F.close()
        quit()

    if "debug" not in scheme and "GEAR" not in scheme:
        raise ValueError(
            "These tests only work for the debug RT scheme.",
            "Compile swift --with-rt=debug",
        )

    if "hydro controlled" in scheme:
        rundata.hydro_controlled_injection = True
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

        newsnap.gas.RTStarIact = Gas["RTDebugStarIact"][:][inds]
        newsnap.gas.RTCallsIactGradient = Gas["RTDebugCallsIactGradient"][:][inds]
        newsnap.gas.RTCallsIactTransport = Gas["RTDebugCallsIactTransport"][:][inds]
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

        Stars = F["PartType4"]
        ids = Stars["ParticleIDs"][:]
        inds = np.argsort(ids)

        newsnap.stars.IDs = ids[inds]
        newsnap.stars.coords = Stars["Coordinates"][:][inds]
        newsnap.stars.h = Stars["SmoothingLengths"][:][inds]

        newsnap.stars.RTHydroIact = Stars["RTDebugHydroIact"][:][inds]
        newsnap.stars.EmissionRateSet = Stars["RTDebugEmissionRateSet"][:][inds]

        newsnap.stars.RadiationEmittedTot = Stars["RTDebugRadEmittedTot"][:][inds]

        snapdata.append(newsnap)

    return snapdata, rundata
