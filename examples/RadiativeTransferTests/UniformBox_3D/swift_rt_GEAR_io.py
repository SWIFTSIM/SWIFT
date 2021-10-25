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
import unyt
import swiftsimio


class RTGasData(object):
    """
    Object to store RT gas particle data of a snapshot
    """

    def __init__(self):

        self.IDs = None
        self.coords = None
        self.h = None

        self.PhotonEnergies = None
        self.PhotonFluxes = None

        return


class RTStarData(object):
    """
    Object to store RT star particle data of a snapshot
    """

    def __init__(self):

        self.IDs = None
        self.coords = None
        self.h = None

        self.InjectedPhotonEnergy = None

        return


class RTSnapData(object):
    """
    Object to store RT snapshot data
    """

    def __init__(self):
        self.snapnr = None
        self.ncells = None
        self.boxsize = None
        self.time = None
        self.has_stars = True

        self.nstars = None
        self.npart = None

        self.stars = RTStarData()
        self.gas = RTGasData()
        return


class Rundata(object):
    """
    Object to store some global swift data
    """

    def __init__(self):
        self.units = None

        self.hydro_controlled_injection = False
        self.use_const_emission_rate = False
        self.has_stars = False  # assume we don't have stars, check while reading in

        self.ngroups = 0  # photon frequency groups
        self.const_emission_rates = None

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

    rundata: Rundata object
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
    firstfile = swiftsimio.load(hdf5files[0])

    try:
        scheme = str(firstfile.metadata.subgrid_scheme["RT Scheme"])
    except KeyError:
        print(
            "These tests only work for the GEAR RT scheme.",
            "Compile swift --with-rt=GEAR_N",
        )
        quit()
    if "GEAR" not in scheme:
        raise ValueError(
            "These tests only work for the GEAR RT scheme.",
            "Compile swift --with-rt=GEAR_N",
        )
        quit()

    ngroups = int(firstfile.metadata.subgrid_scheme["PhotonGroupNumber"])
    rundata.ngroups = ngroups
    rundata.use_const_emission_rate = bool(
        firstfile.metadata.parameters["GEARRT:use_const_emission_rates"]
    )
    rundata.units = firstfile.units

    if rundata.use_const_emission_rate:
        # read emission rate parameter as string
        emissionstr = firstfile.metadata.parameters[
            "GEARRT:star_emission_rates_LSol"
        ].decode("utf-8")
        # clean string up
        if emissionstr.startswith("["):
            emissionstr = emissionstr[1:]
        if emissionstr.endswith("]"):
            emissionstr = emissionstr[:-1]

        # transform string values to floats with unyts
        emissions = emissionstr.split(",")
        emlist = []
        for er in emissions:
            emlist.append(float(er))
        rundata.const_emission_rates = unyt.unyt_array(emlist, unyt.L_Sun)

        if len(rundata.const_emission_rates) != rundata.ngroups:
            print("Got number of emission rates different from number of groups?")
            print(rundata.const_emission_rates, "vs", rundata.ngroups)
            if len(rundata.const_emission_rates) > rundata.ngroups:
                print("Only using first", rundata.ngroups, "emission rates")
                rundata.const_emission_rates = rundata.const_emission_rates[
                    : rundata.ngroups
                ]
            else:
                quit()

    if "hydro controlled" in scheme:
        rundata.hydro_controlled_injection = True

    # -------------------
    # Read in all files
    # -------------------

    for f in hdf5files:
        snapnrstr = f[len(prefix) + 1 : len(prefix) + 5]
        snapnr = int(snapnrstr)

        if skip_snap_zero and snapnr == 0:
            continue
        if skip_last_snap and f == hdf5files[-1]:
            continue

        # Get snap general data
        newsnap = RTSnapData()
        newsnap.snapnr = snapnr

        data = swiftsimio.load(f)
        newsnap.boxsize = data.metadata.boxsize
        newsnap.time = data.metadata.time

        # Get gas data
        Gas = RTGasData()
        Gas.IDs = data.gas.particle_ids
        Gas.coords = data.gas.coordinates
        Gas.h = data.gas.smoothing_lengths
        Gas.PhotonEnergies = swiftsimio.cosmo_array(
            [
                getattr(data.gas.photon_energies, "group" + str(g))
                for g in range(1, rundata.ngroups + 1)
            ]
        ).T

        Gas.PhotonFluxes = swiftsimio.cosmo_array(
            [
                unyt.uvstack(
                    [
                        getattr(data.gas.photon_fluxes, "Group" + str(g) + d)
                        for d in ("X", "Y", "Z")
                    ]
                )
                for g in range(1, rundata.ngroups + 1)
            ],
            data.gas.photon_fluxes.Group1X.units,
        ).T

        newsnap.gas = Gas
        newsnap.npart = Gas.IDs.shape[0]

        #  Get star data
        try:
            Stars = RTStarData()
            Stars.IDs = data.stars.particle_ids
            Stars.coords = data.stars.coordinates
            Stars.h = data.stars.smoothing_lengths
            Stars.InjectedPhotonEnergy = data.stars.rtdebug_injected_photon_energy
            newsnap.stars = Stars
            newsnap.nstars = Stars.IDs.shape[0]
        except AttributeError:
            Stars = RTStarData()
            newsnap.has_stars = False
            newsnap.nstars = 0

        snapdata.append(newsnap)

    for snap in snapdata:
        rundata.has_stars = rundata.has_stars or snap.has_stars

    return snapdata, rundata
