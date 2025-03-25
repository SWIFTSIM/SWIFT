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

import numpy as np
import swiftsimio
import unyt


class RTGasData(object):
    """
    Object to store RT gas particle data of a snapshot
    """

    def __init__(self):
        self.IDs = None
        self.coords = None
        self.h = None
        self.volumes = None

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
        self.has_stars = False
        self.has_star_debug_data = False

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

        self.use_const_emission_rate = False
        self.has_stars = False  # assume we don't have stars, check while reading in
        self.has_star_debug_data = (
            False
        )  # assume we don't have stars, check while reading in

        self.ngroups = 0  # photon frequency groups
        self.const_emission_rates = None
        self.reduced_speed_of_light = -1.0

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
        raise ValueError(
            "These tests only work for the GEAR RT scheme.",
            "Compile swift --with-rt=GEAR_N",
        )
    if "GEAR" not in scheme:
        raise ValueError(
            "These tests only work for the GEAR RT scheme.",
            "Compile swift --with-rt=GEAR_N",
        )

    ngroups = int(firstfile.metadata.subgrid_scheme["PhotonGroupNumber"][0])
    rundata.ngroups = ngroups

    luminosity_model = firstfile.metadata.parameters["GEARRT:stellar_luminosity_model"]
    rundata.use_const_emission_rate = luminosity_model.decode("utf-8") == "const"

    rundata.units = firstfile.units

    if rundata.use_const_emission_rate:
        # read emission rate parameter as string
        emissionstr = firstfile.metadata.parameters[
            "GEARRT:const_stellar_luminosities_LSol"
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
            print(
                "paramfile:",
                len(rundata.const_emission_rates),
                "vs",
                rundata.ngroups,
                "groups",
            )
            if len(rundata.const_emission_rates) > rundata.ngroups:
                rundata.const_emission_rates = rundata.const_emission_rates[
                    : rundata.ngroups
                ]
                print(
                    "Only using first",
                    rundata.ngroups,
                    "emission rates:",
                    rundata.const_emission_rates,
                )
            else:
                quit()
    else:
        print(
            "Didn't detect use of constant stellar emission rates. Proceeding without."
        )

    rundata.reduced_speed_of_light = firstfile.metadata.reduced_lightspeed

    with_mpi = False
    if firstfile.metadata.code["MPI library"] != b"Non-MPI version of SWIFT":
        with_mpi = True
    rundata.with_mpi = with_mpi

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
        masses = data.gas.masses
        densities = data.gas.densities
        # skip potential div by zero
        mask = densities == 0.0
        masses[mask] = 0.0
        densities[mask] = 1.0
        Gas.volumes = masses / densities

        Gas.PhotonEnergies = [
            getattr(data.gas.photon_energies, "group" + str(g + 1))
            for g in range(rundata.ngroups)
        ]
        Gas.PhotonFluxes = [
            unyt.uvstack(
                [
                    getattr(data.gas.photon_fluxes, "Group" + str(g + 1) + d)
                    for d in ("X", "Y", "Z")
                ]
            ).T
            for g in range(rundata.ngroups)
        ]

        newsnap.gas = Gas
        newsnap.npart = Gas.IDs.shape[0]

        #  Get star data
        Stars = RTStarData()
        nstars = 0
        has_stars = False
        has_star_debug_data = False
        try:
            Stars.IDs = data.stars.particle_ids
            Stars.coords = data.stars.coordinates
            Stars.h = data.stars.smoothing_lengths
            nstars = Stars.IDs.shape[0]
            has_stars = True
        except AttributeError:
            pass

        if has_stars:
            try:
                inj = np.atleast_2d(data.stars.rtdebug_injected_photon_energy)
                Stars.InjectedPhotonEnergy = np.reshape(
                    inj, (Stars.IDs.shape[0], ngroups)
                )
                has_star_debug_data = True
            except AttributeError:
                pass

        newsnap.stars = Stars
        newsnap.nstars = nstars
        newsnap.has_stars = has_stars
        newsnap.has_star_debug_data = has_star_debug_data

        snapdata.append(newsnap)

    for snap in snapdata:
        rundata.has_stars = rundata.has_stars or snap.has_stars
        rundata.has_star_debug_data = (
            rundata.has_star_debug_data or snap.has_star_debug_data
        )

    return snapdata, rundata
