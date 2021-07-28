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


# ------------------------------------------------------------------------
# Collection of sanity checks for the 'debug' RT scheme in swift.
# Any run done with the RT debug scheme must pass these tests.
#
# Usage:
#   ./rt_sanity_checks.py
# or
#   ./rt_sanity_checks.py snapshot_basename
#
# where 'snapshot_basename' is the base name of the snapshots that you
# want to work with. (Same argument as Snapshots:basename in the .yml
# file). If 'snapshot_basename' is not given, this script assumes it to
# be 'output'.
#
# The script will then read all the snapshot_basename_XXXX.hdf5 files
# in the directory where you're running the script from
# ------------------------------------------------------------------------


import numpy as np
from sys import argv
from swift_rt_debug_io import get_snap_data


# some behaviour options
skip_snap_zero = True  # skip snap_0000.hdf5
skip_last_snap = True  # skip snap_0max.hdf5
print_diffs = True  # print differences you find
break_on_diff = False  # quit when you find a difference

hydro_controlled_injection = False


if len(argv) > 1:
    file_prefix = argv[1]
else:
    file_prefix = "output"


def check_hydro_sanity(snapdata):
    """
    Sanity checks for hydro variables.
    - injection always done?
    - gradients always done?
    - thermochemistry always done?
    - RT transport calls >= RT gradient calls?
    """

    nsnaps = len(snapdata)
    npart = snapdata[0].gas.coords.shape[0]

    print("Checking hydro")

    # ----------------------------------------------
    # check absolute values of every snapshot
    # ----------------------------------------------
    for snap in snapdata:

        gas = snap.gas
        stars = snap.stars

        # has a particle been called at least once?
        called = (
            gas.RTCallsIactTransportInteraction[:]
            + gas.RTCallsIactGradientInteraction[:]
            + gas.InjectionDone[:]
            + gas.GradientsDone[:]
            + gas.TransportDone[:]
            + gas.ThermochemistryDone[:]
        ) > 0

        # --------------------------------------------------------------
        # check that photons have been updated (ghost1 called)
        # --------------------------------------------------------------

        mask = gas.InjectionDone != 1
        fishy = np.logical_and(mask, called)
        if fishy.any():
            # has particle been active in the meantime?
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            if np.count_nonzero(mask) == npart:
                print("--- WARNING: zero particles finished injection")
            else:
                print(
                    "--- Some photons have injection finished != 1: ",
                    np.count_nonzero(called),
                    "/",
                    npart,
                )
                if print_diffs:
                    print("----- IDs with photons_updated != 1:")
                    print(gas.IDs[fishy], gas.InjectionDone[fishy])

            if break_on_diff:
                quit()

        # --------------------------------------------------------------
        # check that Gradient is finished
        # --------------------------------------------------------------
        mask = gas.GradientsDone != 1
        fishy = np.logical_and(mask, called)
        if fishy.any():
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            if np.count_nonzero(mask) == npart:
                print("---WARNING: zero particles finished gradient")
            else:
                print(
                    "--- Some gradients were finalized != 1",
                    np.count_nonzero(called),
                    "/",
                    npart,
                )
                if print_diffs:
                    print("----- IDs with gradients done != 1:")
                    print(gas.IDs[fishy], gas.GradientsDone[fishy])

            if break_on_diff:
                quit()

        # --------------------------------------------------------------
        # check that transport is finished
        # --------------------------------------------------------------
        mask = gas.TransportDone != 1
        fishy = np.logical_and(mask, called)
        if fishy.any():
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            if np.count_nonzero(mask) == npart:
                print("--- WARNING: zero particles finished transport")
            else:
                print(
                    "--- Some transport was finalised != 1",
                    np.count_nonzero(called),
                    "/",
                    npart,
                )
                if print_diffs:
                    print("----- IDs with transport done != 1:")
                    print(gas.IDs[fishy], gas.TransportDone[fishy])

            if break_on_diff:
                quit()

        # --------------------------------------------------------------
        # check that thermochemistry is finished
        # --------------------------------------------------------------
        mask = gas.ThermochemistryDone != 1
        fishy = np.logical_and(mask, called)
        if fishy.any():
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            if np.count_nonzero(mask) == npart:
                print("--- WARNING: zero particles finished thermochemistry")
            else:
                print(
                    "--- Some thermochemistry done != 1",
                    np.count_nonzero(called),
                    "/",
                    npart,
                )
                if print_diffs:
                    print("----- IDs with Thermochemistry_done != 1:")
                    print(gas.IDs[fishy], gas.ThermochemistryDone[fishy])

            if break_on_diff:
                quit()

        # --------------------------------------------------------------
        # check that number of calls to gradient interactions is
        # at least the number of calls to transport interactions
        # in RT interactions
        # --------------------------------------------------------------
        if (
            gas.RTCallsIactTransportInteraction < gas.RTCallsIactGradientInteraction
        ).any():
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            print(
                "--- Found RT transport calls iact < gradient calls iact:",
                np.count_nonzero(
                    gas.RTCallsIactTransport < gas.RTCallsIactGradientInteraction
                ),
                "/",
                npart,
            )
            if break_on_diff:
                quit()
        if (gas.RTCallsIactTransport < gas.RTCallsIactGradient).any():
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            print(
                "--- Found RT transport calls < gradient calls:",
                np.count_nonzero(
                    gas.RTCallsIactTransport < gas.RTCallsIactGradientInteraction
                ),
                "/",
                npart,
            )
            if break_on_diff:
                quit()

        # --------------------------------------------------------------
        # check that we didn't loose any radiation
        # --------------------------------------------------------------
        sum_gas_tot = gas.RadiationAbsorbedTot.sum()
        if snap.has_stars:
            sum_star_tot = stars.RadiationEmittedTot.sum()
        else:
            sum_star_tot = 0.0

        if sum_gas_tot != sum_star_tot:
            print("- checking hydro sanity pt2; snapshot", snap.snapnr)
            print(
                "--- Total emitted and absorbed radiation not equal: Gas",
                sum_gas_tot,
                "stars",
                sum_star_tot,
                "diff",
                sum_star_tot - sum_gas_tot,
            )
            if break_on_diff:
                quit()

    return


def check_stars_sanity(snapdata):
    """
    Sanity checks for stars variables.
    - total calls keep increasing?
    """

    print("Checking stars")

    nsnaps = len(snapdata)

    # ----------------------------------------------
    #  check consistency of individual snapshots
    # ----------------------------------------------
    for snap in snapdata:
        if not snap.has_stars:
            continue

        this = snap.stars
        nspart = snapdata[0].stars.coords.shape[0]

        if hydro_controlled_injection:
            if (this.EmissionRateSet != 1).any():
                print("- checking stars sanity pt2", snap.snapnr)
                print("--- Emisison Rates not consistent")
                count = 0
                for i in range(nspart):
                    if this.EmissionRateSet[i] != 1:
                        count += 1
                        if print_diffs:
                            print("-----", this.EmissionRateSet[i], "ID", this.IDs[i])

                print("--- count", count, "/", this.EmissionRateSet.shape[0])

                if break_on_diff:
                    quit()

    return


def main():
    """
    Main function to run.
    """

    snapdata, rundata = get_snap_data(
        prefix=file_prefix, skip_snap_zero=skip_snap_zero, skip_last_snap=skip_last_snap
    )
    global hydro_controlled_injection
    hydro_controlled_injection = rundata.hydro_controlled_injection

    check_hydro_sanity(snapdata)
    check_stars_sanity(snapdata)


if __name__ == "__main__":
    main()
