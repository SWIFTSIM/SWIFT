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


# -----------------------------------------------------------------------
# Collection of checks for the 'debug' RT scheme in swift for the
# uniform box test where particles don't move and every time step an
# output file is generated.
#
# Usage:
#   ./rt_uniform_box_checks.py
# or
#   ./rt_uniform_box_checks.py snapshot_basename
#
# where 'snapshot_basename' is the base name of the snapshots that you
# want to work with. (Same argument as Snapshots:basename in the .yml
# file). If 'snapshot_basename' is not given, this script assumes it to
# be 'output'.
#
# The script will then read all the snapshot_basename_XXXX.hdf5 files
# in the directory where you're running the script from
# -----------------------------------------------------------------------


from sys import argv

import numpy as np

from swift_rt_debug_io import get_snap_data

# some behaviour options
skip_snap_zero = True  # skip snap_0000.hdf5
skip_last_snap = True  # skip snap_0max.hdf5
skip_coords = True  # skip coordinates check
skip_sml = True  # skip smoothing lengths check
print_diffs = False  # print differences you find
break_on_diff = False  # quit when you find a difference

# tolerance for a float to be equal
float_comparison_tolerance = 1e-5


if len(argv) > 1:
    file_prefix = argv[1]
else:
    file_prefix = "output"


def check_all_hydro_is_equal(snapdata):
    """
    Check that all the hydro quantities are equal in every snapshot
    (for the relevant quantities of course.)
    """

    ref = snapdata[0]
    npart = ref.gas.coords.shape[0]

    print("checking hydro")

    for compare in snapdata[1:]:

        # Smoothing length ratios
        sml_diff = np.abs(1.0 - ref.gas.h / compare.gas.h)
        # are smoothing lengths within tolerance?
        sml_outside_tolerance = sml_diff > float_comparison_tolerance

        # Coordinates
        if not skip_coords:

            diff = np.abs((ref.gas.coords - compare.gas.coords) / ref.gas.coords)
            if (diff > float_comparison_tolerance).any():
                print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
                print("--- Coordinates vary")
                if print_diffs:
                    for i in range(ref.gas.coords.shape[0]):
                        if (
                            (ref.gas.coords[i] - compare.gas.coords[i])
                            / ref.gas.coords[i]
                        ).any():
                            print(ref.gas.coords[i], "|", compare.gas.coords[i])

                if break_on_diff:
                    quit()

        # Smoothing Lengths
        if not skip_sml:

            diff = np.abs((ref.gas.h - compare.gas.h) / ref.gas.h)
            if (diff > float_comparison_tolerance).any():
                print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
                print("--- Smoothing Lengths vary")
                if print_diffs:
                    for i in range(npart):
                        if ((ref.gas.h[i] - compare.gas.h[i]) / ref.gas.h[i]).any():
                            print(ref.gas.h[i], "|", compare.gas.h[i])

                if break_on_diff:
                    quit()

        # Photon number updates
        if (ref.gas.InjectionDone != compare.gas.InjectionDone).any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            print("--- Calls to InjectionDone")

            if print_diffs:
                for i in range(npart):
                    if ref.gas.InjectionDone[i] != compare.gas.InjectionDone[i]:
                        print(
                            "-----",
                            ref.gas.IDs[i],
                            ref.gas.InjectionDone[i],
                            compare.gas.InjectionDone[i],
                        )

            if break_on_diff:
                quit()

        # Gradient Loop Interaction Calls
        fishy = (
            ref.gas.RTCallsIactGradientInteraction
            != compare.gas.RTCallsIactGradientInteraction
        )
        if fishy.any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            account_for_sml_diff = np.count_nonzero(
                np.logical_and(fishy, sml_outside_tolerance)
            )
            print(
                "--- Calls to iact gradient: count differ: {0:8d} / {1:8d};".format(
                    np.count_nonzero(fishy), npart
                ),
                "After removing ones with acceptable h differences: {0:8d}".format(
                    account_for_sml_diff
                ),
            )

            if account_for_sml_diff > 0:
                if print_diffs:
                    for i in range(npart):
                        if (
                            ref.gas.RTCallsIactGradient[i]
                            != compare.gas.RTCallsIactGradient[i]
                        ):
                            print(
                                "-----",
                                ref.gas.IDs[i],
                                ref.gas.RTCallsIactGradient[i],
                                compare.gas.RTCallsIactGradient[i],
                                sml_outside_tolerance[i],
                                sml_diff[i],
                            )

                if break_on_diff:
                    quit()

        # Transport Loop Interaction Calls
        fishy = (
            ref.gas.RTCallsIactTransportInteraction
            != compare.gas.RTCallsIactTransportInteraction
        )
        if fishy.any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            account_for_sml_diff = np.count_nonzero(
                np.logical_and(fishy, sml_outside_tolerance)
            )
            print(
                "--- Calls to iact transport: count differ: {0:8d} / {1:8d}; ".format(
                    np.count_nonzero(fishy), npart
                ),
                "After removing ones with acceptable h differences: {0:8d}".format(
                    account_for_sml_diff
                ),
            )

            if account_for_sml_diff > 0:
                if print_diffs:
                    for i in range(npart):
                        if (
                            ref.gas.RTCallsIactTransport[i]
                            != compare.gas.RTCallsIactTransport[i]
                        ):
                            print(
                                "-----",
                                ref.gas.IDs[i],
                                ref.gas.RTCallsIactTransport[i],
                                compare.gas.RTCallsIactTransport[i],
                                sml_outside_tolerance[i],
                                sml_diff[i],
                            )

                if break_on_diff:
                    quit()

        # --------------------------------------------------------------
        # Check that every particle called in a previous task is also
        # called in later tasks during one time step.
        # In the code, I check that when a later task is called, the
        # previous one has been executed already. Now check that if
        # it's been called before, it will be called later.
        # --------------------------------------------------------------
        nzs = compare.gas.InjectionDone > 0
        if (compare.gas.GradientsDone[nzs] == 0).any():
            print("Oh no 1")
        if (compare.gas.TransportDone[nzs] == 0).any():
            print("Oh no 2")
        if (compare.gas.ThermochemistryDone[nzs] == 0).any():
            print("Oh no 3")

        # ---------------------------------------------------------------
        # Check numbers of subcycles.
        # ---------------------------------------------------------------
        fishy = ref.gas.nsubcycles != compare.gas.nsubcycles
        if fishy.any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            print(
                "--- Subcycle Calls count differ: {0:8d} / {1:8d}; ".format(
                    np.count_nonzero(fishy), npart
                )
            )
            if not skip_last_snap:
                print(
                    "Note, this might be acceptable behaviour for the final snapshot. You currently aren't skipping it in this check."
                )

            if break_on_diff:
                quit()

    return


def check_all_stars_is_equal(snapdata):
    """
    Check that all the star quantities are equal in every snapshot
    (for the relevant quantities of course.)
    """

    ref = snapdata[0]
    npart = ref.stars.coords.shape[0]

    print("checking stars")

    for compare in snapdata[1:]:

        if not compare.has_stars:
            continue

        # Coordinates
        if not skip_coords:

            diff = np.abs((ref.stars.coords - compare.stars.coords) / ref.stars.coords)
            if (diff > float_comparison_tolerance).any():
                print("- Comparing stars", ref.snapnr, "->", compare.snapnr)
                print("--- Coordinates vary")
                if print_diffs:
                    for i in range(ref.stars.coords.shape[0]):
                        if (
                            (ref.stars.coords[i] - compare.stars.coords[i])
                            / ref.stars.coords[i]
                        ).any():
                            print(ref.stars.coords[i], "|", compare.stars.coords[i])

                if break_on_diff:
                    quit()

        # Smoothing Lengths
        if not skip_sml:

            diff = np.abs((ref.stars.h - compare.stars.h) / ref.stars.h)
            if (diff > float_comparison_tolerance).any():
                print("- Comparing stars", ref.snapnr, "->", compare.snapnr)
                print("--- Smoothing Lengths vary")
                if print_diffs:
                    for i in range(npart):
                        if (
                            (ref.stars.h[i] - compare.stars.h[i]) / ref.stars.h[i]
                        ).any():
                            print(ref.stars.h[i], "|", compare.stars.h[i])

                if break_on_diff:
                    quit()

        # Check all emission rates are set everywhere
        fishy = ref.stars.EmissionRateSet != compare.stars.EmissionRateSet
        if fishy.any():

            print("- Comparing stars", ref.snapnr, "->", compare.snapnr)
            print("--- EmissionRateSet vary")
            if print_diffs:
                for i in range(npart):
                    if ref.stars.EmissionRateSet[i] != compare.stars.EmissionRateSet[i]:
                        print(
                            ref.stars.EmissionRateSet[i],
                            "|",
                            compare.stars.EmissionRateSet[i],
                        )

            if break_on_diff:
                quit()

        # Check all emitted radiation is equal
        fishy = ref.stars.InjectionInteractions != compare.stars.InjectionInteractions
        if fishy.any():

            print("- Comparing stars", ref.snapnr, "->", compare.snapnr)
            print("--- InjectionInteractions vary")
            if print_diffs:
                for i in range(npart):
                    if (
                        ref.stars.InjectionInteractions[i]
                        != compare.stars.InjectionInteractions[i]
                    ):
                        print(
                            ref.stars.InjectionInteractions[i],
                            "|",
                            compare.stars.InjectionInteractions[i],
                        )

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

    check_all_hydro_is_equal(snapdata)
    check_all_stars_is_equal(snapdata)


if __name__ == "__main__":
    main()
