#!/usr/bin/env python3

# -----------------------------------------------------------------------
# Collection of checks for the 'debug' RT scheme in swift for the
# uniform box test where particles don't move and every time step an
# output file is generated.
#
# Usage:
#   ./rt_checks_uniform.py
# or
#   ./rt_checks_uniform.py snapshot_basename
#
# where 'snapshot_basename' is the base name of the snapshots that you
# want to work with. (Same argument as Snapshots:basename in the .yml
# file). If 'snapshot_basename' is not given, this script assumes it to
# be 'output'.
#
# The script will then read all the snapshot_basename_XXXX.hdf5 files
# in the directory where you're running the script from
# -----------------------------------------------------------------------


import numpy as np
from sys import argv
from swift_rt_debug_io import get_snap_data


# some behaviour options
skip_snap_zero = True  # skip snap_0000.hdf5
skip_last_snap = True  # skip snap_0max.hdf5
skip_coords = True  # skip coordinates check
skip_sml = True  # skip smoothing lengths check
print_diffs = True  # print differences you find
break_on_diff = True  # quit when you find a difference

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
        sml_within_tolerance = sml_diff <= float_comparison_tolerance

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

        # Calls to star interactions
        if (ref.gas.RTStarIact != compare.gas.RTStarIact).any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            print("--- Calls to star interactions vary")

            if print_diffs:
                for i in range(npart):
                    if ref.gas.RTStarIact[i] != compare.gas.RTStarIact[i]:
                        print(
                            "-----",
                            ref.gas.IDs[i],
                            ref.gas.RTStarIact[i],
                            compare.gas.RTStarIact[i],
                        )

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

        # Gradient Loop Calls
        fishy = ref.gas.RTCallsIactGradient != compare.gas.RTCallsIactGradient
        if fishy.any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            account_for_sml_diff = np.count_nonzero(
                np.logical_and(fishy, sml_within_tolerance)
            )
            print(
                "--- Calls to iact gradient: count differ: {0:8d} / {1:8d}; After removing ones with acceptable h differences: {2:8d}".format(
                    np.count_nonzero(fishy), npart, account_for_sml_diff
                )
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
                                sml_within_tolerance[i],
                                sml_diff[i],
                            )

                if break_on_diff:
                    quit()

        # Transport Loop Calls
        fishy = ref.gas.RTCallsIactTransport != compare.gas.RTCallsIactTransport
        if fishy.any():
            print("- Comparing hydro", ref.snapnr, "->", compare.snapnr)
            account_for_sml_diff = np.count_nonzero(
                np.logical_and(fishy, sml_within_tolerance)
            )
            print(
                "--- Calls to iact transport: count differ: {0:8d} / {1:8d}; After removing ones with acceptable h differences: {2:8d}".format(
                    np.count_nonzero(fishy), npart, account_for_sml_diff
                )
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
                                sml_within_tolerance[i],
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

            diff = np.abs((ref.gas.h - compare.gas.h) / ref.gas.h)
            if (diff > float_comparison_tolerance).any():
                print("- Comparing stars", ref.snapnr, "->", compare.snapnr)
                print("--- Smoothing Lengths vary")
                if print_diffs:
                    for i in range(npart):
                        if ((ref.gas.h[i] - compare.gas.h[i]) / ref.gas.h[i]).any():
                            print(ref.gas.h[i], "|", compare.gas.h[i])

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
