#!/usr/bin/env python3

# -----------------------------------------------------------------------
# Collection of checks for the 'GEAR' RT scheme in swift for the
# uniform box test where particles don't move and every time step an
# output file is generated. Swift must be compiled with the
# '--enable-debugging-checks' flag.
#
# Usage:
#   ./rt_uniform_box_checks-GEAR.py
# or
#   ./rt_uniform_box_checks-GEAR.py snapshot_basename
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
import unyt
from sys import argv
from swift_rt_GEAR_io import get_snap_data


# some behaviour options
print_diffs = True  # print differences you find
break_on_diff = False  # quit when you find a difference

# tolerance for a float to be equal
float_comparison_tolerance = 1e-5
# tolerance for a float that was summed up over all particles to vary
float_particle_sum_comparison_tolerance = 1e-3


if len(argv) > 1:
    file_prefix = argv[1]
else:
    file_prefix = "output"


def check_injection(snapdata, rundata):
    """
    Do checks related to energy injections.

    snapdata: list of swift_rt_GEAR_io.RTSnapData objects

    rundata: swift_rt_GEAR_io.Rundata object
    """

    print("checking injection")

    if not rundata.use_const_emission_rate:
        print("Sim wasn't run with const emission rates. Skipping those tests.")
        return

    # ----------------------------------------------------------------
    # Check 1: Make sure the right amount of energy has been injected
    # into the gas
    # ----------------------------------------------------------------

    initial_energies = snapdata[0].gas.PhotonEnergies.sum(axis=0)
    initial_time = snapdata[0].time

    emission_rates = rundata.const_emission_rates
    ngroups = rundata.ngroups

    for snap in snapdata:
        dt = snap.time - initial_time
        photon_energies = snap.gas.PhotonEnergies.sum(axis=0)
        injected = snap.nstars * emission_rates * dt
        energies_expected = initial_energies + injected
        diff = np.array(1.0 - energies_expected / photon_energies)

        if (np.abs(diff) > float_particle_sum_comparison_tolerance).any():
            print("Snapshot", snap.snapnr, "Injected Energy Prediction is wrong;")
            if print_diffs:
                for g in range(ngroups):
                    if abs(diff[g]) > float_particle_sum_comparison_tolerance:
                        print("--- group ", g + 1)
                        print("----- diff:           ", diff[g])
                        print("----- photon energies:", photon_energies[g])
                        print("----- expected:       ", energies_expected[g])
                if break_on_diff:
                    quit()


def main():
    """
    Main function to run.
    """

    snapdata, rundata = get_snap_data(
        prefix=file_prefix, skip_snap_zero=False, skip_last_snap=False
    )

    check_injection(snapdata, rundata)

    return


if __name__ == "__main__":
    main()
