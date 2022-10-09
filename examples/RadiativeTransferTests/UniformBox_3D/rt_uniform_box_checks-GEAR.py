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

    emission_rates = rundata.const_emission_rates
    ngroups = rundata.ngroups

    initial_energies = [
        snapdata[0].gas.PhotonEnergies[g].sum(axis=0) for g in range(ngroups)
    ]
    initial_time = snapdata[0].time

    for snap in snapdata:
        dt = snap.time - initial_time
        injected = snap.nstars * emission_rates * dt
        for g in range(ngroups):
            photon_energy = snap.gas.PhotonEnergies[g].sum(axis=0)
            energy_expected = initial_energies[g] + injected[g]
            diff = 1.0 - energy_expected / photon_energy

            if abs(diff) > float_particle_sum_comparison_tolerance:
                print("Snapshot", snap.snapnr, "Injected Energy Prediction is wrong;")
                if print_diffs:
                    print("--- group ", g + 1)
                    print("----- diff:           ", diff)
                    print("----- photon energies:", photon_energy)
                    print("----- expected:       ", energy_expected)
                if break_on_diff:
                    quit()

    return


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
