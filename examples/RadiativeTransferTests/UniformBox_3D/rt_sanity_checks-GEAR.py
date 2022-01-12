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
# Collection of sanity checks for the 'GEAR' RT scheme in swift for the
# Any run done with the 'GEAR' scheme must pass these tests. Swift must
# be compiled with the '--enable-debugging-checks' flag.
#
# Usage:
#   ./rt_sanity_checks-GEAR.py
# or
#   ./rt_sanity_checks-GEAR.py snapshot_basename
#
# where 'snapshot_basename' is the base name of the snapshots that you
# want to work with. (Same argument as Snapshots:basename in the .yml
# file). If 'snapshot_basename' is not given, this script assumes it to
# be 'output'.
#
# The script will then read all the snapshot_basename_XXXX.hdf5 files
# in the directory where you're running the script from
# -----------------------------------------------------------------------


import unyt
import numpy as np
from sys import argv
from swift_rt_GEAR_io import get_snap_data
from matplotlib import pyplot as plt

# some behaviour options
print_diffs = True  # print differences you find
print_additional_information = True
break_on_diff = False  # quit when you find a difference
skip_plots = False  # skip showing plots for diagnosis

# tolerance for a float to be equal
float_comparison_tolerance = 1e-4
# tolerance for a float that was summed up over all particles to vary
float_particle_sum_comparison_tolerance = 5e-4
# tolerance for meshless energy distribution scheme during injeciton comparison
float_psi_comparison_tolerance = 5e-4

# -------------------------------------------------------------------------------


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

    # ----------------------------------------------------------------
    # Check 0: Make sure we don't have NaNs or Infs anywhere
    # ----------------------------------------------------------------
    for snap in snapdata:

        photon_energies = snap.gas.PhotonEnergies
        ok = np.isfinite(photon_energies)
        if not ok.all():
            print("In snapshot", snap.snapnr, ":")
            print("Found NaNs/infs in photon energies:", np.count_nonzero(ok == 0))
            if break_on_diff:
                quit()

        photon_fluxes = snap.gas.PhotonFluxes
        ok = np.isfinite(photon_fluxes)
        if not ok.any():
            print("In snapshot", snap.snapnr, ":")
            print("Found NaNs/infs in photon fluxes:", np.count_nonzero(ok == 0))
            if break_on_diff:
                quit()

        if snap.has_stars:
            injected_energies = snap.stars.InjectedPhotonEnergy
            ok = np.isfinite(injected_energies)
            if not ok.all():
                print("In snapshot", snap.snapnr, ":")
                print(
                    "Found NaNs/infs in star injected energies:",
                    np.count_nonzero(ok == 0),
                )
                if break_on_diff:
                    quit()

    # ----------------------------------------------------------------
    # Check 1: Make sure the right amount of energy has been injected
    # into the gas
    # ----------------------------------------------------------------

    if not rundata.has_stars:
        print("Found no stars in run. Skipping injection tests")
        return

    emission_rates = rundata.const_emission_rates
    ngroups = rundata.ngroups

    initial_energies = np.sum(snapdata[0].gas.PhotonEnergies, axis=0)
    initial_time = snapdata[0].time

    # Check 1a) : sum initial energy + sum injected energy = sum current energy
    # --------------------------------------------------------------------------
    # TODO: this assumes no cosmological expansion

    for snap in snapdata[1:]:
        dt = snap.time - initial_time
        # sum of each group over all particles
        photon_energies = np.sum(snap.gas.PhotonEnergies, axis=0)
        if snap.has_stars:
            # in case we only have 1 star, the sum returns a scalar, so add atleast_1d
            injected_energies = np.atleast_1d(
                np.sum(snap.stars.InjectedPhotonEnergy, axis=0)
            )
        else:
            injected_energies = np.zeros(ngroups, dtype=np.float)

        for g in range(ngroups):
            energy_expected = initial_energies[g] + injected_energies[g]
            if energy_expected != 0.0:
                diff = abs(1.0 - photon_energies[g] / energy_expected)
                if diff > float_particle_sum_comparison_tolerance:
                    print(
                        "Injection Energy Budget is wrong; "
                        + "snapshot {0:d} tolerance {1:.2e}".format(
                            snap.snapnr, float_particle_sum_comparison_tolerance
                        )
                    )
                    if print_diffs:
                        print("--- group:", g + 1)
                        print("----- diff:            ", diff)
                        print("----- photon energies: ", photon_energies[g])
                        print("----- expected:        ", energy_expected)
                    if break_on_diff:
                        quit()

    # Check 1b: compute injection manually and compare
    # --------------------------------------------------------------------------
    #  Stars don't inject energy until the first time when they are active.
    #  If they keep their time step size constant, then once they're active,
    #  they will inject the exact amount of energy that should've been injected
    #  integrated since the beginning of time. From this point on until they're
    #  active again, the total injected energy will be below what's expected.
    #  Make sure that we never exceed the analytically prescribed.
    #  Remember: The reason we have too little injected energy is because we
    #  don't inject any energy during the zeroth time step. We can't, since the
    #  zeroth time step is the one that determines the time step size of the star.

    # TODO: this assumes a constant number of stars. You need to deal with SF
    # TODO: this assumes no cosmological expansion

    upper_boundary_for_plot = []
    snaps_for_1bplot = []

    initial_time = snapdata[0].time
    if snapdata[0].has_stars:
        emission_at_initial_time = snapdata[0].stars.InjectedPhotonEnergy.sum(axis=0)
    else:
        emission_at_initial_time = np.zeros(rundata.ngroups, dtype=np.float) * unyt.erg

    if rundata.use_const_emission_rate and not rundata.hydro_controlled_injection:
        if len(snapdata) <= 2:
            # because it's useless to check only snap_0000
            print("Check 1b: You need at least 2 snapshots to do this particular test")
        else:
            diffs_for_plot = []
            energies_for_plot = []
            found_potential_error = False
            for snap in snapdata[1:]:  # skip snapshot zero
                dt = snap.time - initial_time
                if snap.has_stars:
                    injected_energies = np.atleast_1d(
                        snap.stars.InjectedPhotonEnergy.sum(axis=0)
                        - emission_at_initial_time
                    )
                else:
                    injected_energies = np.zeros(ngroups) * unyt.erg
                energies_expected = snap.nstars * emission_rates * dt
                energies_expected = energies_expected.to(injected_energies.units)
                diff = np.array(injected_energies / energies_expected - 1.0)

                upper_boundary_for_plot.append(energies_expected)
                energies_for_plot.append(injected_energies)
                diffs_for_plot.append(diff)
                snaps_for_1bplot.append(snap.snapnr)

                # diff should be < 0. Allow for some tolerance here
                if (diff > float_psi_comparison_tolerance).any():
                    print(
                        "Injection Energy Prediction upper boundary is wrong; "
                        + "snapshot {0:d} tolerance {1:.2e}".format(
                            snap.snapnr, float_psi_comparison_tolerance
                        )
                    )
                    for g in range(ngroups):
                        #  if energies_expected[g] > injected_energies[g]:
                        print("--- group", g + 1)
                        print("----- injected:", injected_energies[g])
                        print(
                            "----- expected:", energies_expected[g], "should be smaller"
                        )
                        print(
                            "----- ratio   :",
                            (injected_energies[g] / energies_expected[g]),
                        )
                        print("----- diff    :", diff[g], "should be < 0")
                        found_potential_error = True

                        if break_on_diff:
                            quit()

            if not skip_plots and found_potential_error:
                # Make this plot if there are possible errors
                diffs_for_plot = np.array(diffs_for_plot)
                plt.figure()
                for g in range(ngroups):
                    plt.plot(
                        snaps_for_1bplot,
                        diffs_for_plot[:, g],
                        label="group {0:d}".format(g + 1),
                    )
                    plt.plot(
                        [snaps_for_1bplot[0], snaps_for_1bplot[-1]],
                        [0, 0],
                        "k",
                        label="upper boundary",
                    )
                plt.legend()
                plt.xlabel("snapshot")
                plt.ylabel("injected energy / expected energy - 1")
                plt.title(
                    "Difference from expected injected energy - something's fishy"
                )
                plt.show()
                plt.close()

    # --------------------------------
    # Create additional plots?
    # --------------------------------

    if rundata.use_const_emission_rate and not rundata.hydro_controlled_injection:
        if not skip_plots and len(energies_for_plot) > 2:
            # Show me the plot that the injected energy
            # is correctly bounded
            upper_boundary_for_plot = np.array(upper_boundary_for_plot)
            energies_for_plot = np.array(energies_for_plot)

            plt.figure(figsize=(6, 4))
            for g in range(ngroups):
                plt.plot(
                    snaps_for_1bplot,
                    energies_for_plot[:, g],
                    label="group {0:d}".format(g + 1),
                    color="C" + str(g),
                )
                plt.plot(
                    snaps_for_1bplot,
                    upper_boundary_for_plot[:, g],
                    ":",
                    label="group {0:d} upper boundary".format(g + 1),
                    color="C" + str(g),
                )
            plt.legend(fontsize=10)
            plt.xlabel("snapshot")
            plt.ylabel(
                "radiation energy [$"
                + energy_expected.units.latex_representation()
                + "$]"
            )
            plt.title("Injected energies vs expectations")
            plt.show()
            plt.close()

    return


def main():
    """
    Main function to run.
    """
    snapdata, rundata = get_snap_data(prefix=file_prefix)

    check_injection(snapdata, rundata)

    return


if __name__ == "__main__":
    main()
