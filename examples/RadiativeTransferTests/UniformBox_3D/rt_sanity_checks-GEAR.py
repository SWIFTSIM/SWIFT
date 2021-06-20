#!/usr/bin/env python3

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
print_diffs = False  # print differences you find
print_additional_information = True
break_on_diff = False  # quit when you find a difference
skip_plots = False  # skip showing plots for diagnosis

# tolerance for a float to be equal
float_comparison_tolerance = 1e-4
# tolerance for a float that was summed up over all particles to vary
float_particle_sum_comparison_tolerance = 5e-4
# tolerance for meshless energy distribution scheme during injeciton comparison
float_psi_comparison_tolerance = 5e-4


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

        injected_energies = snap.stars.InjectedPhotonEnergy
        ok = np.isfinite(injected_energies)
        if not ok.all():
            print("In snapshot", snap.snapnr, ":")
            print(
                "Found NaNs/infs in star injected energies:", np.count_nonzero(ok == 0)
            )
            if break_on_diff:
                quit()

    # ----------------------------------------------------------------
    # Check 1: Make sure the right amount of energy has been injected
    # into the gas
    # ----------------------------------------------------------------

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
        injected_energies = np.sum(snap.stars.InjectedPhotonEnergy, axis=0)

        for g in range(ngroups):
            energy_expected = initial_energies[g] + injected_energies[g]
            diff = abs(1.0 - photon_energies[g] / energy_expected)
            if diff > float_particle_sum_comparison_tolerance:
                print("Injection Energy Budget is wrong; snapshot", snap.snapnr)
                if print_diffs:
                    print("--- group:", g + 1)
                    print("----- diff:            ", diff)
                    print("----- photon energies: ", photon_energies[g])
                    print("----- expected:        ", energy_expected)
                if break_on_diff:
                    quit()

    # Check 1b) : sum injected energy >= (t_now - t_start * injection_rate)
    # --------------------------------------------------------------------------
    # we may have injected too much energy, because stars inject all the
    # radiation of their entire time step instantaneously. The assumed duration
    # of the time step might not have ended by the time of the writing of the
    # snapshot, so we end up with more energy than a smooth injection would
    # predict.

    # TODO: this assumes a constant number of stars. You need to deal with SF
    # TODO: this assumes no cosmological expansion

    # NOTE: we start with the data from snapshot 1, assuming that all stars have
    # been active at least once by the time the snapshot is written. That isn't
    # necessarily the case if we start at snapshot 0, and this test could find
    # that the energy is missing. (Essentially you'd just need to know the time
    # at which all stars have been active at least once to do this test
    # correctly.)
    initial_time = snapdata[1].time
    emission_at_initial_time = snapdata[1].stars.InjectedPhotonEnergy.sum(axis=0)

    if rundata.use_const_emission_rate and not rundata.hydro_controlled_injection:
        if len(snapdata) <= 2:
            print(
                "Check 1b: You need at least 2 snapshots to", "do this particular test"
            )
        else:
            diffs_for_plot = []
            snaps_for_plot = []
            for snap in snapdata[2:]:
                dt = snap.time - initial_time
                injected_energies = (
                    snap.stars.InjectedPhotonEnergy.sum(axis=0)
                    - emission_at_initial_time
                )
                energies_expected = snap.nstars * emission_rates * dt
                energies_expected = energies_expected.to(injected_energies.units)

                diff = np.array(1.0 - injected_energies / energies_expected)
                diffs_for_plot.append(diff)
                snaps_for_plot.append(snap.snapnr)

                if (np.abs(diff) > float_psi_comparison_tolerance).any():
                    print(
                        "Snapshot", snap.snapnr, "Injected Energy Prediction is wrong;"
                    )
                    for g in range(ngroups):
                        #  if energies_expected[g] > injected_energies[g]:
                        print("--- group", g + 1)
                        print("----- injected:", injected_energies[g])
                        print("----- expected:", energies_expected[g])
                        print(
                            "----- ratio   :",
                            (injected_energies[g] / energies_expected[g]),
                        )
                        print("----- diff    :", diff[g])
                        if diff[g] < 0:
                            print(
                                "----- overshoot; this is expected for discrete star timesteps"
                            )

                        if break_on_diff:
                            quit()

            if not skip_plots:
                diffs_for_plot = np.array(diffs_for_plot)
                plt.figure()
                for g in range(ngroups):
                    plt.plot(
                        snaps_for_plot,
                        diffs_for_plot[:, g],
                        label="group {0:d}".format(g + 1),
                    )
                plt.legend()
                plt.xlabel("snapshot")
                plt.ylabel("1 - injected energy / expected energy")
                plt.title(
                    "Difference from expected injected energy. See comments for interpretation."
                )
                # Here's the comment:
                # If not all stars are updated right before a snapshot, too much energy may have been
                # injected. That is fine, and results in a negative value on the y axis. With
                # increasing snapshot number, the total injected energy should increase, and the
                # energy overshoots should decrease in magnitude, so the curves should converge towards
                # zero. Lastly, float precision summation roundoff errors can add deviations in the order
                # of 10^-3. So you should see:
                #   - decreasing magnitude of curves
                #   - values are negative, save for roundoff errors
                # if that's what you see, you're good.
                plt.show()


def main():
    """
    Main function to run.
    """

    snapdata, rundata = get_snap_data(prefix=file_prefix)

    check_injection(snapdata, rundata)

    return


if __name__ == "__main__":
    main()
