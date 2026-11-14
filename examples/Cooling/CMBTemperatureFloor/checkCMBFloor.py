#!/usr/bin/env python3
"""Check that the gas of a run holds the CMB temperature floor, in every snapshot after the first.

The floor is T_CMB,0 (1 + z) at the redshift z of the snapshot. Official Grackle uses
T_CMB,0 = 2.73 K and geajay 2.725 K: pass --t-cmb-0 for the library of the run.

The snapshots hold no temperature, so it is rebuilt as the cooling library computes it:
T = (gamma - 1) u m_H / (k_B n), with u the physical specific internal energy and n the number of
particles per m_H and per unit mass: HI + HII + e + (HeI + HeII + HeIII) / 4, plus HM and
(H2I + H2II) / 2 when the network has H2, plus the metal mass fraction / 16. Gamma is 5/3. The
H2 correction of the adiabatic index is exact below an H2 to other-particle ratio of 1e-3; above
it the check stops with an error. m_H and k_B are the values of the libraries.

Checks, in this order. Finiteness: every temperature of every snapshot is finite. Floor held: in
every snapshot the lowest temperature is at least the floor times (1 - TOLERANCE - dt_max / 2).
Here dt_max is the largest time-step of the run in ln a (TimeIntegration:dt_max of
used_parameters.yml). The snapshot energy is drifted to the snapshot time and carries up to half
a step of expansion cooling that the floor did not see, so the lag is at most dt_max / 2.
Floor reached: in the last snapshot the highest temperature is at most the floor times
(1 + REACH), so the gas did cool onto the floor. The script exits with status 1 if a check
fails. A run without the floor must fail the floor-held check.

With --expect-below-floor the script checks a run made without the floor. It exits with status 0
only if the floor-held check fails and every temperature of the last snapshot lies below half of
the floor, that is, only if the run is a control that can tell a floor from no floor.
"""

import argparse
import glob
import sys

import h5py
import numpy as np
import yaml

# Float32 round-off of the energy plus the lag of the gas behind the falling floor, relative.
# Half the largest time-step in ln a is added from the run, for the drift of the snapshot energy.
TOLERANCE = 1e-4
# Distance above the floor where the cooled gas must sit, relative.
REACH = 1e-3
# Largest H2 to other-particle ratio where the library adiabatic index is the simple formula.
H2_RATIO_MAX = 1e-3

GAMMA = 5.0 / 3.0
M_H = 1.67262171e-24  # g, as in the libraries
K_B = 1.3806504e-16  # erg/K, as in the libraries
MU_METAL = 16.0


def temperature_from_fields(u_cgs, fractions, metal_fraction):
    """Return the temperature (K) of the gas the way the cooling libraries compute it.

    Parameters
    ----------
    u_cgs : numpy.ndarray
        Physical specific internal energy in cgs.
    fractions : dict
        Mass fractions of the species, keys among HI, HII, e, HeI, HeII, HeIII, HM, H2I, H2II.
    metal_fraction : numpy.ndarray
        Total metal mass fraction used by the cooling.
    """
    helium = fractions["HeI"] + fractions["HeII"] + fractions["HeIII"]
    other = fractions["HI"] + fractions["HII"] + fractions["e"] + 0.25 * helium
    n_h2 = np.zeros_like(other)
    if "H2I" in fractions:
        other = other + fractions["HM"]
        n_h2 = 0.5 * (fractions["H2I"] + fractions["H2II"])
        if np.any(n_h2 > H2_RATIO_MAX * other):
            sys.exit(
                "H2 fraction above the range where the temperature rebuild is exact."
            )
    # Adiabatic index with H2 rotational degrees of freedom (n_h2 = 0 gives GAMMA).
    gamma_eff = 1.0 + (n_h2 + other) / (2.5 * n_h2 + other / (GAMMA - 1.0))
    number = other + n_h2 + metal_fraction / MU_METAL
    return (gamma_eff - 1.0) * u_cgs * M_H / (K_B * number)


def read_snapshot(path):
    """Return the temperatures (K) and the redshift of one snapshot."""
    with h5py.File(path, "r") as f:
        gas = f["PartType0"]
        units = f["Units"].attrs
        to_cgs = (
            float(np.ravel(units["Unit length in cgs (U_L)"])[0]) ** 2
            / float(np.ravel(units["Unit time in cgs (U_t)"])[0]) ** 2
        )
        a = float(np.ravel(f["Cosmology"].attrs["Scale-factor"])[0])
        z = float(np.ravel(f["Cosmology"].attrs["Redshift"])[0])
        energy = gas["InternalEnergies"]
        exponent = float(np.ravel(energy.attrs["a-scale exponent"])[0])
        u_cgs = energy[:].astype(np.float64) * a**exponent * to_cgs
        names = ["HI", "HII", "e", "HeI", "HeII", "HeIII"]
        if "H2I" in gas:
            names += ["HM", "H2I", "H2II"]
        fractions = {n: gas[n][:].astype(np.float64) for n in names}
        metals = gas["SmoothedMetalMassFractions"][:].astype(np.float64)[:, -1]
    return temperature_from_fields(u_cgs, fractions, metals), z


def evaluate(snapshots, t_cmb_0, dt_max=0.0):
    """Return (rows, last_ratio_max) from a list of (temperature, z) pairs."""
    for temperature, _ in snapshots:
        if not np.all(np.isfinite(temperature)):
            count = np.sum(~np.isfinite(temperature))
            return [("finite", False, f"{count} non-finite values")], None
    rows = [("finite", True, f"{len(snapshots)} snapshots")]
    bar = 1.0 - TOLERANCE - 0.5 * dt_max
    ratios = [
        (temperature.min() / (t_cmb_0 * (1.0 + z)), z) for temperature, z in snapshots
    ]
    worst, z_worst = min(ratios)
    rows.append(
        (
            "floor held",
            worst >= bar,
            f"lowest T / floor = {worst:.6g} (at z = {z_worst:.4g}), bar {bar:.6g}",
        )
    )
    temperature, z = snapshots[-1]
    top = temperature.max() / (t_cmb_0 * (1.0 + z))
    rows.append(
        (
            "floor reached",
            top <= 1.0 + REACH,
            f"highest T / floor = {top:.6g} in the last snapshot, bar {1.0 + REACH:.6g}",
        )
    )
    return rows, top


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="Run directory")
    parser.add_argument(
        "--t-cmb-0", type=float, required=True, help="T_CMB,0 of the library (K)"
    )
    parser.add_argument(
        "--expect-below-floor",
        action="store_true",
        help="Check a run made without the floor (see above)",
    )
    opt = parser.parse_args()

    files = sorted(glob.glob(f"{opt.run}/snap/snapshot_*.hdf5"))
    if len(files) < 2:
        sys.exit(f"{opt.run}: need at least two snapshots.")
    snapshots = [read_snapshot(path) for path in files[1:]]
    with open(f"{opt.run}/used_parameters.yml") as stream:
        dt_max = float(yaml.safe_load(stream)["TimeIntegration"]["dt_max"])
    rows, _ = evaluate(snapshots, opt.t_cmb_0, dt_max)
    for name, passed, detail in rows:
        print(f"{name:14s} {'PASS' if passed else 'FAIL'}  {detail}")

    if not opt.expect_below_floor:
        sys.exit(0 if all(passed for _, passed, _ in rows) else 1)

    temperature, z = snapshots[-1]
    floor = opt.t_cmb_0 * (1.0 + z)
    below = bool(np.all(np.isfinite(temperature))) and temperature.max() < 0.5 * floor
    print(
        f"{'below half':14s} {'PASS' if below else 'FAIL'}  "
        f"highest T / floor = {temperature.max() / floor:.6g}, bar 0.5"
    )
    held = dict((n, p) for n, p, _ in rows).get("floor held", True)
    sys.exit(0 if (not held and below) else 1)


if __name__ == "__main__":
    main()
