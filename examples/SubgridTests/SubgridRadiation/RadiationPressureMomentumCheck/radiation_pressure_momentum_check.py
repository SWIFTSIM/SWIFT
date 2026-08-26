###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
"""Tier-1 radiation-pressure momentum-conservation self-check (see README).

Independently recomputes, in Python, the same momentum-injection rate the C
code derives (radiation_get_star_physical_radiation_pressure, src/feedback/
GEAR/radiation.c): p_dot_rad = L_bol / c * (1 + tau_IR). Integrating this over
the run and comparing against the gas's actual total momentum gain checks
that the implementation applies exactly what that formula says -- no missing
factor, no double-counting, no silently-dropped kick. It does NOT validate
that the formula itself is the right physical model (see README/logbook for
that separate question).

Two approximations, both a direct consequence of this test's uniform-box, no
photoionization/cooling/gravity design (see README): the star's local gas
density is taken as the box's mean density (uniform IC, so the SPH-local
value at the star should match it closely), and the Sobolev length term in
Sigma_gas is dropped (zero density gradient by construction). L_bol and
tau_IR are evaluated once from the initial state and assumed constant over
the run -- valid here since stellar winds are off (mass fixed) and the run is
far too short for the gas density near the star to evolve appreciably.
"""

import argparse
import glob

import numpy as np
import swiftsimio as sw
import unyt as u


def blackbody_bolometric_luminosity(mass_msun):
    """Standalone piecewise mass-luminosity fit, kept as a fixed historical
    reference for this example only. It is not synced with the current C
    code: src/feedback/GEAR/radiation.c now reads this quantity from a
    pychem-generated table instead of computing it inline.
    """
    M = float(mass_msun)
    if M < 0.43:
        lum_sol = 0.185 * M**2
    elif M < 2.0:
        lum_sol = M**4
    elif M < 54.0:
        lum_sol = 1.5 * M**3 * np.sqrt(M)
    else:
        lum_sol = 32000.0 * M
    return lum_sol * u.Lsun


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshots",
        type=str,
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots.",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=0.15,
        help="Relative-error pass/fail threshold (default: 0.15).",
    )
    args = parser.parse_args()

    snapshots = sorted(glob.glob(args.snapshots))
    if len(snapshots) < 2:
        raise RuntimeError(
            f"Need at least 2 snapshots, found {len(snapshots)} matching "
            f"{args.snapshots!r}."
        )

    data0 = sw.load(snapshots[0])
    data1 = sw.load(snapshots[-1])

    c = u.physical_constants.speed_of_light

    # --- Star properties at t=0 ---
    M_star = data0.stars.masses[0].to(u.Msun)
    h_star = data0.stars.smoothing_lengths[0].to(u.cm)

    L_bol = blackbody_bolometric_luminosity(M_star.value).to(u.erg / u.s)

    # --- Gas properties at t=0 (uniform box: mean == local-at-star) ---
    box = data0.metadata.boxsize.to(u.cm)
    volume = box[0] * box[1] * box[2]
    # float32 gas masses carry an internal-mass-unit-to-g conversion factor
    # (~2e43) that overflows float32 on its own, well before the physical
    # value does; go via Msun first, whose conversion factor is small.
    M_gas_tot = np.sum(data0.gas.masses).to(u.Msun).to(u.g)
    rho_gas = (M_gas_tot / volume).to(u.g / u.cm**3)

    # GEAR stores the TOTAL metal mass fraction as a synthetic extra "metals"
    # named column (chemistry_get_total_metal_mass_fraction_for_feedback,
    # src/chemistry/GEAR/chemistry.h reads the last of GEAR_CHEMISTRY_ELEMENT_COUNT
    # slots) -- not a mean over the individually-tracked elements (Fe, Mg, O,
    # ...), which do not sum to Z.
    Z_gas = float(np.mean(data0.gas.smoothed_metal_mass_fractions.metals))
    Z_sun = 0.02
    kappa_IR = 10.0 * (u.cm**2 / u.g) * (Z_gas / Z_sun)

    kernel_gamma_wendlandC2_3D = 1.936492  # src/kernel_hydro.h, 3D branch
    Sigma_gas = rho_gas * (h_star * kernel_gamma_wendlandC2_3D)  # Sobolev term ~0
    tau_IR = (kappa_IR * Sigma_gas).to(u.dimensionless).value

    p_dot_rad = (L_bol / c * (1.0 + tau_IR)).to(u.g * u.cm / u.s**2)

    t0 = data0.metadata.time.to(u.s)
    t1 = data1.metadata.time.to(u.s)
    dt = t1 - t0

    p_theory = (p_dot_rad * dt).to(u.Msun * u.km / u.s)

    # --- Actual momentum gained by the gas ---
    # Particles start at rest (Velocities: 0 in the IC), so |m_j v_j| at the
    # final snapshot IS each particle's own momentum kick magnitude. Summing
    # magnitudes (not the vector sum, which cancels by the star's radial
    # symmetry) is the correct comparison: per-particle weights sum to 1
    # (radiation_iact_nonsym_feedback_apply), so Sum_j |Delta p_j| ==
    # p_dot_rad * dt exactly, up to hydro back-reaction over this short
    # window and the box's Poisson density noise.
    v_gas = data1.gas.velocities.to(u.km / u.s)
    m_gas = data1.gas.masses.to(u.Msun)
    speeds = np.linalg.norm(v_gas.value, axis=1) * v_gas.units
    p_actual = np.sum(m_gas * speeds).to(u.Msun * u.km / u.s)

    rel_error = float(abs(p_actual - p_theory) / p_theory)

    print(f"Star mass                    : {M_star:.4f}")
    print(f"L_bol                        : {L_bol:.6e}")
    print(f"Mean gas density             : {rho_gas:.6e}")
    print(f"Star smoothing length        : {h_star:.6e}")
    print(f"Gas metallicity (mass frac.) : {Z_gas:.4f} (Z/Zsun = {Z_gas / Z_sun:.4f})")
    print(f"kappa_IR                     : {kappa_IR:.4f}")
    print(f"tau_IR                       : {tau_IR:.6e}")
    print(f"Elapsed time (t1 - t0)       : {dt:.6e}")
    print(f"Theoretical injected |p|     : {p_theory:.6e}")
    print(f"Actual   Sum_j |m_j v_j|     : {p_actual:.6e}")
    print(f"Relative error               : {rel_error:.2%}")

    if rel_error <= args.tolerance:
        print(f"PASS (tolerance {args.tolerance:.0%})")
    else:
        print(f"FAIL (tolerance {args.tolerance:.0%})")


if __name__ == "__main__":
    main()
