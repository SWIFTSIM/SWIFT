################################################################################
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
################################################################################
"""
Compare the simulated growth of the HII region radius against the classical
analytic Stromgren-sphere solution (equilibrium radius + Spitzer 1978 D-type
expansion).

This is a physics-correctness check, not just a "no crash" check: a
task-graph bug that silently truncates the ionization search radius still
produces a smooth, monotonic, plausible-looking r_hii(t) curve -- it would
just be quantitatively wrong. Comparing against the analytic solution catches
that a raw "did it crash" check cannot.

The ionizing photon rate Q_H(mass) is reimplemented here from
src/feedback/GEAR/radiation.c's
radiation_get_individual_star_ionizing_photon_emission_rate_fit() (and the
radius/luminosity fits it calls), so this script isn't guessing at a
different photon budget than the code itself used. Keep this in sync with
that file if the fit ever changes.

Usage:
    python3 stromgren_analytic_check.py [-s snap/snapshot] [-o out.png]
"""
import argparse
import glob
import os

import numpy as np
from astropy import constants as const
from astropy import units as u

try:
    import swiftsimio as sw
except ImportError:
    sw = None

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Q_H(mass): direct port of radiation_get_individual_star_ionizing_photon_
# emission_rate_fit() and the radius/luminosity/temperature fits it calls
# (src/feedback/GEAR/radiation.c). All in CGS/astropy units here; the C code
# works in SWIFT internal units but the physics is unit-system independent.
# -----------------------------------------------------------------------------
def star_radius(mass_msun):
    """Empirical mass-radius relation (radiation_get_individual_star_radius)."""
    if mass_msun < 1.0:
        return (mass_msun**0.8) * const.R_sun
    elif mass_msun < 8.0:
        return (mass_msun**0.57) * const.R_sun
    else:
        return (mass_msun**0.5) * const.R_sun


def star_luminosity(mass_msun):
    """Empirical mass-luminosity relation (radiation_get_individual_star_luminosity)."""
    if mass_msun < 0.43:
        lum_sol = 0.185 * mass_msun**2
    elif mass_msun < 2.0:
        lum_sol = mass_msun**4
    elif mass_msun < 54.0:
        lum_sol = 1.5 * mass_msun**3.5
    else:
        lum_sol = 32000.0 * mass_msun
    return lum_sol * const.L_sun


def ionizing_photon_rate(mass_msun):
    """
    Q_H, the ionizing photon emission rate for a single star of the given
    mass, from a blackbody fit to its (R, L) -- direct port of
    radiation_get_individual_star_ionizing_photon_emission_rate_fit().
    """
    R = star_radius(mass_msun)
    L = star_luminosity(mass_msun)
    if R <= 0 or L <= 0:
        return 0.0 / u.s

    R_in_Rsun = (R / const.R_sun).decompose().value
    L_in_Lsun = (L / const.L_sun).decompose().value

    T = 5780.0 * (L_in_Lsun / R_in_Rsun**2) ** 0.25 * u.K

    E_threshold = 13.605 * u.eV
    x_0 = (E_threshold / (const.k_B * T)).decompose().value

    if x_0 > 45.0:
        return 0.0 / u.s

    photon_integral_sum = 0.0
    for n in range(1, 6):
        exp_term = np.exp(-n * x_0)
        if exp_term < 1e-10:
            break
        photon_integral_sum += (exp_term / n) * (
            x_0**2 + (2.0 * x_0) / n + 2.0 / n**2
        )

    prefactor = 15.0 / np.pi**4
    N_dot_ion = (L / (const.k_B * T)) * prefactor * photon_integral_sum
    return N_dot_ion.to(1 / u.s)


# -----------------------------------------------------------------------------
# Analytic Stromgren-sphere solution
# -----------------------------------------------------------------------------
# Case-B recombination coefficient at T ~ 1e4 K (src/physical_constants_cgs.h,
# const_caseb_recomb_cgs).
ALPHA_B = 2.6e-13 * u.cm**3 / u.s

# Isothermal sound speed of photoionized (T ~ 1e4 K) hydrogen gas, the
# standard value used in the classical D-type expansion solution.
C_S_IONIZED = 12.85 * u.km / u.s


def stromgren_radius(Q_H, n_H):
    """Equilibrium (Stromgren) radius: R_st = (3 Q_H / (4 pi alpha_B n_H^2))^(1/3)."""
    R_st3 = 3.0 * Q_H / (4.0 * np.pi * ALPHA_B * n_H**2)
    return R_st3.to(u.cm**3) ** (1.0 / 3.0)


def dtype_expansion_radius(t, R_st, c_s=C_S_IONIZED):
    """Spitzer (1978) D-type expansion: R(t) = R_st (1 + 7 c_s t / (4 R_st))^(4/7)."""
    return R_st * (1.0 + 7.0 * c_s * t / (4.0 * R_st)) ** (4.0 / 7.0)


# -----------------------------------------------------------------------------
# Read the simulated r_hii(t) from the snapshots
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
    if sw is None:
        raise RuntimeError("swiftsimio is required to read the snapshots.")

    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    times = []
    r_hii = []
    star_mass_msun = None
    n_H_atom_cc = None
    boxsize_kpc = None

    for f in files:
        data = sw.load(f)
        if boxsize_kpc is None:
            # Smallest box dimension -- the periodic image distance that
            # actually bounds how large r_hii can grow before the HII region
            # starts interacting with its own periodic copies.
            boxsize_kpc = float(np.min(data.metadata.boxsize.to("kpc").value))
        if len(data.stars.hiiregion_radii) == 0:
            continue

        # swiftsimio/unyt quantities can't be combined directly with
        # astropy.units quantities -- convert to plain floats in a fixed
        # unit at this boundary, then re-wrap as astropy Quantities.
        times.append(data.metadata.time.to("Myr").value)
        r_hii.append(float(np.max(data.stars.hiiregion_radii).to("kpc").value))

        if star_mass_msun is None:
            star_mass_msun = float(np.max(data.stars.masses).to("Msun").value)
        if n_H_atom_cc is None:
            # Hydrogen number density from the gas density and the
            # Grackle-reported hydrogen mass fraction (see output.log:
            # "grackle_chemistry_data.HydrogenFractionByMass"). Fall back to
            # the primordial default if it can't be inferred from the data.
            rho_g_cm3 = float(np.mean(data.gas.densities).to("g/cm**3").value)
            X_H = 0.716
            n_H_atom_cc = (
                rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p
            ).to(1 / u.cm**3)

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        star_mass_msun,
        n_H_atom_cc,
        boxsize_kpc * u.kpc,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot-glob",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for the snapshots to read (default: snap/snapshot_*.hdf5)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="stromgren_analytic_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.3,
        help="Relative-error tolerance at the final time for a pass (default: 0.3).",
    )
    args = parser.parse_args()

    t_sim, r_sim, star_mass_msun, n_H, boxsize = read_simulated_r_hii(
        args.snapshot_glob
    )

    if star_mass_msun is None or n_H is None:
        raise RuntimeError("Could not infer star mass / n_H from the snapshots.")

    Q_H = ionizing_photon_rate(star_mass_msun)
    R_st = stromgren_radius(Q_H, n_H)
    box_half_width = 0.5 * boxsize

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print(f"n_H                : {n_H:.4g}")
    print(f"Q_H                : {Q_H:.4g}")
    print(f"Equilibrium R_st   : {R_st.to(u.pc):.4g} = {R_st.to(u.kpc):.4g}")
    print(f"Box half-width     : {box_half_width.to(u.pc):.4g}")
    if R_st > box_half_width:
        print(
            "  WARNING: R_st exceeds the box half-width -- the periodic box "
            "is too small to contain a full infinite-medium Stromgren "
            "sphere. The analytic curve below is an upper bound, not a "
            "value this finite, periodic setup can actually reach; a "
            "mismatch alone does not indicate a code bug."
        )

    r_analytic = dtype_expansion_radius(t_sim, R_st).to(u.kpc)

    # Relative error at the final simulated time.
    rel_error = float(
        (np.abs(r_sim[-1] - r_analytic[-1]) / r_analytic[-1]).decompose()
    )
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"Final time: t={t_sim[-1]:.4g}  r_sim={r_sim[-1]:.4g}  "
        f"r_analytic={r_analytic[-1]:.4g}  rel_error={rel_error:.2%}  [{verdict}]"
    )

    fig, ax = plt.subplots()
    ax.plot(t_sim.to(u.Myr), r_sim.to(u.pc), "o-", label="Simulation ($r_{hii}$)")
    ax.plot(
        t_sim.to(u.Myr),
        r_analytic.to(u.pc),
        "--",
        label="Analytic (Spitzer 1978 D-type)",
    )
    ax.axhline(
        R_st.to(u.pc).value, color="grey", linestyle=":", label="Equilibrium $R_{st}$"
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("HII region radius [pc]")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
