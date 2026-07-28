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
Cosmological counterpart of ../StromgrenSphere/stromgren_analytic_check.py:
compares the simulated HII region growth against the classical Stromgren
equilibrium radius + Spitzer (1978) D-type expansion, evaluated at the
INSTANTANEOUS PHYSICAL density.

Every length/density SWIFT writes to a cosmological snapshot is COMOVING.
This script converts to physical explicitly, per field, using that field's
own "a-scale exponent" HDF5 attribute and the snapshot's Header/Scale-factor
(physical = comoving * a**exponent) -- it does not assume a fixed convention
per field name. That conversion is the entire point of this script: for
params_identity.yml (a~=1 throughout) it should barely matter, but for
params_highz.yml (a_begin=0.1) skipping it would be a ~1000x density error.

Usage:
    python3 cosmo_stromgren_analytic_check.py [-s snap/snapshot] [-o out.png]
"""

import argparse
import glob

import h5py
import numpy as np
from astropy import constants as const
from astropy import units as u

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def a_scale_exponent(dataset):
    """The "a-scale exponent" n such that physical = comoving * a**n for
    this HDF5 dataset (written by SWIFT on every field it outputs)."""
    return float(dataset.attrs["a-scale exponent"][0])


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
    """Q_H, the ionizing photon emission rate for a single star of the given
    mass, from a blackbody fit to its (R, L) -- direct port of
    radiation_get_individual_star_ionizing_photon_emission_rate_fit()."""
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
        photon_integral_sum += (exp_term / n) * (x_0**2 + (2.0 * x_0) / n + 2.0 / n**2)

    prefactor = 15.0 / np.pi**4
    N_dot_ion = (L / (const.k_B * T)) * prefactor * photon_integral_sum
    return N_dot_ion.to(1 / u.s)


# -----------------------------------------------------------------------------
# Analytic Stromgren-sphere solution, evaluated at the run's own applied
# T_ionized (see ../StromgrenSphere/stromgren_analytic_check.py for why this
# departs from the classical flat 1e4 K assumption).
# -----------------------------------------------------------------------------
def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c)."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


def sound_speed_ionized(T_ionized):
    """Isothermal sound speed of photoionized hydrogen gas, scaled from the
    standard 1e4 K reference value (12.85 km/s) to the actual applied
    ionized temperature (c_s ~ sqrt(T) at fixed composition)."""
    return 12.85 * u.km / u.s * np.sqrt(T_ionized / (1e4 * u.K))


def stromgren_radius(Q_H, n_H, alpha_B):
    """Equilibrium (Stromgren) radius: R_st = (3 Q_H / (4 pi alpha_B n_H^2))^(1/3)."""
    R_st3 = 3.0 * Q_H / (4.0 * np.pi * alpha_B * n_H**2)
    return R_st3.to(u.cm**3) ** (1.0 / 3.0)


def dtype_expansion_radius(t, R_st, c_s):
    """Spitzer (1978) D-type expansion: R(t) = R_st (1 + 7 c_s t / (4 R_st))^(4/7)."""
    return R_st * (1.0 + 7.0 * c_s * t / (4.0 * R_st)) ** (4.0 / 7.0)


# -----------------------------------------------------------------------------
# Read the simulated, PHYSICAL r_hii(t) from the snapshots
# -----------------------------------------------------------------------------
def read_simulated_r_hii(snapshot_glob):
    files = sorted(glob.glob(snapshot_glob))
    if not files:
        raise RuntimeError(f"No snapshots found matching {snapshot_glob!r}.")

    times = []
    r_hii = []
    star_mass_msun = None
    n_H_atom_cc = None
    unit_time_to_Myr = None
    time_begin_internal = None

    for f in files:
        with h5py.File(f, "r") as h:
            header = h["Header"]
            a = float(header.attrs["Scale-factor"][0])
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_M = h["Units"].attrs["Unit mass in cgs (U_M)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            if unit_time_to_Myr is None:
                unit_time_to_Myr = (u_t * u.s).to(u.Myr).value

            gas = h["PartType0"]
            stars = h["PartType4"]

            # Header/Time is the ABSOLUTE physical cosmic age at this
            # snapshot's scale factor (src/engine.c: e->time_begin =
            # cosmology->time_begin, itself the age of the Universe at
            # a_begin via the Friedmann integral) -- not elapsed run time,
            # and not a comoving quantity needing an a-factor. Subtract the
            # first snapshot's cosmic age to get elapsed time since a_begin,
            # matching the non-cosmological check's t=0-at-start convention.
            if time_begin_internal is None:
                time_begin_internal = float(header.attrs["Time"][0])
            times.append(
                (float(header.attrs["Time"][0]) - time_begin_internal)
                * unit_time_to_Myr
            )

            boxsize_comoving = header.attrs["BoxSize"][0]  # cubic box
            boxsize_physical_kpc = (
                boxsize_comoving * a ** a_scale_exponent(gas["Coordinates"])
            ) * (u_L * u.cm).to(u.kpc).value

            # r_hii from the star's own HIIRegionRadii (h_hii, comoving) is
            # frozen at the moment each gas particle was tagged and never
            # updates as that particle keeps moving -- see the non-cosmological
            # check for the full rationale. Use the current PHYSICAL position
            # of every gas particle this star has *ever* tagged instead
            # (HIIStarIDs, dimensionless, never reset when a tag expires).
            star_id = int(stars["ParticleIDs"][0])
            hii_star_ids = gas["HIIStarIDs"][...]
            ever_tagged = hii_star_ids == star_id
            if np.any(ever_tagged):
                star_pos_comoving = stars["Coordinates"][0]
                gas_pos_comoving = gas["Coordinates"][...][ever_tagged]
                a_L = a_scale_exponent(gas["Coordinates"])
                star_pos_kpc = (star_pos_comoving * a**a_L) * (u_L * u.cm).to(
                    u.kpc
                ).value
                gas_pos_kpc = (gas_pos_comoving * a**a_L) * (u_L * u.cm).to(u.kpc).value
                # Minimum-image convention in PHYSICAL coordinates, using the
                # PHYSICAL box size at this snapshot's scale factor.
                dx = gas_pos_kpc - star_pos_kpc
                dx -= boxsize_physical_kpc * np.round(dx / boxsize_physical_kpc)
                r_hii.append(float(np.max(np.linalg.norm(dx, axis=1))))
            else:
                hii_radii_comoving = stars["HIIRegionRadii"][...]
                a_L = a_scale_exponent(stars["HIIRegionRadii"])
                r_hii.append(
                    float(np.max(hii_radii_comoving * a**a_L))
                    * (u_L * u.cm).to(u.kpc).value
                )

            if star_mass_msun is None:
                masses = stars["Masses"][...]
                if len(masses) != 1:
                    raise RuntimeError(
                        f"Expected exactly 1 star (this check assumes "
                        f"star_type=single_star, matching the Q_H(mass) fit "
                        f"reimplemented below), found {len(masses)}."
                    )
                star_mass_msun = float(
                    (masses[0] * a ** a_scale_exponent(stars["Masses"]))
                    * (u_M * u.g).to(u.Msun).value
                )

            if n_H_atom_cc is None:
                # Hydrogen number density from the gas PHYSICAL density and
                # the Grackle-reported hydrogen mass fraction (see the
                # non-cosmological check for why this is read, not assumed).
                rho_comoving = gas["Densities"][...]
                a_rho = a_scale_exponent(gas["Densities"])
                rho_g_cm3 = (
                    float(np.mean(rho_comoving * a**a_rho))
                    * ((u_M * u.g) / (u_L * u.cm) ** 3).to(u.g / u.cm**3).value
                )
                X_H_raw = h["Parameters"].attrs.get(
                    "GrackleCooling:HydrogenFractionByMass"
                )
                # Parameters-group attributes are scalar byte-strings (the
                # parsed parameter text), not numeric arrays -- indexing
                # with [0] would silently take the first ASCII byte instead
                # of decoding the value.
                X_H = float(X_H_raw.decode()) if X_H_raw is not None else 0.716
                n_H_atom_cc = (rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p).to(
                    1 / u.cm**3
                )

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        star_mass_msun,
        n_H_atom_cc,
    )


def measure_ionized_temperature_K(files, min_count=20):
    """Median PHYSICAL temperature of the currently-tagged (ionized) gas,
    from the latest snapshot holding at least min_count tagged particles.

    @return (median T [K], source filename), or (None, None) if no
    snapshot holds enough tagged gas.
    """
    for filename in reversed(files):
        with h5py.File(filename, "r") as h:
            gas = h["PartType0"]
            ionized = gas["IsIonizedFlags"][...].astype(bool)
            if np.sum(ionized) < min_count:
                continue
            a = float(h["Header"].attrs["Scale-factor"][0])
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            u_to_cgs_comoving = (
                u_L / u_t
            ) ** 2  # comoving internal energy unit -> erg/g
            a_u = a_scale_exponent(gas["InternalEnergies"])
            if "HI" in gas:
                inv_mu = (
                    gas["HI"][...][ionized]
                    + 2.0 * gas["HII"][...][ionized]
                    + 0.25 * gas["HeI"][...][ionized]
                    + 0.5 * gas["HeII"][...][ionized]
                    + 0.75 * gas["HeIII"][...][ionized]
                )
                mu = 1.0 / np.clip(inv_mu, 1e-10, None)
            else:
                # No species arrays (COOLING_GRACKLE_MODE 0): this selection
                # is fully ionized by construction, so use the fully-ionized
                # primordial mu.
                mu = 0.6
            u_cgs_physical = (
                gas["InternalEnergies"][...][ionized] * a**a_u * u_to_cgs_comoving
            )
            gamma_m1 = 2.0 / 3.0
            T = (
                u_cgs_physical
                * gamma_m1
                * mu
                * const.m_p.cgs.value
                / const.k_B.cgs.value
            )
            return float(np.median(T)), filename
    return None, None


def resolve_T_ionized_K(requested_K, files):
    """Return the T_i [K] the reference is computed at: the measured value
    by default (requested_K is None), or the explicit one with a loud
    warning when it disagrees with the measurement by more than 10%."""
    T_measured_K, source = measure_ionized_temperature_K(files)
    if requested_K is None:
        if T_measured_K is None:
            raise RuntimeError(
                "Could not measure T_i (no snapshot holds enough tagged "
                "gas); pass --T-ionized-K explicitly."
            )
        print(
            f"T_ionized          : {T_measured_K:.4g} K "
            f"(measured from tagged gas, {source})"
        )
        return T_measured_K
    if T_measured_K is not None and abs(T_measured_K / requested_K - 1.0) > 0.1:
        print("\n" + "!" * 72)
        print(
            f"WARNING: --T-ionized-K {requested_K:.4g} K disagrees with the "
            f"temperature this run actually holds its tagged gas at "
            f"({T_measured_K:.4g} K, measured in {source})."
        )
        print(
            "The reference below is computed at a temperature the run did "
            "not use, so the verdict is not meaningful. Omit --T-ionized-K "
            "to use the measured value."
        )
        print("!" * 72 + "\n")
    return requested_K


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
        default="cosmo_stromgren_analytic_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.3,
        help="Relative-error tolerance at the final time for a pass (default: 0.3).",
    )
    parser.add_argument(
        "--T-ionized-K",
        type=float,
        default=None,
        dest="T_ionized_K",
        help="Ionized-gas temperature alpha_B and c_s are evaluated at. "
        "Default: measured from the run's own tagged gas.",
    )
    args = parser.parse_args()
    files = sorted(glob.glob(args.snapshot_glob))
    T_ionized_K = resolve_T_ionized_K(args.T_ionized_K, files)

    T_IONIZED_K = T_ionized_K * u.K
    ALPHA_B = alpha_b_hui_gnedin(T_IONIZED_K)
    C_S_IONIZED = sound_speed_ionized(T_IONIZED_K)

    t_sim, r_sim, star_mass_msun, n_H = read_simulated_r_hii(args.snapshot_glob)

    if star_mass_msun is None or n_H is None:
        raise RuntimeError("Could not infer star mass / n_H from the snapshots.")

    Q_H = ionizing_photon_rate(star_mass_msun)
    R_st = stromgren_radius(Q_H, n_H, ALPHA_B)

    with h5py.File(files[0], "r") as h:
        a0 = float(h["Header"].attrs["Scale-factor"][0])
        u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
        boxsize_comoving0 = h["Header"].attrs["BoxSize"][0]
        a_L = a_scale_exponent(h["PartType0/Coordinates"])
        box_half_width0 = (
            0.5 * boxsize_comoving0 * a0**a_L * (u_L * u.cm).to(u.kpc).value
        ) * u.kpc

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print(f"n_H (physical, t=0): {n_H:.4g}")
    print(f"Q_H                : {Q_H:.4g}")
    print(
        f"T_ionized (applied): {T_IONIZED_K:.4g}  ->  alpha_B = {ALPHA_B:.4g}, "
        f"c_s = {C_S_IONIZED:.4g}"
    )
    print(f"Equilibrium R_st   : {R_st.to(u.pc):.4g} = {R_st.to(u.kpc):.4g}")
    print(f"Physical box half-width (t=0): {box_half_width0.to(u.pc):.4g}")
    if R_st > box_half_width0:
        print(
            "  WARNING: R_st exceeds the physical box half-width -- the "
            "periodic box is too small to contain a full infinite-medium "
            "Stromgren sphere. The analytic curve below is an upper bound."
        )

    r_analytic = dtype_expansion_radius(t_sim, R_st, C_S_IONIZED).to(u.kpc)

    # Same validity window logic as the non-cosmological check: only compare
    # while the analytic radius still fits the (physical) box and the star
    # is alive/tracked (r_sim > 0). box_half_width0 is used throughout since
    # its scale-factor drift over these runs' short a-span is negligible
    # (<2.4% even for params_highz.yml) compared to this check's tolerance.
    valid = (r_analytic <= box_half_width0) & (r_sim > 0)
    box_exceeded_at = t_sim[r_analytic > box_half_width0]
    nonzero_idx = np.where(r_sim > 0)[0]
    dead_mask = np.zeros(len(r_sim), dtype=bool)
    if len(nonzero_idx) > 0:
        dead_mask[nonzero_idx[0] + 1 :] = r_sim[nonzero_idx[0] + 1 :] == 0
    star_dead_at = t_sim[dead_mask]
    if not np.any(valid):
        raise RuntimeError(
            "No snapshot has both an analytic radius within the physical "
            "box and an alive, actively-tracked r_sim > 0 -- no meaningful "
            "comparison point exists in this run."
        )
    last_valid = np.where(valid)[0][-1]
    t_sim_final = t_sim[last_valid]
    r_sim_final = r_sim[last_valid]
    r_analytic_final = r_analytic[last_valid]
    rel_error = float(
        (np.abs(r_sim_final - r_analytic_final) / r_analytic_final).decompose()
    )
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"Comparison time: t={t_sim_final:.4g}  r_sim={r_sim_final:.4g}  "
        f"r_analytic={r_analytic_final:.4g}  rel_error={rel_error:.2%}  [{verdict}]"
    )
    if last_valid < len(t_sim) - 1:
        reason = []
        if len(box_exceeded_at) > 0:
            reason.append(f"box half-width exceeded from t={box_exceeded_at[0]:.4g}")
        if len(star_dead_at) > 0:
            reason.append(f"star inactive/dead from t={star_dead_at[0]:.4g}")
        print(
            f"  (not the final simulated time t={t_sim[-1]:.4g} -- "
            f"{'; '.join(reason)}.)"
        )

    fig, ax = plt.subplots()
    ax.plot(
        t_sim.to(u.Myr), r_sim.to(u.pc), "o-", label="Simulation ($r_{hii}$, physical)"
    )
    ax.plot(
        t_sim.to(u.Myr),
        r_analytic.to(u.pc),
        "--",
        label="Analytic (Spitzer 1978 D-type)",
    )
    ax.axhline(
        R_st.to(u.pc).value, color="grey", linestyle=":", label="Equilibrium $R_{st}$"
    )
    ax.axhline(
        box_half_width0.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        label="Physical box half-width",
    )
    ax.axvline(
        t_sim_final.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.6,
        label="Comparison point",
    )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("HII region radius [pc] (physical)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
