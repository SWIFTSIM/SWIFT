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

The ever-tagged r_hii is gap-rejected, not the max, of ever-tagged gas
radii: HIIStarIDs is a one-way stamp (never cleared on tag expiry), so a
single particle that lapsed early and then advected outward for several
Myr can inflate a plain max() by double digits percent. Discarding
particles by a blind percentile also works, but cuts several
percent into a perfectly smooth, gap-free shell edge where nothing is
actually an outlier; walking down from the top and discarding only points
separated from the bulk by an anomalously large gap keeps the smooth case
exact while still rejecting a genuine straggler. The true max is still
reported alongside it as a secondary diagnostic.

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
        photon_integral_sum += (exp_term / n) * (x_0**2 + (2.0 * x_0) / n + 2.0 / n**2)

    prefactor = 15.0 / np.pi**4
    N_dot_ion = (L / (const.k_B * T)) * prefactor * photon_integral_sum
    return N_dot_ion.to(1 / u.s)


# -----------------------------------------------------------------------------
# Analytic Stromgren-sphere solution
# -----------------------------------------------------------------------------
# This code's own ionized-gas temperature floor is metallicity-dependent
# (min(Delta_u_ionized, u_collisional), src/cooling/grackle/
# cooling_gear_subgrid.h), unlike the classical D-type solution's flat
# 1e4 K assumption. At this example's default Z=0, the applied temperature
# is ~47,500 K (Delta_u_ionized binds; see theory/GEAR/Radiation/,
# Section "Metallicity is not the explanation..."), not 1e4 K -- both
# alpha_B (now temperature-dependent since the Hui & Gnedin 1997 fix,
# src/feedback/GEAR/radiation.c) and the isothermal sound speed driving
# the D-type expansion should be evaluated there, not at the classical
# 1e4 K reference point, or this analytic comparison is systematically
# biased low.
def alpha_b_hui_gnedin(T):
    """Case-B recombination coefficient (Hui & Gnedin 1997, MNRAS 292, 27,
    Appendix A), the same fit now used in
    radiation_get_case_b_recombination_coefficient_cgs
    (src/feedback/GEAR/radiation.c)."""
    lam = 315614.0 / T.to(u.K).value
    return (
        2.753e-14 * lam**1.5 / (1.0 + (lam / 2.740) ** 0.407) ** 2.242 * u.cm**3 / u.s
    )


def sound_speed_ionized(T_ionized):
    """Isothermal sound speed of photoionized hydrogen gas, scaled from the
    standard 1e4 K reference value (12.85 km/s) to the actual applied
    ionized temperature (c_s ~ sqrt(T) for an ideal gas at fixed
    composition)."""
    return 12.85 * u.km / u.s * np.sqrt(T_ionized / (1e4 * u.K))


def stromgren_radius(Q_H, n_H, alpha_B):
    """Equilibrium (Stromgren) radius: R_st = (3 Q_H / (4 pi alpha_B n_H^2))^(1/3)."""
    R_st3 = 3.0 * Q_H / (4.0 * np.pi * alpha_B * n_H**2)
    return R_st3.to(u.cm**3) ** (1.0 / 3.0)


def dtype_expansion_radius(t, R_st, c_s):
    """Spitzer (1978) D-type expansion: R(t) = R_st (1 + 7 c_s t / (4 R_st))^(4/7)."""
    return R_st * (1.0 + 7.0 * c_s * t / (4.0 * R_st)) ** (4.0 / 7.0)


# A point at the top of the sorted ever-tagged distances is rejected as an
# outlier only if the gap to its next-lower neighbor exceeds this many times
# the local point spacing: an isolated advected straggler produces a gap an
# order of magnitude or more above the front's own natural surface
# roughness.
GAP_REJECTION_FACTOR = 5.0
# Fraction of ever-tagged particles, counted in from the front, used to
# estimate the local point spacing scale GAP_REJECTION_FACTOR is measured
# against.
GAP_SCALE_TOP_FRACTION = 0.05
# Below this many ever-tagged particles there are too few points to estimate
# a meaningful local spacing scale; trust the max unmodified.
GAP_MIN_POINTS = 10


def robust_ever_tagged_radius(distances_kpc):
    """Robust "ever-tagged" HII front radius from ever-tagged particle
    distances, via gap-aware outlier rejection.

    HIIStarIDs is a one-way stamp (radiation_reset_part_ionized_tag never
    clears it -- the stamp intentionally marks the swept shell), so a
    single particle whose tag lapsed early and then advected outward for
    several Myr can dominate a plain max() and inflate the reported extent
    by double digits percent. A blind percentile trim removes such
    stragglers too, but also cuts several percent into a
    perfectly smooth, gap-free shell edge where every particle is
    legitimate. Instead, walk down from the largest distance and discard a
    point only if the gap to its next-lower neighbor exceeds
    GAP_REJECTION_FACTOR times the local point spacing (the median spacing
    among the top GAP_SCALE_TOP_FRACTION of points); stop at the first
    point that is not an outlier by this test.

    @param distances_kpc Ever-tagged particle distances from the source (kpc).
    @return (front radius after gap rejection, true max), both in kpc.
    """
    sorted_desc = np.sort(distances_kpc)[::-1]
    true_max = float(sorted_desc[0])
    n = len(sorted_desc)
    if n < GAP_MIN_POINTS:
        return true_max, true_max

    n_top = min(max(2, int(np.ceil(GAP_SCALE_TOP_FRACTION * n))), n - 1)
    top_gaps = -np.diff(sorted_desc[: n_top + 1])
    local_scale = np.median(top_gaps)
    if local_scale <= 0:
        nonzero = top_gaps[top_gaps > 0]
        local_scale = np.median(nonzero) if len(nonzero) else 0.0

    i = 0
    while (
        i < n - 1
        and local_scale > 0
        and (sorted_desc[i] - sorted_desc[i + 1]) > GAP_REJECTION_FACTOR * local_scale
    ):
        i += 1
    return float(sorted_desc[i]), true_max


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
    r_hii_max = []
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

        # r_hii from the star's own HIIRegionRadii (h_hii) is frozen at the
        # moment each gas particle was tagged -- it never updates as that
        # particle keeps moving afterwards. Since ionized gas has substantial
        # outward velocity (confirmed directly: tagged particles show mean
        # v_r several times larger than untagged gas, growing with time),
        # h_hii systematically undercounts the true current extent of the
        # ionized/expanding shell -- by up to ~60-80% at t~2-3 Myr in the
        # e10 check. Instead, use the current position of every gas particle
        # this star has *ever* tagged: HIIStarIDs is set once when a
        # particle is first ionized and, unlike IsIonizedFlags, is never
        # reset when that tag later expires (see
        # radiation_reset_part_ionized_tag, src/feedback/GEAR/radiation.c),
        # so it survives the rebuild-cadence tag/re-tag cycle. A plain max()
        # over these positions is fragile to a single lapsed-tag particle
        # that keeps advecting outward, so robust_ever_tagged_radius gap-
        # rejects it and reports the front position as the primary value,
        # with the true max as a secondary diagnostic.
        star_id = int(data.stars.particle_ids[0])
        ever_tagged = data.gas.hiistar_ids == star_id
        if np.any(ever_tagged):
            star_pos = data.stars.coordinates[0].to("kpc").value
            gas_pos = data.gas.coordinates[ever_tagged].to("kpc").value
            box_kpc = data.metadata.boxsize.to("kpc").value
            # Minimum-image convention: a periodic box means the naive
            # difference can be wrong by a whole box length once gas has
            # wrapped around, which would corrupt the extent with a single
            # spurious huge distance -- match the simulation's own
            # nearest-image distance convention here.
            dx = gas_pos - star_pos
            dx -= box_kpc * np.round(dx / box_kpc)
            distances = np.linalg.norm(dx, axis=1)
            r_robust, r_max = robust_ever_tagged_radius(distances)
            r_hii.append(r_robust)
            r_hii_max.append(r_max)
        else:
            r_now = float(np.max(data.stars.hiiregion_radii).to("kpc").value)
            r_hii.append(r_now)
            r_hii_max.append(r_now)

        if star_mass_msun is None:
            n_stars = len(data.stars.masses)
            if n_stars != 1:
                raise RuntimeError(
                    f"Expected exactly 1 star (this check assumes "
                    f"star_type=single_star, matching the Q_H(mass) fit "
                    f"reimplemented below), found {n_stars}. Re-run with "
                    f"star_type=single_star, n_stars=1, or fix this script's "
                    f"Q_H fit for whatever star_type was actually used before "
                    f"trusting this comparison."
                )
            star_mass_msun = float(np.max(data.stars.masses).to("Msun").value)
        if n_H_atom_cc is None:
            # Hydrogen number density from the gas density and the
            # Grackle-reported hydrogen mass fraction. Read the actual value
            # this run used (GrackleCooling:HydrogenFractionByMass, stored in
            # every snapshot's Parameters group) rather than assuming a
            # fixed constant -- a previous version hardcoded 0.716 here,
            # silently disagreeing with e.g. this example's actual 0.76 and
            # biasing n_H (and R_st ~ n_H^{-2/3}) by a few percent. Only fall
            # back to the primordial default if the run used a different
            # cooling backend that doesn't expose this parameter.
            rho_g_cm3 = float(np.mean(data.gas.densities).to("g/cm**3").value)
            X_H_raw = data.metadata.parameters.get(
                "GrackleCooling:HydrogenFractionByMass"
            )
            X_H = float(X_H_raw) if X_H_raw is not None else 0.716
            n_H_atom_cc = (rho_g_cm3 * u.g / u.cm**3 * X_H / const.m_p).to(1 / u.cm**3)

    return (
        u.Quantity(times, u.Myr),
        u.Quantity(r_hii, u.kpc),
        u.Quantity(r_hii_max, u.kpc),
        star_mass_msun,
        n_H_atom_cc,
        boxsize_kpc * u.kpc,
    )


def measure_ionized_temperature_K(files, min_count=20):
    """Median temperature of the currently-tagged (ionized) gas, from the
    latest snapshot holding at least min_count tagged particles.

    The run's actual T_i depends on build flags
    (IONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K) and
    metallicity; measuring it removes the silent tens-of-percent error of
    computing the reference at a temperature the run never used.

    @return (median T [K], source filename), or (None, None) if no
    snapshot holds enough tagged gas.
    """
    import h5py

    for filename in reversed(files):
        with h5py.File(filename, "r") as h:
            gas = h["PartType0"]
            ionized = gas["IsIonizedFlags"][:].astype(bool)
            if np.sum(ionized) < min_count:
                continue
            u_L = h["Units"].attrs["Unit length in cgs (U_L)"][0]
            u_t = h["Units"].attrs["Unit time in cgs (U_t)"][0]
            u_to_cgs = (u_L / u_t) ** 2  # internal energy unit -> erg/g
            if "HI" in gas:
                inv_mu = (
                    gas["HI"][ionized]
                    + 2.0 * gas["HII"][ionized]
                    + 0.25 * gas["HeI"][ionized]
                    + 0.5 * gas["HeII"][ionized]
                    + 0.75 * gas["HeIII"][ionized]
                )
                mu = 1.0 / np.clip(inv_mu, 1e-10, None)
            else:
                # No species arrays (COOLING_GRACKLE_MODE 0): this selection
                # is fully ionized by construction, so use the fully-ionized
                # primordial mu.
                mu = 0.6
            u_cgs = gas["InternalEnergies"][ionized] * u_to_cgs
            gamma_m1 = 2.0 / 3.0
            T = u_cgs * gamma_m1 * mu * const.m_p.cgs.value / const.k_B.cgs.value
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
        default="stromgren_analytic_check.png",
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
        "Default: measured from the run's own tagged gas, which always "
        "matches whatever build flag/metallicity convention the run used. "
        "An explicit value overrides the measurement but warns loudly if "
        "the two disagree by >10%%.",
    )
    args = parser.parse_args()
    T_ionized_K = resolve_T_ionized_K(
        args.T_ionized_K, sorted(glob.glob(args.snapshot_glob))
    )

    T_IONIZED_K = T_ionized_K * u.K
    ALPHA_B = alpha_b_hui_gnedin(T_IONIZED_K)
    C_S_IONIZED = sound_speed_ionized(T_IONIZED_K)

    t_sim, r_sim, r_sim_max, star_mass_msun, n_H, boxsize = read_simulated_r_hii(
        args.snapshot_glob
    )

    if star_mass_msun is None or n_H is None:
        raise RuntimeError("Could not infer star mass / n_H from the snapshots.")

    Q_H = ionizing_photon_rate(star_mass_msun)
    R_st = stromgren_radius(Q_H, n_H, ALPHA_B)
    box_half_width = 0.5 * boxsize

    print(f"Star mass          : {star_mass_msun:.3f} Msun")
    print(f"n_H                : {n_H:.4g}")
    print(f"Q_H                : {Q_H:.4g}")
    print(
        f"T_ionized (applied): {T_IONIZED_K:.4g}  ->  alpha_B = {ALPHA_B:.4g}, "
        f"c_s = {C_S_IONIZED:.4g}"
    )
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

    r_analytic = dtype_expansion_radius(t_sim, R_st, C_S_IONIZED).to(u.kpc)

    # The comparison is only meaningful while (a) the analytic curve is
    # still predicted to fit inside the box -- beyond that, periodic
    # self-interaction/wrapping makes the simulated and infinite-medium
    # analytic solutions incomparable, not wrong -- and (b) the star is
    # still alive and actively tracked: h_hii resets to 0 once a star
    # dies or ages past HII_max_age (feedback_will_do_feedback),
    # which is bookkeeping, not the region physically vanishing. Picking
    # the LAST snapshot that satisfies both, rather than blindly using
    # the final simulated time, makes this check robust to runs that
    # continue well past either limit (e.g. an acceptance run taken to
    # time_end long after the star has died).
    valid = (r_analytic <= box_half_width) & (r_sim > 0)
    box_exceeded_at = t_sim[r_analytic > box_half_width]
    # r_sim == 0 before the region has ever formed is not "death" -- only
    # flag a zero that follows a nonzero sample as the star going dark.
    nonzero_idx = np.where(r_sim > 0)[0]
    dead_mask = np.zeros(len(r_sim), dtype=bool)
    if len(nonzero_idx) > 0:
        dead_mask[nonzero_idx[0] + 1 :] = r_sim[nonzero_idx[0] + 1 :] == 0
    star_dead_at = t_sim[dead_mask]
    if not np.any(valid):
        raise RuntimeError(
            "No snapshot has both an analytic radius within the box and an "
            "alive, actively-tracked r_sim > 0 -- no meaningful comparison "
            "point exists in this run. Try a bigger boxsize or a higher "
            "gas_density (shrinks R_st ~ n_H^-2/3), or check the run "
            "completed at least one HII rebuild before the star died."
        )
    last_valid = np.where(valid)[0][-1]
    t_sim_final = t_sim[last_valid]
    r_sim_final = r_sim[last_valid]
    r_sim_max_final = r_sim_max[last_valid]
    r_analytic_final = r_analytic[last_valid]
    rel_error = float(
        (np.abs(r_sim_final - r_analytic_final) / r_analytic_final).decompose()
    )
    verdict = "PASS" if rel_error <= args.tol else "FAIL"
    print(
        f"Comparison time: t={t_sim_final:.4g}  r_sim={r_sim_final:.4g}  "
        f"r_analytic={r_analytic_final:.4g}  rel_error={rel_error:.2%}  [{verdict}]"
    )
    print(
        f"  (secondary diagnostic, true max ever-tagged extent: "
        f"{r_sim_max_final:.4g})"
    )
    if last_valid < len(t_sim) - 1:
        reason = []
        if len(box_exceeded_at) > 0:
            reason.append(f"box half-width exceeded from t={box_exceeded_at[0]:.4g}")
        if len(star_dead_at) > 0:
            reason.append(f"star inactive/dead from t={star_dead_at[0]:.4g}")
        print(
            f"  (not the final simulated time t={t_sim[-1]:.4g} -- "
            f"{'; '.join(reason)}; later snapshots are not comparable to "
            f"the infinite-medium analytic solution.)"
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
    ax.axhline(
        box_half_width.to(u.pc).value,
        color="firebrick",
        linestyle=":",
        label="Box half-width",
    )
    ax.axvline(
        t_sim_final.to(u.Myr).value,
        color="black",
        linestyle="-.",
        alpha=0.6,
        label="Comparison point",
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
