/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
/**
 * @file src/feedback/GEAR/radiation_getters.c
 * @brief Star-emission getters for GEAR radiation feedback: IMF-integrated
 * and raw (single-mass) bolometric luminosity, ionization rate and mean
 * excess photon energy, from the 1D (mass-only) and 2D (mass x
 * metallicity) interpolation tables built by radiation_read_data().
 */

/* Config parameters. */
#include <config.h>

/* Include header */
#include "error.h"
#include "inline.h"
#include "interpolation.h"
#include "minmax.h"
#include "radiation.h"

/**
 * @brief Abort with a clear message if #rad's table dimensionality does not
 * match what the calling getter expects.
 *
 * For an expected 2D (mass x metallicity) table, checks #is_active before
 * #is_2d: a #radiation with no table loaded at all (photoionization and
 * radiation pressure both off; radiation_zero_pointers()) also has
 * is_2d = 0, and reporting that case as "holds a mass-only (1D) table"
 * would be actively misleading: there is no table of either
 * dimensionality, not a 1D one. An expected-1D check does not need this:
 * every 1D getter below is only ever reached from an #is_active branch
 * upstream already.
 *
 * @param rad The #radiation model.
 * @param expect_2d Nonzero if the caller needs a mass x metallicity (2D)
 * table; 0 if it needs a mass-only (1D) table.
 * @param caller Name of the calling getter, for the error message.
 */
__attribute__((always_inline)) INLINE static void
radiation_check_dimensionality(const struct radiation *rad, int expect_2d,
                               const char *caller) {
  if (expect_2d) {
    if (!rad->is_active) {
      error("%s called with no radiation table loaded (#rad->is_active = 0).",
            caller);
    }
    if (!rad->is_2d) {
      error(
          "%s requires a mass x metallicity (2D) radiation table: #rad "
          "holds a mass-only (1D) table.",
          caller);
    }
  } else if (rad->is_2d) {
    error(
        "%s has no metallicity argument and cannot read a mass x "
        "metallicity (2D) radiation table.",
        caller);
  }
}

/**
 * @brief Floor a metallicity mass fraction and return its log10, for a 2D
 * getter's log_z argument.
 *
 * Floors at #RADIATION_LOG_FLOOR_CGS (the same 1e-300 floor
 * radiation_read_cgs_array() applies to a table VALUE, reused here for a
 * QUERY metallicity, matching pychem's own convention of flooring before
 * any log10() that could otherwise see an exact 0
 * (PyChemInitTable.parsec.interpolate._safe_log10()/_LOG_FLOOR, also
 * 1e-300): Z=0 (e.g. pristine or Pop III gas) would otherwise hand log10()
 * a genuine 0 and produce -inf, undefined behaviour once #interpolate_2d
 * truncates it to an int. The metallicity axis's own boundary condition
 * (boundary_condition_const, set in radiation_build_tables()) already
 * clamps any out-of-range query to the table's lowest tabulated row
 * regardless of how far below it Z sits, so this floor changes nothing
 * about which row a very-low-Z star lands on; it only keeps the log10()
 * call itself finite. It does not, however, preserve mass-axis
 * interpolation for those stars: #interpolate_2d's out-of-range branch
 * (metallicity below the table's lowest tabulated value, which is exactly
 * what a floored Z=0 query is) returns a single clamped cell rather than
 * blending between neighbouring mass points, so every star at or below
 * the table's lowest tabulated Z is quantized onto the mass grid's own
 * resolution instead of interpolated: a real, if small, discontinuity
 * for the pristine/Pop III population #radiation_get_log_metallicity()
 * exists specifically to route down this path.
 *
 * Computed in double throughout, narrowed to float only on the return:
 * doing the max() in float would first narrow #RADIATION_LOG_FLOOR_CGS
 * (1e-300) to exactly 0.0f (float32's smallest denormal is ~1.4e-45),
 * silently defeating the floor and reintroducing log10(0).
 *
 * @param Z Metallicity mass fraction (may be exactly 0).
 * @return log10(max(Z, #RADIATION_LOG_FLOOR_CGS)).
 */
float radiation_get_log_metallicity(float Z) {
  const double Z_floored = max((double)Z, RADIATION_LOG_FLOOR_CGS);
  return (float)log10(Z_floored);
}

/**
 * @brief Get the IMF-averaged bolometric luminosity per mass.
 *
 * Reads #rad->integrated.luminosities directly (no exponentiation): unlike
 * the raw table, the IMF-integrated table stays in linear value space.
 * See radiation_build_tables()'s doxygen for why pychem's log-log
 * convention does not apply to this SWIFT-side cumulative-integral
 * quantity.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral(const struct radiation *rad,
                                               float log_m1, float log_m2) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  float luminosity_1 = interpolate_1d(&rad->integrated.luminosities, log_m1);
  float luminosity_2 = interpolate_1d(&rad->integrated.luminosities, log_m2);
  return luminosity_2 - luminosity_1;
};

/**
 * @brief Get the non-IMF-integrated bolometric luminosity at a given mass.
 *
 * #rad->raw.luminosities holds log10(luminosity), pychem's own
 * log10(mass)-vs-log10(value) convention (see radiation_build_tables()'s
 * doxygen); the interpolated log-value is exponentiated back here. The
 * narrowing to float happens implicitly on return (exp10() itself returns
 * double); Luminosity is never 0 in this table (unlike Q_H/DotEExcess), so
 * there is no floor-underflow case to worry about here.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return The bolometric luminosity, internal units.
 */
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                          float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.luminosities, log_m));
};

/**
 * @brief Get the IMF-averaged ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param The ionization rate;
 */
double radiation_get_ionization_rate_from_integral(const struct radiation *rad,
                                                   float log_m1, float log_m2) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  return dot_N_ion_2 - dot_N_ion_1;
};

/**
 * @brief Get the non-IMF-integrated ionization rate at a given mass.
 *
 * #rad->raw.dot_N_ion holds log10(dot_N_ion /
 * #RADIATION_DOT_N_ION_TABLE_SCALING in internal units) (see
 * radiation_build_tables()'s doxygen for the log-log storage convention);
 * exp10() undoes the log, exactly recovering the pre-log-transform
 * (already-scaled-down) internal value, then the existing
 * *RADIATION_DOT_N_ION_TABLE_SCALING undoes the scaling as before. The
 * narrowing to `float` below is load-bearing, not cosmetic: near
 * #RADIATION_LOG_FLOOR_CGS (1e-300, applied before the internal-unit
 * conversion divides it further down), float32's underflow floor
 * (~1e-45) is reached well before float64's, so this narrowing is what
 * makes a below-ionization-threshold query reliably return exactly
 * 0.0f; the exact query mass at which that happens depends on the run's
 * own unit system. See radiation_read_cgs_array()'s own doxygen for the
 * reasoning behind #RADIATION_LOG_FLOOR_CGS's specific value; this getter
 * is why it needs to be that extreme.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return The ionization rate, internal units.
 */
double radiation_get_ionization_rate_from_raw(const struct radiation *rad,
                                              float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const float dot_N_ion_scaled =
      (float)exp10(interpolate_1d(&rad->raw.dot_N_ion, log_m));
  return (double)dot_N_ion_scaled * RADIATION_DOT_N_ION_TABLE_SCALING;
};

/**
 * @brief Get the IMF-averaged, Q-weighted mean excess photon energy above
 * the 13.6 eV HI ionization threshold, for a population over a mass window.
 *
 * Ratio of the integrated dot_E_excess and dot_N_ion tables, taken directly
 * on the raw (still /RADIATION_DOT_N_ION_TABLE_SCALING) interpolated
 * values rather than through their public accessors: the scaling constant
 * multiplies both tables identically, so it cancels in the ratio without
 * ever needing to be undone, leaving a result in cgs erg (see
 * #radiation_read_mean_excess_photon_energy_array's own doxygen for why
 * that unit convention holds).
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return Q-weighted mean excess photon energy in cgs erg, or 0 if no
 * ionizing photons are produced over the window (dot_N_ion difference is
 * 0, e.g. no alive ionizing stars).
 */
double radiation_get_mean_excess_photon_energy_HI_from_integral(
    const struct radiation *rad, float log_m1, float log_m2) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1);
  const double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2);
  const double delta_dot_N_ion = dot_N_ion_2 - dot_N_ion_1;

  /* The cumulative table is monotonically non-decreasing in mass, so any
     non-positive difference (no alive ionizing stars in this window, a
     zero-width window, or roundoff noise between two nearly-equal table
     entries) is degenerate. Guard against dividing by it rather than
     testing for exact 0, which a near-cancellation could slip past and
     amplify into a meaningless huge or negative result. */
  if (delta_dot_N_ion <= 0.) return 0.;

  const double dot_E_excess_1 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m1);
  const double dot_E_excess_2 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m2);
  const double delta_dot_E_excess = dot_E_excess_2 - dot_E_excess_1;

  return delta_dot_E_excess / delta_dot_N_ion;
};

/**
 * @brief Get the non-IMF-integrated mean excess photon energy above the
 * 13.6 eV HI ionization threshold, for a single star of a given mass.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_integral on the
 * raw (single-mass) tables: the ratio of the raw dot_E_excess and dot_N_ion
 * tables, taken directly on their still-/RADIATION_DOT_N_ION_TABLE_SCALING
 * values, for the same reason that ratio is exact and comes out in cgs erg
 * there. Both tables now hold log10(value) (see radiation_build_tables()'s
 * doxygen for the log-log storage convention); each is exponentiated
 * back, narrowed to `float`, before the ratio: the same load-bearing
 * float32-underflow narrowing #radiation_get_ionization_rate_from_raw
 * uses to return exactly 0.0f below the table's native ionization
 * threshold, which is what makes the `dot_N_ion <= 0.` guard below
 * actually trigger there instead of dividing by a near-zero double.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return Mean excess photon energy in cgs erg, or 0 if this mass produces
 * no ionizing photons (dot_N_ion(log_m) <= 0).
 */
double radiation_get_mean_excess_photon_energy_HI_from_raw(
    const struct radiation *rad, float log_m) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const double dot_N_ion =
      (float)exp10(interpolate_1d(&rad->raw.dot_N_ion, log_m));

  /* Same degenerate-ratio guard as the integrated getter above. */
  if (dot_N_ion <= 0.) return 0.;

  const double dot_E_excess =
      (float)exp10(interpolate_1d(&rad->raw.dot_E_excess, log_m));
  return dot_E_excess / dot_N_ion;
};

/**
 * @brief Get the non-IMF-integrated bolometric luminosity at a given mass
 * and metallicity, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_luminosities_from_raw exactly, on the 2D table
 * instead of the 1D one. See that getter's own doxygen for the log-log
 * storage convention. Not capped by #main_sequence_lifetime_2d, even
 * though Luminosity is collapsed over the same main-sequence window as
 * Q_H/DotEExcess (pychem's single group-level "time_collapse" attribute
 * applies uniformly, there is no per-field exception): a post-main-sequence
 * star remains genuinely luminous, and truncating radiation pressure (which
 * #radiation_get_star_physical_radiation_pressure() derives from L_bol) to
 * exactly 0 at TAMS would be a worse approximation than over-extending its
 * main-sequence-averaged luminosity. Q_H/DotEExcess are capped instead
 * because an ionizing photon rate that no longer matches the star's real
 * state actively misleads the HII-region budget, whereas an over-extended
 * L_bol only makes radiation pressure conservative.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @return The bolometric luminosity, internal units.
 */
float radiation_get_luminosities_from_raw_2d(const struct radiation *rad,
                                             float log_z, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return (float)exp10(interpolate_2d(&rad->raw.luminosities_2d, log_z, log_m));
};

/**
 * @brief Return whether a star has evolved past MainSequenceLifetime(Z, M)
 * in the 2D ("M,Z") table.
 *
 * Shared by #radiation_get_ionization_rate_from_raw_2d and
 * #radiation_get_mean_excess_photon_energy_HI_from_raw_2d, which both cap
 * their output to exactly 0 once this returns true. See either getter's
 * own doxygen for why.
 *
 * @param rad The #radiation model (must hold an active 2D table).
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's current age, in Myr, ZAMS-anchored; see
 * #radiation_get_ionization_rate_from_raw_2d's own doxygen for why this
 * normalization matters.
 * @return 1 if @p star_age_myr exceeds MainSequenceLifetime(Z, M), 0
 * otherwise.
 */
__attribute__((always_inline)) INLINE static int
radiation_is_past_main_sequence_2d(const struct radiation *rad, float log_z,
                                   float log_m, float star_age_myr) {
  const float ms_lifetime_myr = (float)exp10(
      interpolate_2d(&rad->raw.main_sequence_lifetime_2d, log_z, log_m));
  return star_age_myr > ms_lifetime_myr;
}

/**
 * @brief Get the non-IMF-integrated ionization rate at a given mass,
 * metallicity and stellar age, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_ionization_rate_from_raw exactly on the 2D table,
 * including the load-bearing float32-underflow narrowing documented on
 * that getter; @p star_age_myr adds one gate on top: MainSequenceLifetime
 * (Z, M) is pychem's own documented time-averaging window for the table's
 * Q_H (see #radiation's raw.main_sequence_lifetime_2d doxygen), so past
 * that age the tabulated Q_H no longer describes the star's real state and
 * this returns exactly 0.0 rather than an over-extended value. This cap is
 * distinct from, and independent of, GEAR's own (typically longer) Poirier
 * stellar-lifetime death and #feedback_properties.HII_max_age, both of
 * which still gate this star's feedback independently upstream
 * (feedback_common.c): PARSEC's own main-sequence duration can run several
 * times shorter than Poirier's at some grid points, so this cap can fire
 * well before either of those does.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's current age, in Myr, anchored at the ZAMS
 * (birth of the star as a hydrogen-burning main-sequence object), matching
 * MainSequenceLifetime's own ZAMS-to-TAMS definition. NOT anchored at
 * protostellar formation. A caller deriving this from SWIFT's own star
 * particle age (which starts at particle birth/spawn, including any
 * pre-main-sequence contraction phase the model represents) must confirm
 * that normalization matches before wiring a call site; see pychem's own
 * warning on this in its MainSequenceLifetime dataset description.
 * @return The ionization rate, internal units, or exactly 0.0 if
 * @p star_age_myr exceeds the table's MainSequenceLifetime(Z, M).
 */
double radiation_get_ionization_rate_from_raw_2d(const struct radiation *rad,
                                                 float log_z, float log_m,
                                                 float star_age_myr) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);

  if (radiation_is_past_main_sequence_2d(rad, log_z, log_m, star_age_myr))
    return 0.;

  const float dot_N_ion_scaled =
      (float)exp10(interpolate_2d(&rad->raw.dot_N_ion_2d, log_z, log_m));
  return (double)dot_N_ion_scaled * RADIATION_DOT_N_ION_TABLE_SCALING;
};

/**
 * @brief Get the non-IMF-integrated mean excess photon energy above the
 * 13.6 eV HI ionization threshold, at a given mass, metallicity and
 * stellar age, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_raw exactly on
 * the 2D table, including the ratio-of-raw-tables construction and its
 * degenerate-ratio guard. Gated by the same @p star_age_myr /
 * MainSequenceLifetime(Z, M) cap as
 * #radiation_get_ionization_rate_from_raw_2d, for the same reason
 * (DotEExcess is Q_H * mean excess energy, so it shares Q_H's
 * main-sequence-averaging convention): applied first, before either raw
 * table is read, so a past-main-sequence query short-circuits to 0 without
 * needing the degenerate-ratio guard to catch it.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's current age, in Myr, ZAMS-anchored; see
 * #radiation_get_ionization_rate_from_raw_2d's own doxygen for why this
 * normalization matters.
 * @return Mean excess photon energy in cgs erg, or 0 if this (mass,
 * metallicity, age) produces no ionizing photons (dot_N_ion <= 0, or
 * @p star_age_myr exceeds MainSequenceLifetime(Z, M)).
 */
double radiation_get_mean_excess_photon_energy_HI_from_raw_2d(
    const struct radiation *rad, float log_z, float log_m, float star_age_myr) {

  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);

  if (radiation_is_past_main_sequence_2d(rad, log_z, log_m, star_age_myr))
    return 0.;

  const double dot_N_ion =
      (float)exp10(interpolate_2d(&rad->raw.dot_N_ion_2d, log_z, log_m));

  /* Same degenerate-ratio guard as the 1D getter. */
  if (dot_N_ion <= 0.) return 0.;

  const double dot_E_excess =
      (float)exp10(interpolate_2d(&rad->raw.dot_E_excess_2d, log_z, log_m));
  return dot_E_excess / dot_N_ion;
};

/**
 * @brief Nudge a mass-axis log-mass query strictly inside a 2D IMF-
 * integrated table's own top edge, so an exact-mass_max query
 * deterministically takes the blended (Z-interpolated), not
 * boundary-clamped, branch of #interpolate_2d.
 *
 * A query at @p log_m exactly equal to @p interp's own top edge (`j ==
 * Ny - 1` exactly) takes #interpolate_2d's out-of-range branch, which
 * floor-snaps the Z axis to the nearest tabulated row instead of blending
 * it, unlike the bottom edge (`j == 0`), which is safely in-range. This
 * is not a rare corner case: stellar_evolution_compute_preSN_feedback_
 * spart() clamps its upper mass bound to sm->imf.mass_max, and the default
 * mass_sup_scheme_end_step scheme routinely returns exactly that value for
 * an early-age population (before its first star has died); see Decision
 * #1 in design-radiation-integrated-table-migration.md. The nudge is a
 * #RADIATION_2D_EDGE_EPS relative fraction of the table's own mass-axis
 * span, physically negligible.
 *
 * @param interp The 2D IMF-integrated table the caller is about to query.
 * @param log_m The mass-axis query, in log10.
 * @return @p log_m, or the table's top edge minus a small margin,
 * whichever is smaller.
 */
__attribute__((always_inline)) INLINE static float radiation_nudge_mass_edge_2d(
    const struct interpolation_2d *interp, float log_m) {
  const float cap = interp->ymin + (interp->Ny - 1) * interp->dy *
                                       (1.f - RADIATION_2D_EDGE_EPS);
  return min(log_m, cap);
}

/**
 * @brief Get the IMF-averaged bolometric luminosity per mass, at a given
 * metallicity, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_luminosities_from_integral exactly, blended across
 * the metallicity axis at fixed @p log_z. See that getter's own doxygen
 * for the two-point-subtraction shape, and this file's own #radiation_
 * nudge_mass_edge_2d for the top-edge-boundary fix this 2D getter needs
 * that the 1D one does not. Below the table's own native mass floor (e.g.
 * a PopIII-like table with mass_min above the IMF's own mass_min), the
 * integrated table is flat at exactly 0 (cumulative-from-Mmin, zero at
 * every native floor cell) with no Z-dependence, so the subtraction stays
 * exact there too.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral_2d(const struct radiation *rad,
                                                  float log_z, float log_m1,
                                                  float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.luminosities_2d;
  const float luminosity_1 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m1));
  const float luminosity_2 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m2));
  return luminosity_2 - luminosity_1;
};

/**
 * @brief Get the IMF-averaged ionization rate per mass, at a given
 * metallicity, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_ionization_rate_from_integral exactly. See
 * #radiation_get_luminosities_from_integral_2d's own doxygen for the
 * shared 2D-specific caveats (top-edge nudge, sub-native-floor flat
 * region).
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return The ionization rate, internal units.
 */
double radiation_get_ionization_rate_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.dot_N_ion_2d;
  const double dot_N_ion_1 =
      interpolate_2d(interp, log_z,
                     radiation_nudge_mass_edge_2d(interp, log_m1)) *
      RADIATION_DOT_N_ION_TABLE_SCALING;
  const double dot_N_ion_2 =
      interpolate_2d(interp, log_z,
                     radiation_nudge_mass_edge_2d(interp, log_m2)) *
      RADIATION_DOT_N_ION_TABLE_SCALING;
  return dot_N_ion_2 - dot_N_ion_1;
};

/**
 * @brief Get the IMF-averaged, Q-weighted mean excess photon energy above
 * the 13.6 eV HI ionization threshold, at a given metallicity, from a 2D
 * ("M,Z") table.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_integral
 * exactly, including the degenerate-ratio guard. More load-bearing here
 * than in the 1D case, since Z-axis blending introduces its own roundoff on
 * top of the mass-axis interpolation roundoff the guard already exists to
 * catch. See #radiation_get_luminosities_from_integral_2d's own doxygen
 * for the shared 2D-specific caveats (top-edge nudge, sub-native-floor
 * flat region).
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return Q-weighted mean excess photon energy in cgs erg, or 0 if no
 * ionizing photons are produced over the window.
 */
double radiation_get_mean_excess_photon_energy_HI_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);

  const struct interpolation_2d *dot_N_ion_interp =
      &rad->integrated.dot_N_ion_2d;
  const double dot_N_ion_1 =
      interpolate_2d(dot_N_ion_interp, log_z,
                     radiation_nudge_mass_edge_2d(dot_N_ion_interp, log_m1));
  const double dot_N_ion_2 =
      interpolate_2d(dot_N_ion_interp, log_z,
                     radiation_nudge_mass_edge_2d(dot_N_ion_interp, log_m2));
  const double delta_dot_N_ion = dot_N_ion_2 - dot_N_ion_1;

  /* Same degenerate-ratio guard as the 1D getter. */
  if (delta_dot_N_ion <= 0.) return 0.;

  const struct interpolation_2d *dot_E_excess_interp =
      &rad->integrated.dot_E_excess_2d;
  const double dot_E_excess_1 =
      interpolate_2d(dot_E_excess_interp, log_z,
                     radiation_nudge_mass_edge_2d(dot_E_excess_interp, log_m1));
  const double dot_E_excess_2 =
      interpolate_2d(dot_E_excess_interp, log_z,
                     radiation_nudge_mass_edge_2d(dot_E_excess_interp, log_m2));
  const double delta_dot_E_excess = dot_E_excess_2 - dot_E_excess_1;

  return delta_dot_E_excess / delta_dot_N_ion;
};

/**
 * @brief Get the population-level main-sequence-lifetime-capped upper mass
 * bound, at a given metallicity and population age, from a 2D ("M,Z")
 * table's MainSequenceLifetimeInverse dataset.
 *
 * Implements the 4-step query-time algorithm in design-radiation-
 * integrated-table-migration.md's "Data layout" section: brackets the two
 * native metallicity rows #radiation.longest_ms_lifetime_myr is indexed by
 * (duplicating #interpolate_2d's own index arithmetic against the NATIVE
 * log10(Z) grid stored in #radiation.ms_lifetime_inverse_log_z_min/_step;
 * deliberately not #radiation.raw.main_sequence_lifetime_inverse_2d's own
 * xmin/dx, which describe its OUTPUT-resampled grid, a different, finer
 * index space; see that field's own doxygen), takes the min() of their two
 * longest-tabulated-MS-lifetimes (not a blend: Decision #6 says a blended
 * threshold can admit a query one row's own Excluded mask had already
 * rejected), and gates on both that threshold and #radiation.age_max_myr
 * before trusting the interpolated value.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param star_age_myr The population's age in Myr.
 * @param m_min Returned when no MS-active mass remains at this (Z, age):
 * the caller's own floor mass (sm->imf.mass_min in every call site so far).
 * @return m_MS_end_step, in Msun: the mass whose PARSEC main-sequence
 * lifetime equals @p star_age_myr at @p log_z, or @p m_min if none of the
 * table's masses are still on the main sequence at that age.
 */
float radiation_get_ms_lifetime_inverse_mass_2d(const struct radiation *rad,
                                                float log_z, float star_age_myr,
                                                float m_min) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);

  const float z_idx_f = (log_z - rad->ms_lifetime_inverse_log_z_min) /
                        rad->ms_lifetime_inverse_log_z_step;
  const int midx = max((int)z_idx_f, 0);
  const int z_lo = min(midx, rad->ms_lifetime_inverse_n_metallicity - 1);
  const int z_hi = min(z_lo + 1, rad->ms_lifetime_inverse_n_metallicity - 1);

  const float threshold_myr = min(rad->longest_ms_lifetime_myr[z_lo],
                                  rad->longest_ms_lifetime_myr[z_hi]);

  if (star_age_myr > rad->age_max_myr || star_age_myr > threshold_myr)
    return m_min;

  return (float)exp10(
      interpolate_2d(&rad->raw.main_sequence_lifetime_inverse_2d, log_z,
                     log10f(star_age_myr)));
}

/**
 * @brief Get a single star's bolometric luminosity at a given mass,
 * dispatching on #rad->is_2d between the 1D (mass-only) and 2D (mass x
 * metallicity) raw tables.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @return The bolometric luminosity, internal units.
 */
float radiation_get_star_luminosity(const struct radiation *rad, float log_m,
                                    float log_z) {
  if (rad->is_2d) {
    return radiation_get_luminosities_from_raw_2d(rad, log_z, log_m);
  }
  return radiation_get_luminosities_from_raw(rad, log_m);
}

/**
 * @brief Get a single star's ionization rate at a given mass, dispatching
 * on #rad->is_2d between the 1D (mass-only) and 2D (mass x metallicity) raw
 * tables.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @param star_age_myr The star's current age, in Myr, ZAMS-anchored (see
 * #radiation_get_ionization_rate_from_raw_2d's own doxygen); used only if
 * #rad holds a 2D table.
 * @return The ionization rate, internal units.
 */
double radiation_get_star_ionization_rate(const struct radiation *rad,
                                          float log_m, float log_z,
                                          float star_age_myr) {
  if (rad->is_2d) {
    return radiation_get_ionization_rate_from_raw_2d(rad, log_z, log_m,
                                                     star_age_myr);
  }
  return radiation_get_ionization_rate_from_raw(rad, log_m);
}

/**
 * @brief Get a single star's mean excess photon energy above the 13.6 eV HI
 * ionization threshold at a given mass, dispatching on #rad->is_2d between
 * the 1D (mass-only) and 2D (mass x metallicity) raw tables.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @param star_age_myr The star's current age, in Myr, ZAMS-anchored (see
 * #radiation_get_ionization_rate_from_raw_2d's own doxygen); used only if
 * #rad holds a 2D table.
 * @return Mean excess photon energy in cgs erg.
 */
double radiation_get_star_mean_excess_photon_energy_HI(
    const struct radiation *rad, float log_m, float log_z, float star_age_myr) {
  if (rad->is_2d) {
    return radiation_get_mean_excess_photon_energy_HI_from_raw_2d(
        rad, log_z, log_m, star_age_myr);
  }
  return radiation_get_mean_excess_photon_energy_HI_from_raw(rad, log_m);
}
