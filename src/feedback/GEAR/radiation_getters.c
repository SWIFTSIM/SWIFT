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

#include <math.h>

/**
 * @brief Abort if #rad's table dimensionality does not match what the calling
 * getter expects.
 *
 * For an expected 2D table, #is_active is checked first: a #radiation with no
 * table also has is_2d = 0, and reporting it as 1D would mislead.
 *
 * @param rad The #radiation model.
 * @param expect_2d Nonzero if the caller needs a 2D table, 0 for a 1D table.
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
 * Z=0 (pristine gas) must not reach log10(). Out-of-range Z is clamped to the
 * lowest row as a single cell, with no mass-axis blending, so stars at or
 * below the lowest tabulated Z get a small discontinuity.
 *
 * @param Z Metallicity mass fraction (may be exactly 0).
 * @return log10(max(Z, #RADIATION_LOG_FLOOR_CGS)).
 */
float radiation_get_log_metallicity(float Z) {
  /* Floor in double: in float, 1e-300 would narrow to 0.0f and defeat it. */
  const double Z_floored = max((double)Z, RADIATION_LOG_FLOOR_CGS);
  return (float)log10(Z_floored);
}

/**
 * @brief Get the IMF-averaged bolometric luminosity per mass.
 *
 * Reads #rad->integrated.luminosities in linear value space, with no
 * exponentiation, unlike the raw table. See radiation_build_tables().
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return The bolometric luminosity.
 */
float radiation_get_luminosities_from_integral(const struct radiation *rad,
                                               float log_m1, float log_m2) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  float luminosity_1 = interpolate_1d(&rad->integrated.luminosities, log_m1);
  float luminosity_2 = interpolate_1d(&rad->integrated.luminosities, log_m2);
  return luminosity_2 - luminosity_1;
}

/**
 * @brief Get the non-IMF-integrated bolometric luminosity at a given mass.
 *
 * #rad->raw.luminosities holds log10(luminosity), which is exponentiated back
 * here. Luminosity is never 0, so no underflow case arises.
 *
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return The bolometric luminosity, internal units.
 */
float radiation_get_luminosities_from_raw(const struct radiation *rad,
                                          float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.luminosities, log_m));
}

/**
 * @brief Get the IMF-averaged ionization rate per mass.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return The ionization rate.
 */
double radiation_get_ionization_rate_from_integral(const struct radiation *rad,
                                                   float log_m1, float log_m2) {

  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  double dot_N_ion_1 = interpolate_1d(&rad->integrated.dot_N_ion, log_m1) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  double dot_N_ion_2 = interpolate_1d(&rad->integrated.dot_N_ion, log_m2) *
                       RADIATION_DOT_N_ION_TABLE_SCALING;
  return dot_N_ion_2 - dot_N_ion_1;
}

/**
 * @brief Get the non-IMF-integrated ionization rate at a given mass.
 *
 * #rad->raw.dot_N_ion holds log10(dot_N_ion /
 * #RADIATION_DOT_N_ION_TABLE_SCALING). The narrowing to float is load-bearing:
 * below the ionization threshold the value underflows float32 and returns
 * exactly 0.0f, which depends on the unit system. See
 * radiation_read_cgs_array() for the floor value.
 *
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
}

/**
 * @brief Get the IMF-averaged, Q-weighted mean excess photon energy above the
 * 13.6 eV HI threshold, for a population over a mass window.
 *
 * It is the ratio of the integrated dot_E_excess and dot_N_ion tables, taken on
 * the scaled values: the common scaling cancels, leaving cgs erg.
 *
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

  /* The cumulative table is non-decreasing, so a non-positive difference is
     degenerate (no ionizing stars, zero-width window or roundoff). Test <= 0,
     not == 0, which a near-cancellation could slip past. */
  if (delta_dot_N_ion <= 0.) return 0.;

  const double dot_E_excess_1 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m1);
  const double dot_E_excess_2 =
      interpolate_1d(&rad->integrated.dot_E_excess, log_m2);
  const double delta_dot_E_excess = dot_E_excess_2 - dot_E_excess_1;

  return delta_dot_E_excess / delta_dot_N_ion;
}

/**
 * @brief Get the non-IMF-integrated mean excess photon energy above the 13.6
 * eV HI threshold, for a single star of a given mass.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_integral on the raw
 * tables. Each table is exponentiated and narrowed to float before the ratio,
 * as in #radiation_get_ionization_rate_from_raw, so the `dot_N_ion <= 0.`
 * guard fires below the ionization threshold.
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
}

/**
 * @brief Get the non-IMF-integrated bolometric luminosity at a given mass and
 * metallicity, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_luminosities_from_raw. It is not capped by the
 * main-sequence lifetime: a post-main-sequence star stays luminous, so
 * truncating radiation pressure at TAMS would be worse than over-extending
 * the averaged L_bol.
 *
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
}

/**
 * @brief Return whether a star has evolved past MainSequenceLifetime(Z, M) in
 * the 2D ("M,Z") table. Used to cap Q_H and DotEExcess to 0.
 *
 *
 * @param rad The #radiation model (must hold an active 2D table).
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's age in Myr, ZAMS-anchored.
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
 * Mirrors #radiation_get_ionization_rate_from_raw. Past the table's
 * MainSequenceLifetime(Z, M), the time-averaging window of its Q_H, it returns
 * exactly 0. This cap is independent of Poirier lifetimes and
 * #feedback_properties.HII_max_age, and can fire earlier.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's age in Myr, anchored at the ZAMS, as in
 * MainSequenceLifetime. A caller using SWIFT's particle age (from spawn, with
 * any pre-main-sequence phase) must check the normalization matches.
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
}

/**
 * @brief Get the non-IMF-integrated mean excess photon energy above the 13.6
 * eV HI threshold, at a given mass, metallicity and age, from a 2D ("M,Z")
 * table.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_raw, with the same
 * MainSequenceLifetime cap as #radiation_get_ionization_rate_from_raw_2d,
 * applied first.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @param star_age_myr The star's age in Myr, ZAMS-anchored.
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
}

/**
 * @brief Nudge a log-mass query strictly inside a 2D IMF-integrated table's
 * top edge.
 *
 * A query exactly at the top edge takes #interpolate_2d's out-of-range branch,
 * which snaps the Z axis to the nearest row instead of blending. That case is
 * common: the upper mass bound is clamped to sm->imf.mass_max for an early-age
 * population. The nudge is a #RADIATION_2D_EDGE_EPS fraction of the mass span.
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
 * Mirrors #radiation_get_luminosities_from_integral. Below the table's native
 * mass floor the integrated table is flat at 0, so the subtraction stays
 * exact. The top-edge query is nudged by #radiation_nudge_mass_edge_2d.
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
}

/**
 * @brief Get the IMF-averaged ionization rate per mass, at a given
 * metallicity, from a 2D ("M,Z") table. See
 * #radiation_get_luminosities_from_integral_2d for the 2D caveats.
 *
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
}

/**
 * @brief Get the IMF-averaged, Q-weighted mean excess photon energy above the
 * 13.6 eV HI threshold, at a given metallicity, from a 2D ("M,Z") table.
 *
 * Mirrors #radiation_get_mean_excess_photon_energy_HI_from_integral, including
 * its degenerate-ratio guard. See #radiation_get_luminosities_from_integral_2d
 * for the shared 2D caveats.
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
}

/**
 * @brief Get the upper mass bound still on the main sequence, for a population
 * of a given age and metallicity, from MainSequenceLifetimeInverse.
 *
 * The threshold is the min() of the longest tabulated lifetimes of the two
 * bracketing metallicity rows (not a blend, which could admit a query one row
 * already excluded). It also gates on #radiation.age_max_myr.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param star_age_myr The population's age in Myr.
 * @param m_min Returned when no main-sequence mass remains (the caller's floor
 * mass).
 * @return The mass whose main-sequence lifetime equals @p star_age_myr at
 * @p log_z, in Msun, or @p m_min if none remains.
 */
float radiation_get_main_sequence_lifetime_inverse_mass_2d(
    const struct radiation *rad, float log_z, float star_age_myr, float m_min) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);

  int z_lo, z_hi;
  interpolate_2d_bracket_x(&rad->raw.main_sequence_lifetime_inverse_2d, log_z,
                           &z_lo, &z_hi);

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
 * dispatching on #rad->is_2d.
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
 * @brief Get a single star's ionization rate at a given mass, dispatching on
 * #rad->is_2d.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @param star_age_myr The star's age in Myr, ZAMS-anchored; used only for 2D.
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
 * threshold, dispatching on #rad->is_2d.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @param star_age_myr The star's age in Myr, ZAMS-anchored; used only for 2D.
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

/**
 * @brief Get a single star's photon-number-weighted mean LW photon energy at a
 * given mass, from a 1D table. A diagnostic: no rate reads it.
 *
 * Below pychem's LW mass floor it returns the band midpoint, 12.4 eV in erg, a
 * placeholder (see #radiation.raw.mean_photon_energy_lw).
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return Mean LW photon energy, cgs erg (NOT internal units, and NOT eV).
 */
double radiation_get_mean_photon_energy_lw_from_raw(const struct radiation *rad,
                                                    float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return exp10(interpolate_1d(&rad->raw.mean_photon_energy_lw, log_m));
}

/**
 * @brief Get a single star's mean LW photon energy at a given mass and
 * metallicity, from a 2D table. See
 * #radiation_get_mean_photon_energy_lw_from_raw.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @return Mean LW photon energy, cgs erg.
 */
double radiation_get_mean_photon_energy_lw_from_raw_2d(
    const struct radiation *rad, float log_z, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return exp10(
      interpolate_2d(&rad->raw.mean_photon_energy_lw_2d, log_z, log_m));
}

/**
 * @brief Get a single star's mean LW photon energy, dispatching on #rad->is_2d.
 *
 * Not capped by the main-sequence lifetime: it is a spectral shape, not a
 * rate. The caller gates on the LW luminosity instead.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if @p rad holds a 2D table.
 * @return Mean LW photon energy, cgs erg.
 */
double radiation_get_star_mean_photon_energy_lw(const struct radiation *rad,
                                                float log_m, float log_z) {
  if (rad->is_2d) {
    return radiation_get_mean_photon_energy_lw_from_raw_2d(rad, log_z, log_m);
  }
  return radiation_get_mean_photon_energy_lw_from_raw(rad, log_m);
}

/**
 * @brief Get the photon-number-weighted mean LW photon energy of the population
 * formed between the IMF's mass_min and @p log_m, from a 1D table.
 *
 * It takes one mass bound: the dataset is the ratio Integrated_L_LW /
 * Integrated_Q_LW, and pychem does not export Integrated_Q_LW, so a window mean
 * is not recoverable. The result is intensive: do not rescale it by the star's
 * birth mass.
 *
 * @param rad The #radiation model.
 * @param log_m Upper mass bound of the population, in log.
 * @return Mean LW photon energy, cgs erg.
 */
double radiation_get_mean_photon_energy_lw_from_integral(
    const struct radiation *rad, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return exp10(interpolate_1d(&rad->integrated.mean_photon_energy_lw, log_m));
}

/**
 * @brief Get the mean LW photon energy of the population formed between the
 * IMF's mass_min and @p log_m, at a given metallicity, from a 2D table. See
 * #radiation_get_mean_photon_energy_lw_from_integral.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m Upper mass bound of the population, in log.
 * @return Mean LW photon energy, cgs erg.
 */
double radiation_get_mean_photon_energy_lw_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return exp10(
      interpolate_2d(&rad->integrated.mean_photon_energy_lw_2d, log_z, log_m));
}

/**
 * @brief Get the photospheric effective temperature at a given mass, from a
 * 1D (mass-only) table.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return Effective temperature, internal units.
 */
float radiation_get_teff_from_raw(const struct radiation *rad, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.teff, log_m));
}

/**
 * @brief Get the photospheric effective temperature at a given mass and
 * metallicity, from a 2D ("M,Z") table.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @return Effective temperature, internal units.
 */
float radiation_get_teff_from_raw_2d(const struct radiation *rad, float log_z,
                                     float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return (float)exp10(interpolate_2d(&rad->raw.teff_2d, log_z, log_m));
}

/**
 * @brief Get a single star's effective temperature at a given mass,
 * dispatching on #rad->is_2d. Not capped by the main-sequence lifetime, like
 * #radiation_get_star_luminosity. Valid only when #radiation.has_teff is set;
 * callers must check.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity),
 * used only if #rad holds a 2D table.
 * @return Effective temperature, internal units.
 */
float radiation_get_star_teff(const struct radiation *rad, float log_m,
                              float log_z) {
  if (rad->is_2d) {
    return radiation_get_teff_from_raw_2d(rad, log_z, log_m);
  }
  return radiation_get_teff_from_raw(rad, log_m);
}

/**
 * @brief Get the non-IMF-integrated PE band emission rate at a given mass,
 * from a 1D table. See #radiation_get_luminosities_from_raw.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return PE band emission rate, internal units.
 */
float radiation_get_luminosity_pe_from_raw(const struct radiation *rad,
                                           float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.l_pe, log_m));
}

/**
 * @brief Get the non-IMF-integrated PE band emission rate at a given mass and
 * metallicity, from a 2D table. See #radiation_get_luminosities_from_raw_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @return PE band emission rate, internal units.
 */
float radiation_get_luminosity_pe_from_raw_2d(const struct radiation *rad,
                                              float log_z, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return (float)exp10(interpolate_2d(&rad->raw.l_pe_2d, log_z, log_m));
}

/**
 * @brief Get a single star's PE band emission rate at a given mass,
 * dispatching on #rad->is_2d. Valid only when #radiation.with_ISRF is set;
 * callers must check.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10, used only if #rad holds a 2D
 * table.
 * @return PE band emission rate, internal units.
 */
float radiation_get_star_luminosity_pe(const struct radiation *rad, float log_m,
                                       float log_z) {
  if (rad->is_2d) {
    return radiation_get_luminosity_pe_from_raw_2d(rad, log_z, log_m);
  }
  return radiation_get_luminosity_pe_from_raw(rad, log_m);
}

/**
 * @brief Get the non-IMF-integrated LW band emission rate at a given mass, from
 * a 1D table. See #radiation_get_luminosity_pe_from_raw.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @return Lyman-Werner band emission rate, internal units.
 */
float radiation_get_luminosity_lw_from_raw(const struct radiation *rad,
                                           float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  return (float)exp10(interpolate_1d(&rad->raw.l_lw, log_m));
}

/**
 * @brief Get the non-IMF-integrated LW band emission rate at a given mass and
 * metallicity, from a 2D table. See #radiation_get_luminosity_pe_from_raw_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m The mass in log.
 * @return Lyman-Werner band emission rate, internal units.
 */
float radiation_get_luminosity_lw_from_raw_2d(const struct radiation *rad,
                                              float log_z, float log_m) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  return (float)exp10(interpolate_2d(&rad->raw.l_lw_2d, log_z, log_m));
}

/**
 * @brief Get a single star's LW band emission rate at a given mass. See
 * #radiation_get_star_luminosity_pe.
 *
 * @param rad The #radiation model.
 * @param log_m The mass in log.
 * @param log_z The metallicity in log10, used only if #rad holds a 2D
 * table.
 * @return Lyman-Werner band emission rate, internal units.
 */
float radiation_get_star_luminosity_lw(const struct radiation *rad, float log_m,
                                       float log_z) {
  if (rad->is_2d) {
    return radiation_get_luminosity_lw_from_raw_2d(rad, log_z, log_m);
  }
  return radiation_get_luminosity_lw_from_raw(rad, log_m);
}

/**
 * @brief Get the IMF-averaged PE band emission rate per mass, from a 1D table.
 * Valid only when #radiation.with_ISRF is set.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return PE band emission rate per Msun of stars formed, internal units.
 */
float radiation_get_luminosity_pe_from_integral(const struct radiation *rad,
                                                float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const float l_pe_1 = interpolate_1d(&rad->integrated.l_pe, log_m1);
  const float l_pe_2 = interpolate_1d(&rad->integrated.l_pe, log_m2);
  return l_pe_2 - l_pe_1;
}

/**
 * @brief Get the IMF-averaged PE band emission rate per mass, at a given
 * metallicity, from a 2D table. See
 * #radiation_get_luminosities_from_integral_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return PE band emission rate per Msun of stars formed, internal units.
 */
float radiation_get_luminosity_pe_from_integral_2d(const struct radiation *rad,
                                                   float log_z, float log_m1,
                                                   float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.l_pe_2d;
  const float l_pe_1 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m1));
  const float l_pe_2 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m2));
  return l_pe_2 - l_pe_1;
}

/**
 * @brief Get the IMF-averaged LW band emission rate per mass, from a 1D table.
 * See #radiation_get_luminosity_pe_from_integral.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return Lyman-Werner band emission rate per Msun of stars formed,
 * internal units.
 */
float radiation_get_luminosity_lw_from_integral(const struct radiation *rad,
                                                float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const float l_lw_1 = interpolate_1d(&rad->integrated.l_lw, log_m1);
  const float l_lw_2 = interpolate_1d(&rad->integrated.l_lw, log_m2);
  return l_lw_2 - l_lw_1;
}

/**
 * @brief Get the IMF-averaged LW band emission rate per mass, at a given
 * metallicity, from a 2D table. See
 * #radiation_get_luminosity_pe_from_integral_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return Lyman-Werner band emission rate per Msun of stars formed,
 * internal units.
 */
float radiation_get_luminosity_lw_from_integral_2d(const struct radiation *rad,
                                                   float log_z, float log_m1,
                                                   float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.l_lw_2d;
  const float l_lw_1 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m1));
  const float l_lw_2 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m2));
  return l_lw_2 - l_lw_1;
}

/**
 * @brief Get the IMF-averaged PE band lower-edge spectral photon rate (scaled
 * by E_lo^2 and unit-converted like #l_pe), from a 1D table.
 *
 * It uses the same difference pattern as
 * #radiation_get_luminosity_pe_from_integral, so lambda_E(PE) - 1 never mixes
 * a windowed numerator with a whole-population denominator.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return E_lo(PE)^2 * (IMF-averaged dQ/dE at the PE edge), internal power
 * units.
 */
float radiation_get_luminosity_edge_pe_from_integral(
    const struct radiation *rad, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const float l_edge_pe_1 = interpolate_1d(&rad->integrated.l_edge_pe, log_m1);
  const float l_edge_pe_2 = interpolate_1d(&rad->integrated.l_edge_pe, log_m2);
  return l_edge_pe_2 - l_edge_pe_1;
}

/**
 * @brief Get the IMF-averaged PE band lower-edge spectral photon rate, at a
 * given metallicity, from a 2D table. See
 * #radiation_get_luminosity_pe_from_integral_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return E_lo(PE)^2 * (IMF-averaged dQ/dE at the PE edge), internal power
 * units.
 */
float radiation_get_luminosity_edge_pe_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.l_edge_pe_2d;
  const float l_edge_pe_1 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m1));
  const float l_edge_pe_2 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m2));
  return l_edge_pe_2 - l_edge_pe_1;
}

/**
 * @brief Get the IMF-averaged LW band lower-edge spectral photon rate, from a
 * 1D table. See #radiation_get_luminosity_edge_pe_from_integral.
 *
 * @param rad The #radiation model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return E_lo(LW)^2 * (IMF-averaged dQ/dE at the LW edge), internal power
 * units.
 */
float radiation_get_luminosity_edge_lw_from_integral(
    const struct radiation *rad, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/0, __func__);
  const float l_edge_lw_1 = interpolate_1d(&rad->integrated.l_edge_lw, log_m1);
  const float l_edge_lw_2 = interpolate_1d(&rad->integrated.l_edge_lw, log_m2);
  return l_edge_lw_2 - l_edge_lw_1;
}

/**
 * @brief Get the IMF-averaged LW band lower-edge spectral photon rate, at a
 * given metallicity, from a 2D table. See
 * #radiation_get_luminosity_edge_pe_from_integral_2d.
 *
 * @param rad The #radiation model.
 * @param log_z The metallicity in log10 (see #radiation_get_log_metallicity).
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @return E_lo(LW)^2 * (IMF-averaged dQ/dE at the LW edge), internal power
 * units.
 */
float radiation_get_luminosity_edge_lw_from_integral_2d(
    const struct radiation *rad, float log_z, float log_m1, float log_m2) {
  radiation_check_dimensionality(rad, /*expect_2d=*/1, __func__);
  const struct interpolation_2d *interp = &rad->integrated.l_edge_lw_2d;
  const float l_edge_lw_1 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m1));
  const float l_edge_lw_2 = interpolate_2d(
      interp, log_z, radiation_nudge_mass_edge_2d(interp, log_m2));
  return l_edge_lw_2 - l_edge_lw_1;
}
