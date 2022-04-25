/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STARS_H
#define SWIFT_EAGLE_STARS_H

#include "exp10.h"

#include <float.h>

/**
 * @brief Computes the time-step length of a given star particle from stars
 * physics
 *
 * @param sp Pointer to the s-particle data.
 * @param stars_properties Properties of the stars model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static float stars_compute_timestep(
    const struct spart* const sp, const struct stars_props* stars_properties,
    const int with_cosmology, const struct cosmology* cosmo,
    const double time) {

  /* Background star particles have no time-step limits */
  if (sp->birth_time == -1.) {
    return FLT_MAX;
  }

  /* Star age (in internal units) */
  double star_age;
  if (with_cosmology) {

    /* Deal with rounding issues */
    if (sp->birth_scale_factor >= cosmo->a) {
      star_age = 0.;
    } else {
      star_age = cosmology_get_delta_time_from_scale_factors(
          cosmo, sp->birth_scale_factor, cosmo->a);
    }
  } else {
    star_age = time - sp->birth_time;
  }

  /* What age category are we in? */
  if (star_age > stars_properties->age_threshold_unlimited) {
    return FLT_MAX;
  } else if (star_age > stars_properties->age_threshold) {
    return stars_properties->max_time_step_old;
  } else {
    return stars_properties->max_time_step_young;
  }
}

/**
 * @brief Returns the age of a star in internal units
 *
 * @param sp The star particle.
 * @param cosmo The cosmological model.
 * @param time The current time (in internal units).
 * @param with_cosmology Are we running with cosmological integration?
 */
__attribute__((always_inline)) INLINE static double stars_compute_age(
    const struct spart* sp, const struct cosmology* cosmo, double time,
    const int with_cosmology) {

  if (with_cosmology) {
    return cosmology_get_delta_time_from_scale_factors(
        cosmo, (double)sp->birth_scale_factor, cosmo->a);
  } else {
    return time - (double)sp->birth_time;
  }
}

/**
 * @brief Prepares a s-particle for its interactions
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_init_spart(
    struct spart* sp) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    sp->ids_ngbs_density[i] = -1;
  sp->num_ngb_density = 0;
#endif

  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->N_density = 0;
  sp->N_density_exact = 0;
  sp->rho = 0.f;
  sp->rho_exact = 0.f;
  sp->n = 0.f;
  sp->n_exact = 0.f;
  sp->inhibited_exact = 0;
#endif
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param stars_properties Properties of the stars model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param scale_factor The current scale-factor (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp, const struct stars_props* stars_properties,
    const int with_cosmology, const double scale_factor, const double time) {

  sp->time_bin = 0;
  sp->f_E = -1.f;
  sp->count_since_last_enrichment = -1;
  sp->number_of_heating_events = 0.;
  sp->number_of_SNII_events = 0;

  if (stars_properties->overwrite_birth_time)
    sp->birth_time = stars_properties->spart_first_init_birth_time;
  if (stars_properties->overwrite_birth_density)
    sp->birth_density = stars_properties->spart_first_init_birth_density;
  if (stars_properties->overwrite_birth_temperature)
    sp->birth_temperature =
        stars_properties->spart_first_init_birth_temperature;

  if (with_cosmology)
    sp->last_enrichment_time = scale_factor;
  else
    sp->last_enrichment_time = time;

  stars_init_spart(sp);
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void stars_predict_extra(
    struct spart* restrict sp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sp The particle.
 */
__attribute__((always_inline)) INLINE static void stars_reset_predicted_values(
    struct spart* restrict sp) {}

/**
 * @brief Finishes the calculation of (non-gravity) forces acting on stars
 *
 * Multiplies the forces and accelerations by the appropiate constants
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_end_feedback(
    struct spart* sp) {}

/**
 * @brief Kick the additional variables
 *
 * @param sp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void stars_kick_extra(
    struct spart* sp, float dt) {}

/**
 * @brief Finishes the calculation of density on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_density(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->density.wcount *= h_inv_dim;
  sp->density.wcount_dh *= h_inv_dim_plus_one;

#ifdef SWIFT_STARS_DENSITY_CHECKS
  sp->rho *= h_inv_dim;
  sp->n *= h_inv_dim;
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_spart_has_no_neighbours(
    struct spart* restrict sp, const struct cosmology* cosmo) {

  warning(
      "Star particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      sp->id, sp->h, sp->density.wcount);

  /* Re-set problematic values */
  sp->density.wcount = 0.f;
  sp->density.wcount_dh = 0.f;
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_acceleration(
    struct spart* restrict p) {

#ifdef DEBUG_INTERACTIONS_STARS
  p->num_ngb_feedback = 0;
#endif
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on star, therefore no need to use it.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stars_reset_feedback(
    struct spart* restrict p) {

#ifdef DEBUG_INTERACTIONS_STARS
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    p->ids_ngbs_feedback[i] = -1;
  p->num_ngb_feedback = 0;
#endif
}

/**
 * @brief Compute the luminosities of a particles in different bands.
 *
 * @param sp The particle.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param cosmo The #cosmology object.
 * @param time The current physical time (internal units).
 * @param phys_const The physical constants in internal units.
 * @param props The #stars_props of that run.
 * @param luminosities (return) The luminosity in each band.
 */
INLINE static void stars_get_luminosities(
    const struct spart* sp, const int with_cosmology,
    const struct cosmology* cosmo, const double time,
    const struct phys_const* phys_const, const struct stars_props* props,
    float luminosities[luminosity_bands_count]) {

  const int count_Z = eagle_stars_lum_tables_N_Z;
  const int count_ages = eagle_stars_lum_tables_N_ages;

  /* Get star properties (all in internal units */
  const float Z =
      chemistry_get_star_total_metal_mass_fraction_for_luminosity(sp);
  const float mass = sp->mass_init;
  float age;
  if (with_cosmology)
    age = cosmology_get_delta_time_from_scale_factors(
        cosmo, sp->birth_scale_factor, cosmo->a);
  else
    age = time - sp->birth_time;

  /* Convert to the units of the tables */
  const float mass_Msun = mass / phys_const->const_solar_mass;
  const float age_Gyr = age / phys_const->const_year / 1e9;

  for (int i = 0; i < (int)luminosity_bands_count; ++i) {

    /* Log things */
    float log10_Z = log10(Z + FLT_MIN);
    float log10_age_Gyr = log10(age_Gyr + FLT_MIN);

    /* Clip the input */
    log10_Z = max(log10_Z, props->lum_tables_Z[i][0]);
    log10_Z = min(log10_Z, props->lum_tables_Z[i][count_Z - 1]);
    log10_age_Gyr = max(log10_age_Gyr, props->lum_tables_ages[i][0]);
    log10_age_Gyr =
        min(log10_age_Gyr, props->lum_tables_ages[i][count_ages - 1]);

    /* Get index along the interpolation axis */
    int Z_index = 0, age_index = 0;
    for (int j = 0; j < count_Z - 1; ++j) {
      if (log10_Z >= props->lum_tables_Z[i][j]) ++Z_index;
    }
    for (int j = 0; j < count_ages - 1; ++j) {
      if (log10_age_Gyr >= props->lum_tables_ages[i][j]) ++age_index;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (age_index == 0) error("Invalid age index!");
    if (Z_index == 0) error("Invalid Z index!");
#endif

    const float* array = props->lum_tables_luminosities[i];

    /* 2D interpolation */

    const float f_11 = array[(Z_index - 1) * count_ages + (age_index - 1)];
    const float f_12 = array[(Z_index - 0) * count_ages + (age_index - 1)];
    const float f_21 = array[(Z_index - 1) * count_ages + (age_index - 0)];
    const float f_22 = array[(Z_index - 0) * count_ages + (age_index - 0)];

    const float x_diff1 = props->lum_tables_ages[i][age_index] - log10_age_Gyr;
    const float x_diff2 =
        log10_age_Gyr - props->lum_tables_ages[i][age_index - 1];
    const float x_diff3 = props->lum_tables_ages[i][age_index] -
                          props->lum_tables_ages[i][age_index - 1];

    const float y_diff1 = props->lum_tables_Z[i][Z_index] - log10_Z;
    const float y_diff2 = log10_Z - props->lum_tables_Z[i][Z_index - 1];
    const float y_diff3 =
        props->lum_tables_Z[i][Z_index] - props->lum_tables_Z[i][Z_index - 1];

    const float f_1 = (f_11 * x_diff1 + f_21 * x_diff2) / x_diff3;
    const float f_2 = (f_12 * x_diff1 + f_22 * x_diff2) / x_diff3;

    const float log10_f = (y_diff1 * f_1 + y_diff2 * f_2) / y_diff3;

    /* Final conversion */
    luminosities[i] = exp10f(log10_f) * mass_Msun * props->lum_tables_factor;
  }
}

#endif /* SWIFT_EAGLE_STARS_H */
