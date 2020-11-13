/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#endif /* SWIFT_EAGLE_STARS_H */
