/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_GEAR_STARS_H
#define SWIFT_GEAR_STARS_H

#include "minmax.h"

#include <float.h>

/**
 * @brief Computes the time-step length of a given star particle from star
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

  /* Note (01.08.2024):
     The following lines come from the EAGLE model. However, in GEAR, a star
     particle can now be a stellar population, a single star on a part of a
     population (continuous IMF).
     The two populations types are similar in behaviour and I think we can use
     the code below as-is. However, the single stars must be treated in a
     different way. They only have one SN feedback. Thus, they can be 'old' only
     after this explosion. The question is how to do that...
     Maybe we should flag the stars as being old or young? The discrete stars
     then will only become old after their SN feedback.

     Notice however that this require knowledge of the star_type, which is
     currently in the feedback module. I planned to move this from the feedback
     to the stars module (it makes more sense to be with the stars), but this
     require care because the sink module also need to know about this
     star_type. Moving the star type to the stars module simplifies the above
     since we would know the type of the star. Hence, no flag is needed.

     The purpose of taking care about the discrete stars is to avoid them
     having small timesteps when they are already dead. They won't do anything
     anymore so they don't need to be waken up often.
  */

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
__attribute__((always_inline)) INLINE static float stars_compute_age(
    const struct spart* sp, const struct cosmology* cosmo, double time,
    const int with_cosmology) {

  if (with_cosmology) {
    const double birth = sp->birth_scale_factor;
    return cosmology_get_delta_time_from_scale_factors(
        cosmo, min(birth, cosmo->a), cosmo->a);
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
}

/**
 * @brief Initialises the s-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sp The particle to act upon.
 * @param stars_properties Properties of the stars model.
 * @param with_cosmology Are we running a cosmological simulation?
 * @param scale_factor The current scale factor.
 * @param time The current time.
 */
__attribute__((always_inline)) INLINE static void stars_first_init_spart(
    struct spart* sp, const struct stars_props* stars_properties,
    const int with_cosmology, const double scale_factor, const double time) {

  sp->time_bin = 0;

  if (stars_properties->overwrite_birth_time)
    sp->birth_time = stars_properties->spart_first_init_birth_time;

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

#endif /* SWIFT_GEAR_STARS_H */
