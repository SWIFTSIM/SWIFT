/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_BASIC_SIDM_H
#define SWIFT_BASIC_SIDM_H

#include <float.h>

/* Local includes */
#include "sidm_part.h"
#include "sidm_properties.h"

/**
 * @brief Initialises the si-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sip The particle to act upon
 * @param sidm_properties The properties of the SIDM model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param scale_factor The current scale-factor (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void sidm_first_init_sipart(
    struct sipart *sip, const struct sidm_props *sidm_properties) {
  sip->time_bin = 0;
}

/**
 * @brief Prepares an si-particle for its interactions
 *
 * @param sip The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_sipart(
    struct sipart *sip) {

  /* TODO - add SIDM debugging */
  // #ifdef DEBUG_INTERACTIONS_SIDM
  //   for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SIDM; ++i)
  //     sip->ids_ngbs_density[i] = -1;
  //   sip->num_ngb_density = 0;
  // #endif

  sip->density.wcount = 0.f;
  sip->density.wcount_dh = 0.f;
  sip->rho = 0.f;
  sip->density.rho_dh = 0.f;

#ifdef SWIFT_SIDM_DENSITY_CHECKS
  sip->N_density = 0;
  sip->N_density_exact = 0;
  sip->rho_exact = 0.f;
  sip->n = 0.f;
  sip->n_exact = 0.f;
  sip->inhibited_exact = 0;
#endif
}

/**
 * @brief Finishes the calculation of density on SIDM
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_end_density(
    struct sipart *sip) {

  /* Some smoothing length multiples. */
  const float h = sip->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Final operation on the density (add self-contribution). */
  sip->rho += sip->mass * kernel_root;
  sip->density.rho_dh -= hydro_dimension * sip->mass * kernel_root;
  sip->density.wcount += kernel_root;
  sip->density.wcount_dh -= hydro_dimension * kernel_root;

  /* Finish the calculation by inserting the missing h-factors */
  sip->rho *= h_inv_dim;
  sip->density.rho_dh *= h_inv_dim_plus_one;
  sip->density.wcount *= h_inv_dim;
  sip->density.wcount_dh *= h_inv_dim_plus_one;

#ifdef SWIFT_SIDM_DENSITY_CHECKS
  sip->n += kernel_root;
  sip->n *= h_inv_dim;
#endif
}

/**
 * @brief Sets all particle fields to sensible values when the #sipart has 0
 * ngbs.
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_sipart_has_no_neighbours(
    struct sipart *restrict sip, const struct cosmology *cosmo) {

  warning(
      "SIDM particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      sip->id, sip->h, sip->density.wcount);

  /* Some smoothing length multiples. */
  const float h = sip->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  sip->density.wcount = kernel_root * h_inv_dim;
  sip->density.wcount_dh = 0.f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sip The #sipart.
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void sidm_predict_extra(
    struct sipart *restrict sip, float dt_drift) {

  /* TODO: We may need to add something here later like evolving rho */
}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param sip The si-particle of interest
 */
__attribute__((always_inline)) INLINE static float sidm_get_comoving_density(
    const struct sipart *restrict sip) {

  return sip->rho;
}

/**
 * @brief Kick the additional variables
 *
 * @param sip The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sidm_kick_extra(
    struct sipart *sip, float dt) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sip The particle.
 */
__attribute__((always_inline)) INLINE static void sidm_reset_predicted_values(
    struct sipart *restrict sip) {}

/**
 * @brief Computes the time-step of a given SIDM particle.
 *
 * @param sip Pointer to the si-particle data.
 * @param sidm_properties Properties of the SIDM model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 * @param time_base The time base.
 */
__attribute__((always_inline)) INLINE static float sidm_compute_timestep(
    const struct sipart *const sip, const struct sidm_props *sidm_properties,
    const int with_cosmology, const struct cosmology *cosmo, const double time,
    const double time_base) {

  return FLT_MAX;
}

/**
 * @brief Operations performed when an SIDM particle gets removed from the
 * simulation volume.
 *
 * @param sip The si-particle.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void sidm_remove_sipart(
    const struct sipart *sip, const double time) {}
#endif /* SWIFT_BASIC_SIDM_H */
