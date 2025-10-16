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
    struct sipart* sip, const struct sidm_props* sidm_properties) {
  sip->time_bin = 0;
}

/**
 * @brief Prepares an si-particle for its interactions
 *
 * @param sip The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_sipart(
    struct sipart* sip) {

  /* TODO - add SIDM debugging */
  // #ifdef DEBUG_INTERACTIONS_SIDM
  //   for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_SIDM; ++i)
  //     sip->ids_ngbs_density[i] = -1;
  //   sip->num_ngb_density = 0;
  // #endif

  sip->density.wcount = 0.f;
  sip->density.wcount_dh = 0.f;
}

/**
 * @brief Finishes the calculation of density on SIDM
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_end_density(
    struct sipart* sip, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sip->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sip->density.wcount *= h_inv_dim;
  sip->density.wcount_dh *= h_inv_dim_plus_one;
}

/**
 * @brief Sets all particle fields to sensible values when the #sipart has 0
 * ngbs.
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_sipart_has_no_neighbours(
    struct sipart* restrict sip, const struct cosmology* cosmo) {

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

#endif /* SWIFT_BASIC_SIDM_H */
