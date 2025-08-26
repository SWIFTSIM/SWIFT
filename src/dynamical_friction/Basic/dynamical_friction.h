/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DYNAMICAL_FRICTION_BASIC_H
#define SWIFT_DYNAMICAL_FRICTION_BASIC_H

#include "cosmology.h"
#include "error.h"
#include "feedback_properties.h"
#include "hydro_properties.h"
#include "part.h"
#include "units.h"

/**
 * @brief Prepares a s-particle for its df interactions from dm
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void df_from_dm_init_spart(
    struct spart* sp) {

  sp->df_data.density_dm.wcount = 0.f;
  sp->df_data.density_dm.wcount_dh = 0.f;

}

/**
 * @brief Prepares a s-particle for its df interactions from stars
 *
 * @param sp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void df_from_stars_init_spart(
    struct spart* sp) {

  sp->df_data.density_stars.wcount = 0.f;
  sp->df_data.density_stars.wcount_dh = 0.f;

}


/**
 * @brief Initialises the s-particles df props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in
 *
 * @param sp The particle to act upon.
 * @param df_props The properties of the df model.
 */
__attribute__((always_inline)) INLINE static void df_from_dm_first_init_spart(
    struct spart* sp, const struct stars_props* props) {

  /* Set the initial DF smoothing to be equal to the hydro smoothing just so we have something */
  sp->df_data.h_dm = sp->h;

  df_from_dm_init_spart(sp);
}

/**
 * @brief Initialises the s-particles df props for the first time
 *
 * This function is called only once just after the ICs have been
 * read in
 *
 * @param sp The particle to act upon.
 * @param df_props The properties of the df model.
 */
__attribute__((always_inline)) INLINE static void df_from_stars_first_init_spart(
    struct spart* sp, const struct stars_props* props) {

  /* Set the initial DF smoothing to be equal to the hydro smoothing just so we have something */
  sp->df_data.h_stars = sp->h;

  df_from_stars_init_spart(sp);
}


/**
 * @brief Finishes the calculation of density from DM on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_df_from_dm(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->df_data.h_dm;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->df_data.density_dm.wcount *= h_inv_dim;
  sp->df_data.density_dm.wcount_dh *= h_inv_dim_plus_one;
}


/**
 * @brief Finishes the calculation of density from stars on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_df_from_stars(
    struct spart* sp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = sp->df_data.h_stars;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->df_data.density_stars.wcount *= h_inv_dim;
  sp->df_data.density_stars.wcount_dh *= h_inv_dim_plus_one;
}


#endif /* SWIFT_DYNAMICAL_FRICTION_BASIC_H */
