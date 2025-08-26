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

  sp->df_data.dm_df_a[0] = 0.f;
  sp->df_data.dm_df_a[1] = 0.f;
  sp->df_data.dm_df_a[2] = 0.f;

  sp->df_data.dm_v_medium[0] = 0.f;
  sp->df_data.dm_v_medium[1] = 0.f;
  sp->df_data.dm_v_medium[2] = 0.f;

  sp->df_data.dm_ngb_mass = 0.f;
  sp->df_data.dm_rho = 0.f;
  sp->df_data.dm_sigma2 = 0.f;

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
    struct spart* sp, const struct cosmology* cosmo, const struct gravity_props* gprops, 
    const struct phys_const* phys_const) {

  /* Some smoothing length multiples. */
  const float h = sp->df_data.h_dm;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* Finish the calculation by inserting the missing h-factors */
  sp->df_data.density_dm.wcount *= h_inv_dim;
  sp->df_data.density_dm.wcount_dh *= h_inv_dim_plus_one;

  /* Finish the density calculation */
  sp->df_data.dm_rho *= h_inv_dim;

  /* Finish the velocity calculation */
  sp->df_data.dm_v_medium[0] *= h_inv_dim / sp->df_data.dm_rho;
  sp->df_data.dm_v_medium[1] *= h_inv_dim / sp->df_data.dm_rho;
  sp->df_data.dm_v_medium[2] *= h_inv_dim / sp->df_data.dm_rho;
    
  /* Finish the velocity dispersion calculation */
  sp->df_data.dm_sigma2 *= h_inv_dim / sp->df_data.dm_rho;

  float gc_v2 = sp->df_data.dm_v_medium[0] * sp->df_data.dm_v_medium[0] + 
                sp->df_data.dm_v_medium[1] * sp->df_data.dm_v_medium[1] + 
                sp->df_data.dm_v_medium[2] * sp->df_data.dm_v_medium[2];

  /* Calculate (actual) gas velocity dispersion. Currently, the variable
  * 'sigma2' holds <v^2> instead. */
  sp->df_data.dm_sigma2 -= gc_v2;

  /* Invert the sign to get the velocity of the GC relative to the medium */
  float gc_v[3];
  gc_v[0] = -1. * sp->df_data.dm_v_medium[0];
  gc_v[1] = -1. * sp->df_data.dm_v_medium[1];
  gc_v[2] = -1. * sp->df_data.dm_v_medium[2];

  float gc_v_mag = sqrt(gc_v2);

  /* Gather some physical constants (all in internal units) */
  const double G = phys_const->const_newton_G;

  /* Use a fixed value of eta for now. Use 6 as we'll start with the Chen prescription */
  /* This will probably be too strong as we're also using the velocity dispersion term in b_min from Anna */

  const float df_eta = 6.f;

  /* Coulomb logarithm */
  const float b_max = df_eta * gravity_get_softening(sp->gpart, gprops);

  /* Deflection radius, with Genina+24 velocity dispersion term */
  float b_90 = G * sp->mass / (gc_v2 + ((2./3.) * sp->df_data.dm_sigma2));

  // Gc stuff for b_min - ignore
  // /* Twice the fixed GC stellar half-mass radius (rather than the Schwarzchild radius as used for BHs) */
  // const float r_gc = 2. * sp->rh[0]; // Standalone GC information is always in the 0th element
  // float b_min = max(b_90, r_gc);

  /* Make sure b_min is never larger than b_max - this can happen if the velocity is very low */
  /* When b_min >= b_max, we are in the regime where DF is fully resolved, so the correction is zero. */
  float b_min = min(b_90,b_max);

  double log_lambda = log(b_max/b_min);

  /* From Anna's method - save for later */
  // double rho_slow = ngb_rho * ngb_mass_slow / ngb_mass_total;
  // double prefactor = -4.f * M_PI * G * G * sp->mass / (gc_v2 * gc_v_mag); // -4pi G^2 M/ v^3

  /* MAXWELLIAN APPROXIMATION - does not require a second loop over neighbours to find slow particles. */

  const float gc_v_ov_sigma = gc_v_mag/sqrt(sp->df_data.dm_sigma2);
  double f_maxwell_integral = erf(gc_v_ov_sigma) - (2. * gc_v_ov_sigma * exp(-1. * gc_v_ov_sigma * gc_v_ov_sigma) / sqrt(M_PI));
  double prefactor = -4.f * M_PI * G * G * sp->mass / (gc_v2 * gc_v_mag); // -4pi G^2 M / v^3

  sp->df_data.dm_df_a[0] = prefactor * sp->df_data.dm_rho * f_maxwell_integral * log_lambda * gc_v[0];
  sp->df_data.dm_df_a[1] = prefactor * sp->df_data.dm_rho * f_maxwell_integral * log_lambda * gc_v[1];
  sp->df_data.dm_df_a[2] = prefactor * sp->df_data.dm_rho * f_maxwell_integral * log_lambda * gc_v[2];

  message("DM DF: star pid=%lld, h=%e, ngb_mass=%e, ngb_rho=%e, prefactor=%e, f_maxwell=%e, log_lambda=%e, b_max=%e, b_min=%e, b_90=%e, v^2=%e sigma^2=%e, gc_v_rel=[%f, %f, %f], gc_v_abs=[%f, %f, %f], df_a=[%f, %f, %f]",
            sp->id, h, sp->df_data.dm_ngb_mass, sp->df_data.dm_rho, prefactor, f_maxwell_integral, log_lambda, b_max, b_min, b_90, gc_v2, sp->df_data.dm_sigma2, gc_v[0], gc_v[1], gc_v[2], sp->v[0],sp->v[1],sp->v[2], sp->df_data.dm_df_a[0],sp->df_data.dm_df_a[1],sp->df_data.dm_df_a[2]);



}


/**
 * @brief Finishes the calculation of density from stars on stars
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void stars_end_df_from_stars(
    struct spart* sp, const struct cosmology* cosmo, const struct gravity_props* gprops, 
    const struct phys_const* phys_const) {

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
