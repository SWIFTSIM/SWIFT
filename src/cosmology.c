/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
 *  @file cosmology.c
 *  @brief Functions relating cosmological parameters
 */

/* This object's header. */
#include "cosmology.h"

/* Some standard headers */
#include <math.h>

/**
 * @brief Update the cosmological parameters to the current simulation time.
 *
 * @param c The #cosmology struct.
 * @param e The #engine containing information about the simulation time.
 */
void cosmology_update(struct cosmology *c, const struct engine *e) {

  /* Get scale factor */
  const double a = c->a_begin * exp(e->ti_current * e->timeBase);
  const double a_inv = 1. / a;
  c->a = a;
  c->a_inv = a_inv;

  /* Redshift */
  c->z  = a_inv - 1.;

  /* E(z) */
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_l = c->Omega_lambda;
  const double Omega_k = c->Omega_k;

  const double E_z2 = Omega_r * a_inv * a_inv * a_inv * a_inv
                    + Omega_m * a_inv * a_inv * a_inv
                    + Omega_k * a_inv * a_inv
                    + Omega_l;

  const double E_z = sqrt(E_z2);

  /* H(z) */
  c->H = c->H0 * E_z;
}


void cosmology_init(const struct swift_params* params,
		    const struct unit_system* us,
		    const struct phys_const* phys_const,
		    struct cosmology *c) {

  /* Read in the cosmological parameters */
  c->Omega_m = parser_get_param_double(params, "Cosmology:Omega_m");
  c->Omega_r = parser_get_param_double(params, "Cosmology:Omega_r");
  c->Omega_lambda = parser_get_param_double(params, "Cosmology:Omega_lambda");
  c->Omega_b = parser_get_param_double(params, "Cosmology:Omega_b");
  c->h = parser_get_param_double(params, "Cosmology:h");
  c->a_begin = parser_get_param_double(params, "Cosmology:a_begin");
  c->a_end = parser_get_param_double(params, "Cosmology:a_end");

  /* Construct derived quantities */
  c->Omega_k = 1. - (c->Omega_m + c->Omega_r + c->Omega_lambda);
  const double km = 1.e5 / units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  const double s = 1. / units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  message("km=%e", km);
  const double H0_cgs = 100. * c->h * km / s / (1.e6 * phys_const->const_parsec); /* s^-1 */
  message("s=%e", s);
  c->H0 = H0_cgs;

  /* Set remaining variables to invalid values */
  c->H = -1.;
  c->a = -1.;
  c->a_inv = -1;
  c->z = -1.;
}


void cosmology_print(const struct cosmology *c) {

  message("Density parameters: [O_m, O_l, O_b, O_k, O_r] = [%f, %f, %f, %f, %f]",
	  c->Omega_m, c->Omega_lambda, c->Omega_b, c->Omega_k, c->Omega_r);
  message("Hubble constant: h = %f, H_0 = %e U_t^(-1) (internal units)",
	  c->h, c->H0);

}


