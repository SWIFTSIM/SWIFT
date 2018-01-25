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

/* Local headers */
#include "inline.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#endif

const int cosmology_table_length = 10000;
const size_t GSL_workspace_size = 10000;

/**
 * @brief Compute \f$ E(z) \f$.
 */
static INLINE double E(double Or, double Om, double Ok, double Ol,
                       double a_inv) {

  return sqrt(Or * a_inv * a_inv * a_inv * a_inv + Om * a_inv * a_inv * a_inv +
              Ok * a_inv * a_inv + Ol);
}

double cosmology_get_time(const struct cosmology *c, double a) {

  const double log_a = log(a);

  /* Position in the table */
  const double x =
      ((log_a - c->log_a_begin) / (c->log_a_end - c->log_a_begin)) *
      cosmology_table_length;
  const int i =
      ((int)x >= cosmology_table_length) ? cosmology_table_length - 1 : (int)x;

  if (i <= 1)
    return c->time_interp_table_offset + x * c->time_interp_table[0];
  else
    return c->time_interp_table_offset + c->time_interp_table[i - 1] +
           (c->time_interp_table[i] - c->time_interp_table[i - 1]) * (x - i);
}

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
  c->z = a_inv - 1.;

  /* E(z) */
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, a_inv);

  /* H(z) */
  c->H = c->H0 * E_z;

  /* Time since Big Bang */
  c->time = cosmology_get_time(c, a);
}

void cosmology_init_tables(struct cosmology *c) {

  const ticks tic = getticks();

#ifdef HAVE_LIBGSL

  /* Retrieve some constants */
  const double a_begin = c->a_begin;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double H0 = c->H0;

  /* Allocate memory for the interpolation tables */
  c->drift_fac_interp_table = malloc(cosmology_table_length * sizeof(double));
  c->grav_kick_fac_interp_table =
      malloc(cosmology_table_length * sizeof(double));
  c->time_interp_table = malloc(cosmology_table_length * sizeof(double));

  /* Define the integrands we need */

  /* dt / a^2 */
  double drift_integrand(double a, void *param) {

    const double a_inv = 1. / a;
    const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, a_inv);
    const double H = H0 * E_z;

    return (1. / H) * a_inv * a_inv * a_inv;
  }

  /* dt / a */
  double gravity_kick_integrand(double a, void *param) {

    const double a_inv = 1. / a;
    const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, a_inv);
    const double H = H0 * E_z;

    return (1. / H) * a_inv * a_inv;
  }

  /* dt */
  double time_integrand(double a, void *param) {

    const double a_inv = 1. / a;
    const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, a_inv);
    const double H = H0 * E_z;

    return (1. / H) * a_inv;
  }

  /* Prepare a table of scale factors for the integral bounds */
  const double delta_a =
      (c->log_a_end - c->log_a_begin) / cosmology_table_length;
  double *a_table = malloc(cosmology_table_length * sizeof(double));
  for (int i = 0; i < cosmology_table_length; i++)
    a_table[i] = exp(c->log_a_begin + delta_a * (i + 1));

  /* Initalise the GSL workspace */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(GSL_workspace_size);

  double result, abserr;

  /* Integrate the drift factor \int_{a_begin}^{a_table[i]} dt/a^2 */
  gsl_function F = {&drift_integrand, NULL};
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->drift_fac_interp_table[i] = result;
  }

  /* Integrate the kick factor \int_{a_begin}^{a_table[i]} dt/a */
  F.function = &gravity_kick_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->grav_kick_fac_interp_table[i] = result;
  }

  /* Integrate the time \int_{a_begin}^{a_table[i]} dt */
  F.function = &time_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->time_interp_table[i] = result;
  }

  /* Integrate the time \int_{0}^{a_begin} dt */
  gsl_integration_qag(&F, 0., a_begin, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->time_interp_table_offset = result;

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);
  free(a_table);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif

  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());
}

void cosmology_init(const struct swift_params *params,
                    const struct unit_system *us,
                    const struct phys_const *phys_const, struct cosmology *c) {

  /* Read in the cosmological parameters */
  c->Omega_m = parser_get_param_double(params, "Cosmology:Omega_m");
  c->Omega_r = parser_get_param_double(params, "Cosmology:Omega_r");
  c->Omega_lambda = parser_get_param_double(params, "Cosmology:Omega_lambda");
  c->Omega_b = parser_get_param_double(params, "Cosmology:Omega_b");
  c->h = parser_get_param_double(params, "Cosmology:h");

  /* Read the start and end of the simulation */
  c->a_begin = parser_get_param_double(params, "Cosmology:a_begin");
  c->a_end = parser_get_param_double(params, "Cosmology:a_end");
  c->log_a_begin = log(c->a_begin);
  c->log_a_end = log(c->a_end);

  /* Construct derived quantities */

  /* Curvature density (for closure) */
  c->Omega_k = 1. - (c->Omega_m + c->Omega_r + c->Omega_lambda);

  /* Hubble constant in internal units */
  const double km = 1.e5 / units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  const double H0_cgs =
      100. * c->h * (km / (1.e6 * phys_const->const_parsec)); /* s^-1 */
  c->H0 = H0_cgs * units_cgs_conversion_factor(us, UNIT_CONV_TIME);

  /* Set remaining variables to invalid values */
  c->H = -1.;
  c->a = -1.;
  c->a_inv = -1;
  c->z = -1.;
  c->time = -1.;

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  cosmology_init_tables(c);
}

void cosmology_print(const struct cosmology *c) {

  message(
      "Density parameters: [O_m, O_l, O_b, O_k, O_r] = [%f, %f, %f, %f, %f]",
      c->Omega_m, c->Omega_lambda, c->Omega_b, c->Omega_k, c->Omega_r);
  message("Hubble constant: h = %f, H_0 = %e U_t^(-1) (internal units)", c->h,
          c->H0);
}
