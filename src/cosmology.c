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
#include "adiabatic_index.h"
#include "inline.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#endif

/*! Number of values stored in the cosmological interpolation tables */
const int cosmology_table_length = 10000;

/*! Size of the GSL workspace */
const size_t GSL_workspace_size = 100000;

/**
 * @brief Returns the interpolated value from a table.
 *
 * Uses linear interpolation.
 *
 * @brief table The table of value to interpolate from (should be of length
 * cosmology_table_length).
 * @brief x The value to interpolate at.
 * @brief x_min The mininum of the range of x.
 * @brief x_max The maximum of the range of x.
 */
static INLINE double interp_table(double *table, double x, double x_min,
                                  double x_max) {

  const double xx = ((x - x_min) / (x_max - x_min)) * cosmology_table_length;
  const int i = (int)xx;
  const int ii = (i >= cosmology_table_length) ? cosmology_table_length - 1 : i;

  if (ii <= 1)
    return table[0] * xx;
  else
    return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
}

/**
 * @brief Computes the dark-energy equation of state at a given scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 *
 * @param a The current scale-factor
 * @param w_0 The equation of state parameter at z=0
 * @param w_a The equation of state evolution parameter
 */
static INLINE double cosmology_dark_energy_EoS(double a, double w_0,
                                               double w_a) {

  return w_0 + w_a * (1. - a);
}

/**
 * @brief Compute \f$ E(z) \f$.
 */
static INLINE double E(double Or, double Om, double Ok, double Ol, double w,
                       double a_inv) {

  return sqrt(Or * a_inv * a_inv * a_inv * a_inv + Om * a_inv * a_inv * a_inv +
              Ok * a_inv * a_inv + Ol * pow(a_inv, 3. * (1. + w)));
}

/**
 * @brief Returns the time (in internal units) since Big Bang at a given
 * scale-factor.
 *
 * @param c The current #cosmology.
 * @param a Scale-factor of interest.
 */
double cosmology_get_time_since_big_bang(const struct cosmology *c, double a) {

  /* Time between a_begin and a */
  const double delta_t =
      interp_table(c->time_interp_table, log(a), c->log_a_begin, c->log_a_end);

  return c->time_interp_table_offset + delta_t;
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
  c->a3_inv = a_inv * a_inv * a_inv;

  /* Redshift */
  c->z = a_inv - 1.;

  /* Dark-energy equation of state */
  c->w = cosmology_dark_energy_EoS(a, c->w_0, c->w_a);

  /* E(z) */
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w = c->w;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w, a_inv);

  /* H(z) */
  c->H = c->H0 * E_z;

  /* Expansion rate */
  c->a_dot = c->H * c->a;

  /* Time */
  c->time = cosmology_get_time_since_big_bang(c, a);
  c->lookback_time = c->universe_age_at_present_day - c->time;
}

/**
 * @brief Computes \f$ dt / a^2 \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double drift_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double w = cosmology_dark_energy_EoS(a, w_0, w_a);
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w, a_inv);
  const double H = H0 * E_z;

  return (1. / H) * a_inv * a_inv * a_inv;
}

/**
 * @brief Computes \f$ dt / a \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double gravity_kick_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double w = cosmology_dark_energy_EoS(a, w_0, w_a);
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w, a_inv);
  const double H = H0 * E_z;

  return (1. / H) * a_inv * a_inv;
}

/**
 * @brief Computes \f$ dt / a^{3(\gamma - 1) + 1} \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double hydro_kick_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double w = cosmology_dark_energy_EoS(a, w_0, w_a);
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w, a_inv);
  const double H = H0 * E_z;

  /* Note: we can't use the pre-defined pow_gamma_xxx() function as
     as we need double precision accuracy for the GSL routine. */
  return (1. / H) * pow(a_inv, 3. * hydro_gamma_minus_one) * a_inv;
}

/**
 * @brief Computes \f$ dt \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double time_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double a_inv = 1. / a;
  const double w = cosmology_dark_energy_EoS(a, w_0, w_a);
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w, a_inv);
  const double H = H0 * E_z;

  return (1. / H) * a_inv;
}

/**
 * @brief Initialise the interpolation tables for the integrals.
 */
void cosmology_init_tables(struct cosmology *c) {

  const ticks tic = getticks();

#ifdef HAVE_LIBGSL

  /* Retrieve some constants */
  const double a_begin = c->a_begin;

  /* Allocate memory for the interpolation tables */
  c->drift_fac_interp_table = malloc(cosmology_table_length * sizeof(double));
  c->grav_kick_fac_interp_table =
      malloc(cosmology_table_length * sizeof(double));
  c->hydro_kick_fac_interp_table =
      malloc(cosmology_table_length * sizeof(double));
  c->time_interp_table = malloc(cosmology_table_length * sizeof(double));

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
  gsl_function F = {&drift_integrand, c};
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

  /* Integrate the kick factor \int_{a_begin}^{a_table[i]} dt/a^(3(g-1)+1) */
  F.function = &hydro_kick_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->hydro_kick_fac_interp_table[i] = result;
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

  /* Integrate the time \int_{0}^{1} dt */
  gsl_integration_qag(&F, 0., 1, 0, 1.0e-13, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->universe_age_at_present_day = result;

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);
  free(a_table);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif

  message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
          clocks_getunit());
}

/**
 * @brief Initialises the #cosmology from the values read in the parameter file.
 *
 * @param params The parsed values.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in the current system of units.
 * @param c The #cosmology to initialise.
 */
void cosmology_init(const struct swift_params *params,
                    const struct unit_system *us,
                    const struct phys_const *phys_const, struct cosmology *c) {

  /* Read in the cosmological parameters */
  c->Omega_m = parser_get_param_double(params, "Cosmology:Omega_m");
  c->Omega_r = parser_get_opt_param_double(params, "Cosmology:Omega_r", 0.);
  c->Omega_lambda = parser_get_param_double(params, "Cosmology:Omega_lambda");
  c->Omega_b = parser_get_param_double(params, "Cosmology:Omega_b");
  c->w_0 = parser_get_opt_param_double(params, "Cosmology:w_0", -1.);
  c->w_a = parser_get_opt_param_double(params, "Cosmology:w_a", 0.);
  c->h = parser_get_param_double(params, "Cosmology:h");

  /* Read the start and end of the simulation */
  c->a_begin = parser_get_param_double(params, "Cosmology:a_begin");
  c->a_end = parser_get_param_double(params, "Cosmology:a_end");
  c->log_a_begin = log(c->a_begin);
  c->log_a_end = log(c->a_end);

  /* Construct derived quantities */

  /* Curvature density (for closure) */
  c->Omega_k = 1. - (c->Omega_m + c->Omega_r + c->Omega_lambda);

  /* Dark-energy equation of state */
  c->w = cosmology_dark_energy_EoS(c->a_begin, c->w_0, c->w_a);

  /* Hubble constant in internal units */
  const double km = 1.e5 / units_cgs_conversion_factor(us, UNIT_CONV_LENGTH);
  const double H0_cgs =
      100. * c->h * (km / (1.e6 * phys_const->const_parsec)); /* s^-1 */
  c->H0 = H0_cgs * units_cgs_conversion_factor(us, UNIT_CONV_TIME);
  c->Hubble_time = 1. / c->H0;

  /* Set remaining variables to invalid values */
  c->H = -1.;
  c->a = -1.;
  c->a_inv = -1;
  c->z = -1.;
  c->a_dot = -1.;
  c->time = -1.;
  c->universe_age_at_present_day = -1.;

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  cosmology_init_tables(c);
}

/**
 * @brief Prints the #cosmology model to stdout.
 */
void cosmology_print(const struct cosmology *c) {

  message(
      "Density parameters: [O_m, O_l, O_b, O_k, O_r] = [%f, %f, %f, %f, %f]",
      c->Omega_m, c->Omega_lambda, c->Omega_b, c->Omega_k, c->Omega_r);
  message("Dark energy equation of state: w_0=%f w_a=%f", c->w_0, c->w_a);
  message("Hubble constant: h = %f, H_0 = %e U_t^(-1)", c->h, c->H0);
  message("Hubble time: 1/H0 = %e U_t", c->Hubble_time);
  message("Universe age at present day: %e U_t", c->universe_age_at_present_day);
}

/**
 * @brief Computes the cosmology factor that enters the drift operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^2 \f$ using the interpolation table.
 *
 * @param e The #engine.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_drift_factor(const struct engine *e,
                                  integertime_t ti_start,
                                  integertime_t ti_end) {

  const struct cosmology *c = e->cosmology;

  const double a_start = c->log_a_begin + ti_start * e->timeBase;
  const double a_end = c->log_a_begin + ti_end * e->timeBase;

  const double int_start = interp_table(c->drift_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->drift_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the gravity kick operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a \f$ using the interpolation table.
 *
 * @param e The #engine.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_grav_kick_factor(const struct engine *e,
                                      integertime_t ti_start,
                                      integertime_t ti_end) {

  const struct cosmology *c = e->cosmology;

  const double a_start = c->log_a_begin + ti_start * e->timeBase;
  const double a_end = c->log_a_begin + ti_end * e->timeBase;

  const double int_start = interp_table(c->grav_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->grav_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the hydro kick operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a \f$ using the interpolation table.
 *
 * @param e The #engine.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_hydro_kick_factor(const struct engine *e,
                                       integertime_t ti_start,
                                       integertime_t ti_end) {

  const struct cosmology *c = e->cosmology;

  const double a_start = c->log_a_begin + ti_start * e->timeBase;
  const double a_end = c->log_a_begin + ti_end * e->timeBase;

  const double int_start = interp_table(c->hydro_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->hydro_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}
