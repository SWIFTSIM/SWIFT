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
#include "common_io.h"
#include "inline.h"
#include "restart.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#endif

/*! Number of values stored in the cosmological interpolation tables */
const int cosmology_table_length = 10000;

#ifdef HAVE_LIBGSL
/*! Size of the GSL workspace */
const size_t GSL_workspace_size = 100000;
#endif

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
 * @brief Computes the integral of the dark-energy equation of state
 * up to a scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 * and compute \f$ \tilde{w}(a) = \int_0^a\frac{1 + w(z)}{1+z}dz \f$.
 *
 * @param a The current scale-factor.
 * @param w0 The equation of state parameter at z=0
 * @param wa The equation of state evolution parameter
 */
static INLINE double w_tilde(double a, double w0, double wa) {
  return (a - 1.) * wa - (1. + w0 + wa) * log(a);
}

/**
 * @brief Compute \f$ E(z) \f$.
 */
static INLINE double E(double Or, double Om, double Ok, double Ol, double w0,
                       double wa, double a) {
  const double a_inv = 1. / a;
  return sqrt(Or * a_inv * a_inv * a_inv * a_inv + Om * a_inv * a_inv * a_inv +
              Ok * a_inv * a_inv + Ol * exp(3. * w_tilde(a, w0, wa)));
}

/**
 * @brief Returns the time (in internal units) since Big Bang at a given
 * scale-factor.
 *
 * @param c The current #cosmology.
 * @param a Scale-factor of interest.
 */
double cosmology_get_time_since_big_bang(const struct cosmology *c, double a) {

#ifdef SWIFT_DEBUG_CHECKS
  if (a < c->a_begin) error("Error a can't be smaller than a_begin");
#endif

  /* Time between a_begin and a */
  const double delta_t =
      interp_table(c->time_interp_table, log(a), c->log_a_begin, c->log_a_end);

  return c->time_interp_table_offset + delta_t;
}

/**
 * @brief Update the cosmological parameters to the current simulation time.
 *
 * @param c The #cosmology struct.
 * @param phys_const The physical constants in the internal units.
 * @param ti_current The current (integer) time.
 */
void cosmology_update(struct cosmology *c, const struct phys_const *phys_const,
                      integertime_t ti_current) {

  /* Save the previous state */
  c->z_old = c->z;
  c->a_old = c->a;

  /* Get scale factor and powers of it */
  const double a = c->a_begin * exp(ti_current * c->time_base);
  const double a_inv = 1. / a;
  c->a = a;
  c->a_inv = a_inv;
  c->a2_inv = a_inv * a_inv;
  c->a3_inv = a_inv * a_inv * a_inv;
  c->a_factor_internal_energy =
      pow(a, -3. * hydro_gamma_minus_one);          /* a^{3*(1-gamma)} */
  c->a_factor_pressure = pow(a, -3. * hydro_gamma); /* a^{-3*gamma} */
  c->a_factor_sound_speed =
      pow(a, -1.5 * hydro_gamma_minus_one); /* a^{3*(1-gamma)/2} */
  c->a_factor_grav_accel = a_inv * a_inv;   /* 1 / a^2 */
  c->a_factor_hydro_accel =
      pow(a, -3. * hydro_gamma + 2.); /* 1 / a^(3*gamma - 2) */
  c->a_factor_mu =
      pow(a, 0.5 * (3. * hydro_gamma - 5.)); /* a^{(3*gamma - 5) / 2} */
  c->a_factor_Balsara_eps =
      pow(a, 0.5 * (1. - 3. * hydro_gamma)); /* a^{(1 - 3*gamma) / 2} */

  /* Redshift */
  c->z = a_inv - 1.;

  /* Dark-energy equation of state */
  c->w = cosmology_dark_energy_EoS(a, c->w_0, c->w_a);

  /* E(z) */
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w0 = c->w_0;
  const double wa = c->w_a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w0, wa, a);

  /* H(z) */
  c->H = c->H0 * E_z;

  /* Expansion rate */
  c->a_dot = c->H * c->a;

  /* Critical density */
  c->critical_density =
      3. * c->H * c->H / (8. * M_PI * phys_const->const_newton_G);

  /* Time-step conversion factor */
  c->time_step_factor = c->H;

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
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
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
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
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
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  /* Note: we can't use the pre-defined pow_gamma_xxx() function as
     as we need double precision accuracy for the GSL routine. */
  return (1. / H) * pow(a_inv, 3. * hydro_gamma_minus_one) * a_inv;
}

/**
 * @brief Computes \f$a dt\f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double hydro_kick_corr_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_r = c->Omega_r;
  const double Omega_m = c->Omega_m;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;

  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return 1. / H;
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
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (1. / H) * a_inv;
}

/**
 * @brief Initialise the interpolation tables for the integrals.
 */
void cosmology_init_tables(struct cosmology *c) {

#ifdef HAVE_LIBGSL

  /* Retrieve some constants */
  const double a_begin = c->a_begin;

  /* Allocate memory for the interpolation tables */
  c->drift_fac_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));
  c->grav_kick_fac_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));
  c->hydro_kick_fac_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));
  c->hydro_kick_corr_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));
  c->time_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));
  c->scale_factor_interp_table =
      (double *)malloc(cosmology_table_length * sizeof(double));

  /* Prepare a table of scale factors for the integral bounds */
  const double delta_a =
      (c->log_a_end - c->log_a_begin) / cosmology_table_length;
  double *a_table = (double *)malloc(cosmology_table_length * sizeof(double));
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

  /* Integrate the kick correction factor \int_{a_begin}^{a_table[i]} a dt */
  F.function = &hydro_kick_corr_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->hydro_kick_corr_interp_table[i] = result;
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

  /* Update the times */
  c->time_begin = cosmology_get_time_since_big_bang(c, c->a_begin);
  c->time_end = cosmology_get_time_since_big_bang(c, c->a_end);

  /*
   * Inverse t(a)
   */

  const double delta_t = (c->time_end - c->time_begin) / cosmology_table_length;

  /* index in the time_interp_table */
  int i_a = 0;

  for (int i_time = 0; i_time < cosmology_table_length; i_time++) {
    /* Current time
     * time_interp_table = \int_a_begin^a => no need of time_begin */
    double time_interp = delta_t * (i_time + 1);

    /* Find next time in time_interp_table */
    while (i_a < cosmology_table_length &&
           c->time_interp_table[i_a] <= time_interp) {
      i_a++;
    }

    /* Find linear interpolation scaling */
    double scale = 0;
    if (i_a != cosmology_table_length) {
      scale = time_interp - c->time_interp_table[i_a - 1];
      scale /= c->time_interp_table[i_a] - c->time_interp_table[i_a - 1];
    }

    scale += i_a;

    /* Compute interpolated scale factor */
    double log_a = c->log_a_begin + scale * (c->log_a_end - c->log_a_begin) /
                                        cosmology_table_length;
    c->scale_factor_interp_table[i_time] = exp(log_a) - c->a_begin;
  }

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);
  free(a_table);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif
}

/**
 * @brief Initialises the #cosmology from the values read in the parameter file.
 *
 * @param params The parsed values.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in the current system of units.
 * @param c The #cosmology to initialise.
 */
void cosmology_init(struct swift_params *params, const struct unit_system *us,
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
  c->time_base = (c->log_a_end - c->log_a_begin) / max_nr_timesteps;
  c->time_base_inv = 1. / c->time_base;

  /* If a_begin == a_end we hang */

  if (c->a_begin >= c->a_end)
    error("a_begin must be strictly before (and not equal to) a_end");

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

  /* Critical density at present day */
  c->critical_density_0 =
      3. * c->H0 * c->H0 / (8. * M_PI * phys_const->const_newton_G);

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  cosmology_init_tables(c);

  /* Set remaining variables to alid values */
  cosmology_update(c, phys_const, 0);

  /* Update the times */
  c->time_begin = cosmology_get_time_since_big_bang(c, c->a_begin);
  c->time_end = cosmology_get_time_since_big_bang(c, c->a_end);

  /* Initialise the old values to a valid state */
  c->a_old = c->a_begin;
  c->z_old = 1. / c->a_old - 1.;
}

/**
 * @brief Initialise the #cosmology for non-cosmological time-integration
 *
 * Essentially sets all constants to 1 or 0.
 *
 * @param c The #cosmology to initialise.
 */
void cosmology_init_no_cosmo(struct cosmology *c) {

  c->Omega_m = 0.;
  c->Omega_r = 0.;
  c->Omega_k = 0.;
  c->Omega_lambda = 0.;
  c->Omega_b = 0.;
  c->w_0 = 0.;
  c->w_a = 0.;
  c->h = 1.;
  c->w = -1.;

  c->a_begin = 1.;
  c->a_end = 1.;
  c->log_a_begin = 0.;
  c->log_a_end = 0.;

  c->H = 0.;
  c->H0 = 0.;
  c->a = 1.;
  c->z = 0.;
  c->a_inv = 1.;
  c->a2_inv = 1.;
  c->a3_inv = 1.;
  c->a_factor_internal_energy = 1.;
  c->a_factor_pressure = 1.;
  c->a_factor_sound_speed = 1.;
  c->a_factor_mu = 1.;
  c->a_factor_Balsara_eps = 1.;
  c->a_factor_hydro_accel = 1.;
  c->a_factor_grav_accel = 1.;

  c->a_old = 1.;
  c->z_old = 0.;

  c->critical_density = 0.;
  c->critical_density_0 = 0.;

  c->time_step_factor = 1.;

  c->a_dot = 0.;
  c->time = 0.;
  c->universe_age_at_present_day = 0.;
  c->Hubble_time = 0.;
  c->lookback_time = 0.;

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->hydro_kick_corr_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  c->scale_factor_interp_table = NULL;

  c->time_begin = 0.;
  c->time_end = 0.;
}

/**
 * @brief Computes the cosmology factor that enters the drift operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^2 \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_drift_factor(const struct cosmology *c,
                                  integertime_t ti_start,
                                  integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

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
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_grav_kick_factor(const struct cosmology *c,
                                      integertime_t ti_start,
                                      integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->grav_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->grav_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the hydro kick operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^{3(gamma - 1)} \f$ using the
 * interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_hydro_kick_factor(const struct cosmology *c,
                                       integertime_t ti_start,
                                       integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->hydro_kick_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->hydro_kick_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the hydro kick correction
 * operator for the meshless schemes (GIZMO-MFV).
 *
 * Computes \f$ \int_{a_start}^{a_end} a dt \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_corr_kick_factor(const struct cosmology *c,
                                      integertime_t ti_start,
                                      integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->hydro_kick_corr_interp_table,
                                        a_start, c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->hydro_kick_corr_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Computes the cosmology factor that enters the thermal variable kick
 * operator.
 *
 * Computes \f$ \int_{a_start}^{a_end} dt/a^2 \f$ using the interpolation table.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start of the drift.
 * @param ti_end the (integer) time of the end of the drift.
 */
double cosmology_get_therm_kick_factor(const struct cosmology *c,
                                       integertime_t ti_start,
                                       integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double a_start = c->log_a_begin + ti_start * c->time_base;
  const double a_end = c->log_a_begin + ti_end * c->time_base;

  const double int_start = interp_table(c->drift_fac_interp_table, a_start,
                                        c->log_a_begin, c->log_a_end);
  const double int_end = interp_table(c->drift_fac_interp_table, a_end,
                                      c->log_a_begin, c->log_a_end);

  return int_end - int_start;
}

/**
 * @brief Compute the cosmic time (in internal units) between two points
 * on the integer time line.
 *
 * @param c The current #cosmology.
 * @param ti_start the (integer) time of the start.
 * @param ti_end the (integer) time of the end.
 */
double cosmology_get_delta_time(const struct cosmology *c,
                                integertime_t ti_start, integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_start) error("ti_end must be >= ti_start");
#endif

  const double log_a_start = c->log_a_begin + ti_start * c->time_base;
  const double log_a_end = c->log_a_begin + ti_end * c->time_base;

  /* Time between a_begin and a_start */
  const double t1 = interp_table(c->time_interp_table, log_a_start,
                                 c->log_a_begin, c->log_a_end);

  /* Time between a_begin and a_end */
  const double t2 = interp_table(c->time_interp_table, log_a_end,
                                 c->log_a_begin, c->log_a_end);

  return t2 - t1;
}

/**
 * @brief Compute scale factor from time since big bang (in internal units).
 *
 * @param c The current #cosmology.
 * @param t time since the big bang
 * @return The scale factor.
 */
double cosmology_get_scale_factor(const struct cosmology *c, double t) {
  /* scale factor between time_begin and t */
  const double a =
      interp_table(c->scale_factor_interp_table, t, c->time_interp_table_offset,
                   c->universe_age_at_present_day);
  return a + c->a_begin;
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
  message("Universe age at present day: %e U_t",
          c->universe_age_at_present_day);
}

void cosmology_clean(struct cosmology *c) {

  free(c->drift_fac_interp_table);
  free(c->grav_kick_fac_interp_table);
  free(c->hydro_kick_fac_interp_table);
  free(c->hydro_kick_corr_interp_table);
  free(c->time_interp_table);
  free(c->scale_factor_interp_table);
}

#ifdef HAVE_HDF5
void cosmology_write_model(hid_t h_grp, const struct cosmology *c) {

  io_write_attribute_d(h_grp, "a_beg", c->a_begin);
  io_write_attribute_d(h_grp, "a_end", c->a_end);
  io_write_attribute_d(h_grp, "time_beg [internal units]", c->time_begin);
  io_write_attribute_d(h_grp, "time_end [internal units]", c->time_end);
  io_write_attribute_d(h_grp, "Universe age [internal units]", c->time);
  io_write_attribute_d(h_grp, "Lookback time [internal units]",
                       c->lookback_time);
  io_write_attribute_d(h_grp, "h", c->h);
  io_write_attribute_d(h_grp, "H0 [internal units]", c->H0);
  io_write_attribute_d(h_grp, "H [internal units]", c->H);
  io_write_attribute_d(h_grp, "Hubble time [internal units]", c->Hubble_time);
  io_write_attribute_d(h_grp, "Omega_m", c->Omega_m);
  io_write_attribute_d(h_grp, "Omega_r", c->Omega_r);
  io_write_attribute_d(h_grp, "Omega_b", c->Omega_b);
  io_write_attribute_d(h_grp, "Omega_k", c->Omega_k);
  io_write_attribute_d(h_grp, "Omega_lambda", c->Omega_lambda);
  io_write_attribute_d(h_grp, "w_0", c->w_0);
  io_write_attribute_d(h_grp, "w_a", c->w_a);
  io_write_attribute_d(h_grp, "w", c->w);
  io_write_attribute_d(h_grp, "Redshift", c->z);
  io_write_attribute_d(h_grp, "Scale-factor", c->a);
  io_write_attribute_d(h_grp, "Critical density [internal units]",
                       c->critical_density);
}
#endif

/**
 * @brief Write a cosmology struct to the given FILE as a stream of bytes.
 *
 * @param cosmology the struct
 * @param stream the file stream
 */
void cosmology_struct_dump(const struct cosmology *cosmology, FILE *stream) {
  restart_write_blocks((void *)cosmology, sizeof(struct cosmology), 1, stream,
                       "cosmology", "cosmology function");
}

/**
 * @brief Restore a cosmology struct from the given FILE as a stream of
 * bytes.
 *
 * @param enabled whether cosmology is enabled.
 * @param cosmology the struct
 * @param stream the file stream
 */
void cosmology_struct_restore(int enabled, struct cosmology *cosmology,
                              FILE *stream) {
  restart_read_blocks((void *)cosmology, sizeof(struct cosmology), 1, stream,
                      NULL, "cosmology function");

  /* Re-initialise the tables if using a cosmology. */
  if (enabled) cosmology_init_tables(cosmology);
}
