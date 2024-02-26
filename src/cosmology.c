/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "align.h"
#include "common_io.h"
#include "inline.h"
#include "memuse.h"
#include "minmax.h"
#include "restart.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#endif

/*! Number of values stored in the cosmological interpolation tables */
const int cosmology_table_length = 30000;

#ifdef HAVE_LIBGSL
/*! Size of the GSL workspace */
const size_t GSL_workspace_size = 100000;
#endif

/**
 * @brief Returns the interpolated value from a table.
 *
 * Uses linear interpolation.
 *
 * @param table The table of value to interpolate from (should be of length
 * cosmology_table_length).
 * @param x The value to interpolate at.
 * @param x_min The mininum of the range of x.
 * @param x_max The maximum of the range of x.
 */
static INLINE double interp_table(const double *table, const double x,
                                  const double x_min, const double x_max) {

  const double xx =
      ((x - x_min) / (x_max - x_min)) * ((double)cosmology_table_length);

  const int i = (int)xx;
  const int ii = min(cosmology_table_length - 1, i);

  /* Indicate that the whole array is aligned on boundaries */
  swift_align_information(double, table, SWIFT_STRUCT_ALIGNMENT);

  if (ii < 1)
    return table[0] * xx;
  else
    return table[ii - 1] + (table[ii] - table[ii - 1]) * (xx - ii);
}

/**
 * @brief Invert a function y(a) which is tabulated at intervals in log(a)
 *
 * The function to invert must be monotonically increasing and is
 * assumed to be zero at a=a_begin.
 *
 * @param y_table Input array with cosmology_table_length elements. Element i
 * contains the value of y at log(a)=log(a_begin)+delta_log_a*(i+1).
 * @param log_a_begin Log of expansion factor at the start of the interval
 * @param delta_y Interval in y at which to tabulate a-a_begin in a_table
 * @param delta_log_a Interval in log(a) at which the function is tabulated in
 *        y_table
 * @param a_table Output array with cosmology_table_length elements. Element i
 *        contains the value of a-a_begin at which y=delta_y*(i+1).
 *
 */
#ifdef HAVE_LIBGSL
static void invert_table(const double *y_table, const double log_a_begin,
                         const double delta_y, const double delta_log_a,
                         double *a_table) {

  int i_a = 0;
  for (int i_y = 0; i_y < cosmology_table_length; i_y++) {

    double y_interp = delta_y * (i_y + 1);

    /* Find next y in tabulated y(a) */
    while (i_a < cosmology_table_length && y_table[i_a] <= y_interp) {
      i_a++;
    }

    /* Find y values we're interpolating between */
    double scale = 0.0;
    if (i_a == 0) {
      /* We're interpolating between y=0 and the first tabulated point  */
      double y1 = 0.0;
      double y2 = y_table[i_a];
      scale = (y_interp - y1) / (y2 - y1) + i_a;
    } else if ((i_a > 0) && (i_a < cosmology_table_length)) {
      /* We're interpolating between two tabulated points in the array */
      double y1 = y_table[i_a - 1];
      double y2 = y_table[i_a];
      scale = (y_interp - y1) / (y2 - y1) + i_a;
    } else if (i_a == cosmology_table_length) {
      /* This happens when y_interp equals the final tabulated point */
      scale = i_a;
    } else {
      error("Interpolating function to invert outside tabulated range!");
    }

    /* Compute log(a) at this point */
    const double log_a = log_a_begin + scale * delta_log_a;

    /* Store value of a-a_begin corresponding to y=y_interp */
    a_table[i_y] = exp(log_a) - exp(log_a_begin);
  }
}
#endif /* HAVE_LIBGSL */

/**
 * @brief Computes the dark-energy equation of state at a given scale-factor a.
 *
 * We follow the convention of Linder & Jenkins, MNRAS, 346, 573, 2003
 *
 * @param a The current scale-factor
 * @param w_0 The equation of state parameter at z=0
 * @param w_a The equation of state evolution parameter
 */
__attribute__((const)) static INLINE double cosmology_dark_energy_EoS(
    const double a, const double w_0, const double w_a) {

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
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 */
__attribute__((const)) static INLINE double w_tilde(const double a,
                                                    const double w0,
                                                    const double wa) {
  return (a - 1.) * wa - (1. + w0 + wa) * log(a);
}

/**
 * @brief Compute \f$ E(z) \f$.
 *
 * @param Omega_r The radiation density parameter \f$ \Omega_r \f$.
 * @param Omega_m The matter density parameter \f$ \Omega_m \f$.
 * @param Omega_k The curvature density parameter \f$ \Omega_k \f$.
 * @param Omega_l The cosmological constant density parameter \f$ \Omega_\Lambda
 * \f$.
 * @param w0 The equation of state parameter at z=0.
 * @param wa The equation of state evolution parameter.
 * @param a The current scale-factor.
 */
__attribute__((const)) static INLINE double E(
    const double Omega_r, const double Omega_m, const double Omega_k,
    const double Omega_l, const double w0, const double wa, const double a) {

  const double a_inv = 1. / a;

  return sqrt(Omega_r * a_inv * a_inv * a_inv * a_inv + /* Radiation */
              Omega_m * a_inv * a_inv * a_inv +         /* Matter */
              Omega_k * a_inv * a_inv +                 /* Curvature */
              Omega_l * exp(3. * w_tilde(a, w0, wa)));  /* Lambda */
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

  /* Update the neutrino density */
  c->Omega_nu = cosmology_get_neutrino_density(c, a);

  /* E(z) */
  const double Omega_r = c->Omega_r + c->Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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

  /* Mean density */
  c->mean_density = c->critical_density_0 * c->a3_inv;

  /* Mean matter density */
  c->mean_density_Omega_m = c->mean_density * Omega_m;

  /* Mean baryonic density */
  c->mean_density_Omega_b = c->mean_density * c->Omega_b;

  /* Over-density threshold for virialization
   * Fitting function from Bryan & Norman, 1998, ApJ, 495, 1, 80-99
   * Equation 6. */
  const double x = Omega_m * c->a3_inv / (E_z * E_z) - 1.;
  c->overdensity_BN98 = 18. * M_PI * M_PI + 82. * x - 39 * x * x;

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
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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
 * @brief Computes \f$ c dt / a \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double comoving_distance_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
  const double Omega_k = c->Omega_k;
  const double Omega_l = c->Omega_lambda;
  const double w_0 = c->w_0;
  const double w_a = c->w_a;
  const double H0 = c->H0;
  const double const_speed_light_c = c->const_speed_light_c;
  const double a_inv = 1. / a;
  const double E_z = E(Omega_r, Omega_m, Omega_k, Omega_l, w_0, w_a, a);
  const double H = H0 * E_z;

  return (const_speed_light_c / H) * a_inv * a_inv;
}

/**
 * @brief Computes \f$ dt / a^{3(\gamma - 1) + 1} \f$ for the current cosmology.
 *
 * @param a The scale-factor of interest.
 * @param param The current #cosmology.
 */
double hydro_kick_integrand(double a, void *param) {

  const struct cosmology *c = (const struct cosmology *)param;
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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
  const double Omega_nu = cosmology_get_neutrino_density(c, a);
  const double Omega_r = c->Omega_r + Omega_nu;
  const double Omega_m = c->Omega_cdm + c->Omega_b;
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
 * @brief Evaluates the neutrino density momentum integrand
 * \f$ x^2 \sqrt{x^2 + y^2} / (1+e^x) \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift.
 *
 * This is used to evaluate the integral on (0, 1).
 *
 * @param x The momentum integration variable
 * @param param Neutrino mass y scaled by temperature at redshift of interest.
 * @return The integrand evaluated at x
 */
double neutrino_density_integrand(double x, void *param) {
  double y = *(double *)param;
  double numerator = x * x * hypot(x, y);

  /* Handle overflows */
  if (x > 20 + log(numerator)) {
    return numerator * exp(-x);
  }

  return numerator / (1.0 + exp(x));
}

/**
 * @brief Evaluates the transformed neutrino density momentum integrand
 * \f$ w^{-4} \sqrt{w^{-2} + y^2} / (1+e^{-w}) \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift.
 *
 * This is used to evaluate the integral on (1, infinity).
 *
 * @param w The transformed momentum integration variable w=1/x
 * @param param Neutrino mass y scaled by temperature at redshift of interest.
 * @return The integrand evaluated at w
 */
double neutrino_density_integrand_transformed(double w, void *param) {
  return neutrino_density_integrand(1. / w, param) / (w * w);
}

#ifdef HAVE_LIBGSL

/**
 * @brief Performs the neutrino density momentum integral
 * \f$ \int_0^\infty x^2 \sqrt{x^2 + y^2} / (1+e^x) dx \f$, where
 * \f$ y = a*M_nu/(kb*T) \f$ is the neutrino mass scaled with redshift,
 * without pre-factors.
 *
 * @param space The GSL working space
 * @param y Neutrino mass y scaled by temperature at redshift of interest.
 * @return The integral evaluated at y
 */
double neutrino_density_integrate(gsl_integration_workspace *space, double y) {
  double intermediate, abserr;

  double result = 0;

  gsl_function F1 = {&neutrino_density_integrand, &y};
  gsl_function F2 = {&neutrino_density_integrand_transformed, &y};
  /* Integrate between 0 and 1 */
  gsl_integration_qag(&F1, 0.0, 1.0, 0, 1.0e-13, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &intermediate, &abserr);

  result += intermediate;
  /* Integrate between 1 and infinity */
  gsl_integration_qag(&F2, 0.0, 1.0, 0, 1.0e-13, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &intermediate, &abserr);
  result += intermediate;

  return result;
}

#endif

/**
 * @brief Find a time when all neutrinos are still relativistic. Store the
 * starting, mid, and end points of the neutrino interpolation tables in c.
 *
 * @param c The cosmology structure
 * @param tol Tolerance in density integral
 */
void neutrino_find_relativistic_redshift(struct cosmology *c, double tol) {

  /* Find the largest neutrino mass */
  double M_max_eV = c->M_nu_eV[0];
  for (int i = 1; i < c->N_nu; i++) {
    M_max_eV = fmax(M_max_eV, c->M_nu_eV[i]);
  }

  /* A safe starting time when neutrinos are relativistic */
  double a_safe = 0.5 * tol / M_max_eV;

  /* Dont start the early table later than the simulation */
  a_safe = fmin(a_safe, 0.9 * c->a_begin);

  /* Start the late table just before the start of the simulation */
  double a_midpoint = 0.99 * c->a_begin;

  /* End the late table today (a=1) or at a_end, whichever is later */
  double a_final = fmax(1.0, c->a_end);

  /* Integrate early table on (a_start, a_mid) and late table on (a_mid, a_f) */
  c->log_a_long_begin = log(a_safe);
  c->log_a_long_mid = log(a_midpoint);
  c->log_a_long_end = log(a_final);
}

/**
 * @brief Initialise the neutrino density interpolation tables (early and late).
 */
void cosmology_init_neutrino_tables(struct cosmology *c) {

  /* Skip if we have no massive neutrinos */
  if (c->N_nu == 0) return;

#ifdef HAVE_LIBGSL

  /* Allocate memory for longer interpolation tables */
  if (swift_memalign("cosmo.table", (void **)&c->neutrino_density_early_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->neutrino_density_late_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");

  /* Initalise the GSL workspace */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(GSL_workspace_size);

  /* Find a safe redshift to start the neutrino density interpolation table */
  neutrino_find_relativistic_redshift(c, 1e-7);

  const double pre_factor = 15. * pow(c->T_nu_0 * M_1_PI / c->T_CMB_0, 4);
  const double early_delta_a =
      (c->log_a_long_mid - c->log_a_long_begin) / cosmology_table_length;
  const double late_delta_a =
      (c->log_a_long_end - c->log_a_long_mid) / cosmology_table_length;

  double result;

  /* Fill the early neutrino density table between (a_long_begin, a_long_mid) */
  for (int i = 0; i < cosmology_table_length; i++) {
    double O_nu = 0.;
    double a = exp(c->log_a_long_begin + early_delta_a * (i + 1));

    /* Integrate the FD distribtution for each species */
    for (int j = 0; j < c->N_nu; j++) {
      double y = a * c->M_nu_eV[j] / c->T_nu_0_eV;
      result = neutrino_density_integrate(space, y);
      O_nu += c->deg_nu[j] * result * pre_factor * c->Omega_g;
    }

    c->neutrino_density_early_table[i] = O_nu;
  }

  /* Fill the late neutrino density table between (a_long_mid, a_long_end) */
  for (int i = 0; i < cosmology_table_length; i++) {
    double O_nu = 0.;
    double a = exp(c->log_a_long_mid + late_delta_a * (i + 1));

    /* Integrate the FD distribtution for each species */
    for (int j = 0; j < c->N_nu; j++) {
      double y = a * c->M_nu_eV[j] / c->T_nu_0_eV;
      result = neutrino_density_integrate(space, y);
      O_nu += c->deg_nu[j] * result * pre_factor * c->Omega_g;
    }

    c->neutrino_density_late_table[i] = O_nu;
  }

  /* Free the workspace */
  gsl_integration_workspace_free(space);

#else

  error("Code not compiled with GSL. Can't compute cosmology integrals.");

#endif
}

/**
 * @brief Initialise the interpolation tables for the integrals.
 */
void cosmology_init_tables(struct cosmology *c) {

#ifdef HAVE_LIBGSL

  /* Retrieve some constants */
  const double a_begin = c->a_begin;
  const double a_end = c->a_end;

  /* Allocate memory for the interpolation tables */
  if (swift_memalign("cosmo.table", (void **)&c->drift_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->grav_kick_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->hydro_kick_fac_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->hydro_kick_corr_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->time_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->scale_factor_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign("cosmo.table", (void **)&c->comoving_distance_interp_table,
                     SWIFT_STRUCT_ALIGNMENT,
                     cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");
  if (swift_memalign(
          "cosmo.table", (void **)&c->comoving_distance_inverse_interp_table,
          SWIFT_STRUCT_ALIGNMENT, cosmology_table_length * sizeof(double)) != 0)
    error("Failed to allocate cosmology interpolation table");

  /* Prepare a table of scale factors for the integral bounds */
  const double delta_a =
      (c->log_a_end - c->log_a_begin) / cosmology_table_length;
  double *a_table = (double *)swift_malloc(
      "cosmo.table", cosmology_table_length * sizeof(double));
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

  /* Integrate the time \int_{0}^{a_end} dt */
  gsl_integration_qag(&F, 0., a_end, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->time_interp_table_max = result;

  /* Integrate the time \int_{0}^{1} dt */
  gsl_integration_qag(&F, 0., 1, 0, 1.0e-13, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->universe_age_at_present_day = result;

  /* Integrate the comoving distance \int_{a_begin}^{a_table[i]} c dt/a */
  F.function = &comoving_distance_integrand;
  for (int i = 0; i < cosmology_table_length; i++) {
    gsl_integration_qag(&F, a_begin, a_table[i], 0, 1.0e-10, GSL_workspace_size,
                        GSL_INTEG_GAUSS61, space, &result, &abserr);

    /* Store result */
    c->comoving_distance_interp_table[i] = result;
  }

  /* Integrate the comoving distance \int_{a_begin}^{1.0} c dt/a */
  F.function = &comoving_distance_integrand;
  gsl_integration_qag(&F, a_begin, 1.0, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->comoving_distance_interp_table_offset = result;

  /* Integrate the comoving distance \int_{a_begin}^{a_end} c dt/a */
  F.function = &comoving_distance_integrand;
  gsl_integration_qag(&F, a_begin, a_end, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);
  c->comoving_distance_start_to_end = result;

  /* Update the times */
  c->time_begin = cosmology_get_time_since_big_bang(c, c->a_begin);
  c->time_end = cosmology_get_time_since_big_bang(c, c->a_end);

  /* Interval in log(a) at which the time and comoving distance functions are
   * tabulated */
  const double delta_log_a =
      (c->log_a_end - c->log_a_begin) / cosmology_table_length;

  /* Tabulate inverted t(a) function */
  const double delta_t = (c->time_end - c->time_begin) / cosmology_table_length;
  invert_table(c->time_interp_table, c->log_a_begin, delta_t, delta_log_a,
               c->scale_factor_interp_table);

  /* Tabulate inverted comoving distance function */
  const double r_begin = cosmology_get_comoving_distance(c, a_begin);
  const double r_end = cosmology_get_comoving_distance(c, a_end);
  const double delta_r = (r_begin - r_end) / cosmology_table_length;
  invert_table(c->comoving_distance_interp_table, c->log_a_begin, delta_r,
               delta_log_a, c->comoving_distance_inverse_interp_table);

  /* Free the workspace and temp array */
  gsl_integration_workspace_free(space);
  swift_free("cosmo.table", a_table);

#ifdef SWIFT_DEBUG_CHECKS

  const int n = 1000 * cosmology_table_length;
  double max_error_time = 0;
  double max_error_distance = 0;

  for (int i = 0; i < n; i += 1) {

    double a_check, frac_error;

    /* Choose value of expansion factor for check */
    const double dloga = (c->log_a_end - c->log_a_begin) / (n - 1);
    double a = exp(c->log_a_begin + dloga * i);
    a = fmax(a, c->a_begin);
    a = fmin(a, c->a_end);

    /* Verify that converting expansion factor to time and back recovers the
     * original value */
    const double t = cosmology_get_time_since_big_bang(c, a);
    a_check = cosmology_get_scale_factor(c, t);
    frac_error = fabs(a_check / a - 1.0);
    if (frac_error > max_error_time) max_error_time = frac_error;

    /* Verify that converting expansion factor to comoving distance and back
     * recovers the original value */
    const double r = cosmology_get_comoving_distance(c, a);
    a_check = cosmology_scale_factor_at_comoving_distance(c, r);
    frac_error = fabs(a_check / a - 1.0);
    if (frac_error > max_error_distance) max_error_distance = frac_error;
  }

  message("Max fractional error in a to age of universe round trip = %16.8e\n",
          max_error_time);
  message(
      "Max fractional error in a to comoving distance round trip = %16.8e\n",
      max_error_distance);

#endif /* SWIFT_DEBUG_CHECKS */

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

  /* Check first for outdated parameter files still giving Omega_m */
  const double test_Omega_m =
      parser_get_opt_param_double(params, "Cosmology:Omega_m", -1.);
  if (test_Omega_m != -1.)
    error(
        "Parameter file contains Cosmology:Omega_m. This is deprecated. Please "
        "specify Omega_cdm (the cold dark matter density parameter) and "
        "optionally neutrino parameters.\nIf that simulation did not use "
        "neutrinos then the new Omega_cdm parameter should just be (old) "
        "Omega_m - Omega_b.");

  /* Read in the cosmological parameters */
  c->Omega_cdm = parser_get_param_double(params, "Cosmology:Omega_cdm");
  c->Omega_r = parser_get_opt_param_double(params, "Cosmology:Omega_r", 0.);
  c->Omega_lambda = parser_get_param_double(params, "Cosmology:Omega_lambda");
  c->Omega_b = parser_get_param_double(params, "Cosmology:Omega_b");
  c->w_0 = parser_get_opt_param_double(params, "Cosmology:w_0", -1.);
  c->w_a = parser_get_opt_param_double(params, "Cosmology:w_a", 0.);
  c->h = parser_get_param_double(params, "Cosmology:h");

  /* Neutrino temperature (inferred from T_CMB_0 if not specified) */
  c->T_nu_0 = parser_get_opt_param_double(params, "Cosmology:T_nu_0", 0.);

  /* Number of ultra-relativistic (massless) and massive neutrino species */
  c->N_ur = parser_get_opt_param_double(params, "Cosmology:N_ur", 0.);
  c->N_nu = parser_get_opt_param_int(params, "Cosmology:N_nu", 0);

  /* Make sure that the cosmological parameters are not overdetermined */
  if (c->Omega_r != 0. && c->N_ur != 0.) {
    error("Cannot use both Cosmology:Omega_r and Cosmology:N_ur.");
  }

  /* If there are massive neutrinos, read the masses and degeneracies */
  if (c->N_nu > 0) {
    c->M_nu_eV = (double *)swift_malloc("Mnu", c->N_nu * sizeof(double));
    c->deg_nu = (double *)swift_malloc("degnu", c->N_nu * sizeof(double));

    /* Set default values */
    for (int i = 0; i < c->N_nu; i++) {
      c->M_nu_eV[i] = 0.0;
      c->deg_nu[i] = 1.0;
    }

    parser_get_opt_param_double_array(params, "Cosmology:M_nu_eV", c->N_nu,
                                      c->M_nu_eV);
    parser_get_opt_param_double_array(params, "Cosmology:deg_nu", c->N_nu,
                                      c->deg_nu);
  }

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

  /* Some constants */
  const double cc = phys_const->const_speed_light_c;
  const double rho_c3_on_4sigma = c->critical_density_0 * cc * cc * cc /
                                  (4. * phys_const->const_stefan_boltzmann);

  /* Store speed of light in internal units */
  c->const_speed_light_c = phys_const->const_speed_light_c;

  /* Handle neutrinos only if present */
  if (c->N_ur == 0. && c->N_nu == 0) {
    /* Infer T_CMB_0 from Omega_r */
    c->T_CMB_0 = pow(c->Omega_r * rho_c3_on_4sigma, 1. / 4.);
    c->T_CMB_0_K =
        c->T_CMB_0 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

    c->Omega_g = c->Omega_r;
    c->T_nu_0 = 0.;
    c->T_nu_0_eV = 0.;
    c->Omega_ur = 0.;
    c->Omega_nu_0 = 0.;
    c->Omega_nu = 0.;
    c->N_eff = 0.;
    c->deg_nu_tot = 0.;

    c->neutrino_density_early_table = NULL;
    c->neutrino_density_late_table = NULL;

  } else {

    /* Infer T_CMB_0 from Omega_r if the latter is specified */
    if (c->Omega_r != 0) {
      c->T_CMB_0 = pow(c->Omega_r * rho_c3_on_4sigma, 1. / 4.);
    } else {
      c->T_CMB_0 = phys_const->const_T_CMB_0;
    }

    /* Approximate the neutrino temperature if unspecified */
    const double decoupling_factor = cbrt(4. / 11);
    const double dec_4 = pow(decoupling_factor, 4);
    if (c->T_nu_0 == 0.) {
      c->T_nu_0 = c->T_CMB_0 * decoupling_factor;
    }

    /* Unit conversions */
    c->T_CMB_0_K =
        c->T_CMB_0 / units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
    c->T_nu_0_eV = c->T_nu_0 * phys_const->const_boltzmann_k /
                   phys_const->const_electron_volt;

    /* Set photon density and compute radiation density if necessary */
    if (c->Omega_r != 0.) {

      c->Omega_g = c->Omega_r;
      c->Omega_ur = 0.;

    } else {

      /* Infer CMB density from the temperature */
      c->Omega_g = pow(c->T_CMB_0, 4) / rho_c3_on_4sigma;

      /* Compute the density of ultra-relativistic fermionic species */
      c->Omega_ur = c->N_ur * (7. / 8.) * dec_4 * c->Omega_g;

      /* Compute the total radiation density */
      c->Omega_r = c->Omega_g + c->Omega_ur;
    }

    /* Compute effective number of relativistic species at early times */
    double N_nu_tot_deg = 0.;
    for (int i = 0; i < c->N_nu; i++) {
      N_nu_tot_deg += c->deg_nu[i];
    }
    c->N_eff = c->N_ur + N_nu_tot_deg * pow(c->T_nu_0 / c->T_CMB_0, 4) / dec_4;
    c->deg_nu_tot = N_nu_tot_deg;

    /* Initialise the neutrino density interpolation tables if necessary */
    c->neutrino_density_early_table = NULL;
    c->neutrino_density_late_table = NULL;
    cosmology_init_neutrino_tables(c);

    /* Retrieve the present-day total density due to massive neutrinos */
    c->Omega_nu_0 = cosmology_get_neutrino_density(c, 1);
    c->Omega_nu = c->Omega_nu_0;  // will be updated
  }

  /* Curvature density (for closure) */
  const double Omega_m = c->Omega_cdm + c->Omega_b;
  c->Omega_k = 1. - (Omega_m + c->Omega_r + c->Omega_lambda + c->Omega_nu_0);

  /* Initialise the interpolation tables */
  c->drift_fac_interp_table = NULL;
  c->grav_kick_fac_interp_table = NULL;
  c->hydro_kick_fac_interp_table = NULL;
  c->time_interp_table = NULL;
  c->time_interp_table_offset = 0.;
  cosmology_init_tables(c);

  /* Set remaining variables to valid values */
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

  c->Omega_cdm = 0.;
  c->Omega_r = 0.;
  c->Omega_nu = 0.;
  c->Omega_k = 0.;
  c->Omega_lambda = 0.;
  c->Omega_b = 0.;
  c->w_0 = 0.;
  c->w_a = 0.;
  c->h = 1.;
  c->w = -1.;

  c->Omega_ur = 0.;
  c->Omega_g = 0.;
  c->T_nu_0 = 0.;
  c->T_nu_0_eV = 0.;
  c->N_nu = 0;
  c->N_ur = 0.;
  c->N_eff = 0.;
  c->deg_nu_tot = 0.;

  c->a_begin = 1.;
  c->a_end = 1.;
  c->log_a_begin = 0.;
  c->log_a_end = 0.;
  c->log_a_long_begin = 0.;
  c->log_a_long_mid = 0.;
  c->log_a_long_end = 0.;

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
  c->mean_density = 0.;
  c->mean_density_Omega_m = 0;
  c->mean_density_Omega_b = 0;
  c->overdensity_BN98 = 0.;
  c->T_CMB_0 = 0.;
  c->T_CMB_0_K = 0.;

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
  c->neutrino_density_early_table = NULL;
  c->neutrino_density_late_table = NULL;
  c->time_interp_table_offset = 0.;
  c->scale_factor_interp_table = NULL;
  c->comoving_distance_interp_table = NULL;
  c->comoving_distance_inverse_interp_table = NULL;

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
                                  const integertime_t ti_start,
                                  const integertime_t ti_end) {

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
                                      const integertime_t ti_start,
                                      const integertime_t ti_end) {

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
                                       const integertime_t ti_start,
                                       const integertime_t ti_end) {

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
                                      const integertime_t ti_start,
                                      const integertime_t ti_end) {

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
                                       const integertime_t ti_start,
                                       const integertime_t ti_end) {

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
                                const integertime_t ti_start,
                                const integertime_t ti_end) {

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
 * @brief Compute the comoving distance to the specified scale factor
 *
 * @param c The current #cosmology.
 * @param a The scale factor
 */
double cosmology_get_comoving_distance(const struct cosmology *c,
                                       const double a) {

#ifdef SWIFT_DEBUG_CHECKS
  if (a < c->a_begin) error("a must be >= a_begin");
  if (a > c->a_end) error("a must be <= a_end");
#endif

  const double log_a = log(a);

  /* Comoving distance from a_begin to a */
  const double dist = interp_table(c->comoving_distance_interp_table, log_a,
                                   c->log_a_begin, c->log_a_end);

  /* Subtract dist from comoving distance from a_begin to a=1 */
  return c->comoving_distance_interp_table_offset - dist;
}

/**
 * @brief Compute scale factor from a comoving distance (in internal units).
 *
 * @param c The current #cosmology.
 * @param r The comoving distance
 * @return The scale factor.
 */
double cosmology_scale_factor_at_comoving_distance(const struct cosmology *c,
                                                   double r) {

  /* Get comoving distance from a_begin to a corresponding to input r */
  const double r_interp = c->comoving_distance_interp_table_offset - r;

  const double a =
      interp_table(c->comoving_distance_inverse_interp_table, r_interp, 0.0,
                   c->comoving_distance_start_to_end);
  return a + c->a_begin;
}

/**
 * @brief Compute neutrino density parameter Omega_nu at the given scale-factor
 * This is the effective present day value, i.e. must be multiplied by (1+z)^4
 *
 * @param c The current #cosmology.
 * @param a The scale factor
 * @return The density parameter
 */
double cosmology_get_neutrino_density(const struct cosmology *c, double a) {

  if (c->N_nu == 0) return 0.;

  const double log_a = log(a);

  if (log_a < c->log_a_long_begin)
    return c->neutrino_density_early_table[0];
  else if (log_a < c->log_a_long_mid)
    return interp_table(c->neutrino_density_early_table, log_a,
                        c->log_a_long_begin, c->log_a_long_mid);
  else
    return interp_table(c->neutrino_density_late_table, log_a,
                        c->log_a_long_mid, c->log_a_long_end);
}

/**
 * @brief Compute the cosmic time (in internal units) between two scale factors
 *
 * @param c The current #cosmology.
 * @param a_start the starting scale factor
 * @param a_end the ending scale factor
 */
double cosmology_get_delta_time_from_scale_factors(const struct cosmology *c,
                                                   const double a_start,
                                                   const double a_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (a_end < a_start) error("a_end must be >= a_start");
  if (a_end < c->a_begin) error("Error a_end can't be smaller than a_begin");
#endif

  const double log_a_start = log(a_start);
  const double log_a_end = log(a_end);

  /* Time between a_begin and a_start */
  const double t1 = interp_table(c->time_interp_table, log_a_start,
                                 c->log_a_begin, c->log_a_end);

  /* Time between a_begin and a_end */
  const double t2 = interp_table(c->time_interp_table, log_a_end,
                                 c->log_a_begin, c->log_a_end);

  return t2 - t1;
}

/**
 * @brief Compute the time corresponding to the timebase interval at the current
 * redshift
 *
 * This function is slow as it performs the actual integral.
 *
 * @param c The comology model.
 * @param ti_current The current point on the time-line.
 *
 * @return The time corresponding to c->time_base in internal time units.
 */
double cosmology_get_timebase(struct cosmology *c,
                              const integertime_t ti_current) {

#ifdef HAVE_LIBGSL

  /* We are going to average over a number of time-bins
   * as in some cases the smallest bin is too small for double accuracy */
  const int num_bins = 24;

  /* Range in log(a) */
  const double log_a_start = c->log_a_begin + (ti_current)*c->time_base;
  const double log_a_end =
      c->log_a_begin + (ti_current + (1LL << num_bins)) * c->time_base;

  /* Range in (a) */
  const double a_start = exp(log_a_start);
  const double a_end = exp(log_a_end);

  /* Initalise the GSL workspace and function */
  gsl_integration_workspace *space =
      gsl_integration_workspace_alloc(GSL_workspace_size);
  gsl_function F = {&time_integrand, c};

  /* Perform the integral */
  double result, abserr;
  gsl_integration_qag(&F, a_start, a_end, 0, 1.0e-10, GSL_workspace_size,
                      GSL_INTEG_GAUSS61, space, &result, &abserr);

  result /= (1LL << num_bins);

  /* Free the workspace */
  gsl_integration_workspace_free(space);

  return result;

#else
  error("Code not compiled with GSL. Can't compute cosmology integrals.");
  return 0.;
#endif
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
                   c->time_interp_table_max);
  return a + c->a_begin;
}

/**
 * @brief Prints the #cosmology model to stdout.
 */
void cosmology_print(const struct cosmology *c) {

  const double Omega_m = c->Omega_cdm + c->Omega_b;
  message(
      "Density parameters: [O_m, O_l, O_b, O_k, O_r] = [%f, %f, %f, %f, %f]",
      Omega_m, c->Omega_lambda, c->Omega_b, c->Omega_k, c->Omega_r);
  message(
      "Additional density parameters: [O_nu_0, O_cdm, O_ur, O_g] = [%f, "
      "%f, %f, %f]",
      c->Omega_nu_0, c->Omega_cdm, c->Omega_ur, c->Omega_g);
  message("Dark energy equation of state: w_0=%f w_a=%f", c->w_0, c->w_a);
  message("Hubble constant: h = %f, H_0 = %e U_t^(-1)", c->h, c->H0);
  message("Hubble time: 1/H0 = %e U_t", c->Hubble_time);
  message("CMB temperature at z=0 implied by cosmology: T_CMB = %e U_T",
          c->T_CMB_0);
  message("Neutrino temperature at z=0: T_nu = %e U_T", c->T_nu_0);
  message("Numbers of relatistic species: [N_nu, N_ur, N_eff] = [%d, %f, %f]",
          c->N_nu, c->N_ur, c->N_eff);
  /* Print neutrino masses and degeneracies */
  if (c->N_nu > 0) {
    char neutrino_mass_string[10 * c->N_nu];
    char neutrino_deg_string[10 * c->N_nu];
    for (int i = 0; i < c->N_nu; i++) {
      sprintf(neutrino_mass_string + i * 10, "%.2e  ", c->M_nu_eV[i]);
      sprintf(neutrino_deg_string + i * 10, "%.2e  ", c->deg_nu[i]);
    }
    message("Neutrino masses: %seV", neutrino_mass_string);
    message("Neutrino degeneracies: %s", neutrino_deg_string);
  }
  message("Universe age at present day: %e U_t",
          c->universe_age_at_present_day);
}

void cosmology_clean(struct cosmology *c) {

  swift_free("cosmo.table", c->drift_fac_interp_table);
  swift_free("cosmo.table", c->grav_kick_fac_interp_table);
  swift_free("cosmo.table", c->hydro_kick_fac_interp_table);
  swift_free("cosmo.table", c->hydro_kick_corr_interp_table);
  swift_free("cosmo.table", c->time_interp_table);
  swift_free("cosmo.table", c->scale_factor_interp_table);
  swift_free("cosmo.table", c->comoving_distance_interp_table);
  swift_free("cosmo.table", c->comoving_distance_inverse_interp_table);
  if (c->N_nu > 0) {
    swift_free("cosmo.table", c->neutrino_density_early_table);
    swift_free("cosmo.table", c->neutrino_density_late_table);
    swift_free("Mnu", c->M_nu_eV);
    swift_free("degnu", c->deg_nu);
  }
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
  io_write_attribute_d(h_grp, "Omega_m", c->Omega_cdm + c->Omega_b);
  io_write_attribute_d(h_grp, "Omega_r", c->Omega_r);
  io_write_attribute_d(h_grp, "Omega_b", c->Omega_b);
  io_write_attribute_d(h_grp, "Omega_k", c->Omega_k);
  io_write_attribute_d(h_grp, "Omega_lambda", c->Omega_lambda);
  io_write_attribute_d(h_grp, "Omega_nu", c->Omega_nu);
  io_write_attribute_d(h_grp, "Omega_nu_0", c->Omega_nu_0);
  io_write_attribute_d(h_grp, "Omega_ur", c->Omega_ur);
  io_write_attribute_d(h_grp, "Omega_cdm", c->Omega_cdm);
  io_write_attribute_d(h_grp, "Omega_g", c->Omega_g);
  io_write_attribute_d(h_grp, "T_nu_0 [internal units]", c->T_nu_0);
  io_write_attribute_d(h_grp, "T_nu_0 [eV]", c->T_nu_0_eV);
  io_write_attribute_d(h_grp, "N_eff", c->N_eff);
  io_write_attribute_d(h_grp, "N_ur", c->N_ur);
  io_write_attribute_i(h_grp, "N_nu", c->N_nu);
  if (c->N_nu > 0) {
    io_write_attribute(h_grp, "M_nu_eV", DOUBLE, c->M_nu_eV, c->N_nu);
    io_write_attribute(h_grp, "deg_nu", DOUBLE, c->deg_nu, c->N_nu);
  }
  io_write_attribute_d(h_grp, "deg_nu_tot", c->deg_nu_tot);
  io_write_attribute_d(h_grp, "T_CMB_0 [internal units]", c->T_CMB_0);
  io_write_attribute_d(h_grp, "T_CMB_0 [K]", c->T_CMB_0_K);
  io_write_attribute_d(h_grp, "w_0", c->w_0);
  io_write_attribute_d(h_grp, "w_a", c->w_a);
  io_write_attribute_d(h_grp, "w", c->w);
  io_write_attribute_d(h_grp, "Redshift", c->z);
  io_write_attribute_d(h_grp, "Scale-factor", c->a);
  io_write_attribute_d(h_grp, "Critical density [internal units]",
                       c->critical_density);
  io_write_attribute_d(h_grp,
                       "Critical density at redshift zero [internal units]",
                       c->critical_density_0);
  io_write_attribute_d(h_grp, "Mean matter density [internal units]",
                       c->mean_density_Omega_m);
  io_write_attribute_d(h_grp, "Mean baryonic density [internal units]",
                       c->mean_density_Omega_b);
  io_write_attribute_d(h_grp, "Virial overdensity (BN98)", c->overdensity_BN98);
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

  /* Also store the neutrino mass and degeneracy arrays if necessary */
  if (cosmology->N_nu > 0) {
    restart_write_blocks((double *)cosmology->M_nu_eV, sizeof(double),
                         cosmology->N_nu, stream, "cosmology->M_nu_eV",
                         "neutrino masses eV");
    restart_write_blocks((double *)cosmology->deg_nu, sizeof(double),
                         cosmology->N_nu, stream, "cosmology->deg_nu",
                         "neutrino degeneracies");
  }
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

  /* Restore the neutrino mass and degeneracy arrays if necessary */
  if (cosmology->N_nu > 0) {
    cosmology->M_nu_eV =
        (double *)swift_malloc("Mnu", cosmology->N_nu * sizeof(double));
    restart_read_blocks((double *)cosmology->M_nu_eV, sizeof(double),
                        cosmology->N_nu, stream, NULL, "neutrino masses eV");
    cosmology->deg_nu =
        (double *)swift_malloc("degnu", cosmology->N_nu * sizeof(double));
    restart_read_blocks((double *)cosmology->deg_nu, sizeof(double),
                        cosmology->N_nu, stream, NULL, "neutrino degeneracies");
  }

  /* Re-initialise the tables if using a cosmology. */
  if (enabled) {
    cosmology_init_neutrino_tables(cosmology);
    cosmology_init_tables(cosmology);
  }
}
