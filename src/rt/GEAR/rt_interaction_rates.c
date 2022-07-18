/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#include "rt_interaction_rates.h"

#include "rt_blackbody.h"
#include "rt_properties.h"

#include <gsl/gsl_integration.h>

/**
 * @brief compute the chosen spectrum for a given frequency nu.
 * This function is intended to be used to integrate the spectra
 * over a frequency range in order to obtain averaged ionization
 * cross sections.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
__attribute__((always_inline)) INLINE static double
rt_interaction_rates_get_spectrum(const double nu, void *params) {

  struct rt_spectrum_integration_params *pars =
      (struct rt_spectrum_integration_params *)params;

  if (pars->spectrum_type == 0) {
    /* Constant spectrum */
    if (nu <= pars->const_stellar_spectrum_frequency_max) {
      return 1.;
    } else {
      return 0.;
    }
  } else if (pars->spectrum_type == 1) {
    /* Blackbody spectrum */
    const double T = pars->T;
    const double kB = pars->kB;
    const double h_planck = pars->h_planck;
    const double c = pars->c;
    return blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
  } else {
    error("Unknown stellar spectrum type selected: %d", pars->spectrum_type);
    return 0.;
  }
}

/**
 * Spectrum function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_integrand(double nu, void *params) {
  return rt_interaction_rates_get_spectrum(nu, params);
}

/**
 * Spectrum times cross section function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_times_sigma_integrand(double nu, void *params) {
  struct rt_spectrum_integration_params *p =
      (struct rt_spectrum_integration_params *)params;
  const double E = nu * p->h_planck;
  const double sigma =
      photoionization_cross_section(E, p->species, p->cs_params);
  const double J = rt_interaction_rates_get_spectrum(nu, params);
  return J * sigma;
}

/**
 * Spectrum times cross section divided by h*nu function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
double spectrum_times_sigma_over_hnu_integrand(double nu, void *params) {
  struct rt_spectrum_integration_params *p =
      (struct rt_spectrum_integration_params *)params;
  const double E = nu * p->h_planck;
  const double sigma =
      photoionization_cross_section(E, p->species, p->cs_params);
  const double J = rt_interaction_rates_get_spectrum(nu, params);
  return J * sigma / E;
}

/**
 * Integrate a function from nu_start to nu_stop with GSL routines
 *
 * @param function function to integrate
 * @param nu_start lower boundary of the integral
 * @param nu_stop upper boundary of the integral
 * @param params spectrum integration params.
 * */
double rt_interaction_rates_integrate_gsl(
    double (*function)(double, void *), double nu_start, double nu_stop,
    int npoints, struct rt_spectrum_integration_params *params) {

  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npoints);
  double result, error;

  F.function = function;
  F.params = (void *)params;
  gsl_integration_qags(&F, nu_start, nu_stop, /*espabs=*/0., /*epsrel=*/1e-7,
                       npoints, w, &result, &error);

  return result;
}

/**
 * @brief allocate and pre-compute the averaged cross sections
 * for each photon group and ionizing species.
 *
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param us internal units struct
 **/
void rt_interaction_rates_init(struct rt_props *restrict rt_props,
                               const struct phys_const *restrict phys_const,
                               const struct unit_system *restrict us) {

  /* Allocate the space to store the (cross section) integrals */
  /* --------------------------------------------------------- */
  double **cse = malloc(RT_NGROUPS * sizeof(double *));
  double **csn = malloc(RT_NGROUPS * sizeof(double *));
  for (int group = 0; group < RT_NGROUPS; group++) {
    cse[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
    csn[group] = malloc(RT_NIONIZING_SPECIES * sizeof(double));
  }

  double integral_E[RT_NGROUPS];
  double integral_sigma_E[RT_NGROUPS][RT_NIONIZING_SPECIES];
  double integral_sigma_E_over_hnu[RT_NGROUPS][RT_NIONIZING_SPECIES];

  /* Prepare parameter struct for integration functions */
  /* -------------------------------------------------- */
  struct rt_spectrum_integration_params integration_params = {0,  0,  0., 0.,
                                                              0., 0., 0., NULL};
  integration_params.spectrum_type = rt_props->stellar_spectrum_type;
  integration_params.const_stellar_spectrum_frequency_max =
      rt_props->const_stellar_spectrum_max_frequency;

  /* Grab constants and conversions in cgs */
  const double T_cgs = rt_props->stellar_spectrum_blackbody_T *
                       units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);
  integration_params.T = T_cgs;

  const float dimension_kB[5] = {1, 2, -2, 0, -1}; /* [g cm^2 s^-2 K^-1] */
  const double kB_to_cgs =
      units_general_cgs_conversion_factor(us, dimension_kB);
  const double kB_cgs = phys_const->const_boltzmann_k * kB_to_cgs;
  integration_params.kB = kB_cgs;

  const float dimension_h[5] = {1, 2, -1, 0, 0}; /* [g cm^2 s^-1] */
  const double h_to_cgs = units_general_cgs_conversion_factor(us, dimension_h);
  const double h_planck_cgs = phys_const->const_planck_h * h_to_cgs;
  integration_params.h_planck = h_planck_cgs;

  const double c_cgs = phys_const->const_speed_light_c *
                       units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
  integration_params.c = c_cgs;

  struct rt_photoion_cs_parameters cs_params = rt_init_photoion_cs_params_cgs();
  integration_params.cs_params = &cs_params;

  /* Compute integrals */
  /* ----------------- */

  /* Get start and end points of the integrals */
  double nu_stop_final;
  if (rt_props->stellar_spectrum_type == 0) {
    nu_stop_final = rt_props->const_stellar_spectrum_max_frequency;
  } else if (rt_props->stellar_spectrum_type == 1) {
    nu_stop_final = 10. * blackbody_peak_frequency(integration_params.T,
                                                   integration_params.kB,
                                                   integration_params.h_planck);
  } else {
    nu_stop_final = -1.;
    error("Unknown stellar spectrum type %d", rt_props->stellar_spectrum_type);
  }

  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++)
    nu_start[group] = rt_props->photon_groups[group];
  for (int group = 0; group < RT_NGROUPS - 1; group++)
    nu_stop[group] = rt_props->photon_groups[group + 1];

  nu_start[0] = 1e-20; /* don't start at exactly 0 to avoid unlucky divisions */
  nu_stop[RT_NGROUPS - 1] = nu_stop_final;

  for (int group = 0; group < RT_NGROUPS; group++) {
    /* This is independent of species. */
    integral_E[group] = rt_interaction_rates_integrate_gsl(
        spectrum_integrand, nu_start[group], nu_stop[group],
        RT_INTEGRAL_NPOINTS, &integration_params);

    for (int species = 0; species < RT_NIONIZING_SPECIES; species++) {
      integration_params.species = species;
      integral_sigma_E[group][species] = rt_interaction_rates_integrate_gsl(
          spectrum_times_sigma_integrand, nu_start[group], nu_stop[group],
          RT_INTEGRAL_NPOINTS, &integration_params);
      integral_sigma_E_over_hnu[group][species] =
          rt_interaction_rates_integrate_gsl(
              spectrum_times_sigma_over_hnu_integrand, nu_start[group],
              nu_stop[group], RT_INTEGRAL_NPOINTS, &integration_params);
    }
  }

  /* Now compute the actual average cross sections */
  /* --------------------------------------------- */
  for (int group = 0; group < RT_NGROUPS; group++) {
    for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
      if (integral_E[group] > 0.) {
        cse[group][spec] = integral_sigma_E[group][spec] / integral_E[group];
        csn[group][spec] =
            integral_sigma_E_over_hnu[group][spec] / integral_E[group];
      } else {
        /* No radiation = no interaction */
        cse[group][spec] = 0.;
        csn[group][spec] = 0.;
      }
    }
  }

  /* Store the results */
  /* ----------------- */
  rt_props->energy_weighted_cross_sections = cse;
  rt_props->number_weighted_cross_sections = csn;
}
