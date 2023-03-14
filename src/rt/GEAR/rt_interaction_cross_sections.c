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

#include "rt_interaction_cross_sections.h"

#include "error.h"
#include "rt_blackbody.h"
#include "rt_properties.h"
#include "rt_species.h"

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
  /* Keep this function in the .c file so you don't have to
   * include the spectra .h files like rt_blackbody.h anywhere
   * else */

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
static double spectrum_integrand(double nu, void *params) {
  return rt_interaction_rates_get_spectrum(nu, params);
}

/**
 * Spectrum function divided by photon energy h*nu to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
static double spectrum_over_hnu_integrand(double nu, void *params) {
  struct rt_spectrum_integration_params *p =
      (struct rt_spectrum_integration_params *)params;
  const double E = nu * p->h_planck;
  const double J = rt_interaction_rates_get_spectrum(nu, params);
  if (E > 0.) {
    return J / E;
  } else {
    return 0.;
  }
}

/**
 * Spectrum times cross section function to be integrated.
 * This function is called by the GSL integrator.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
static double spectrum_times_sigma_integrand(double nu, void *params) {
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
static double spectrum_times_sigma_over_hnu_integrand(double nu, void *params) {
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
static double rt_cross_sections_integrate_gsl(
    double (*function)(double, void *), double nu_start, double nu_stop,
    int npoints, struct rt_spectrum_integration_params *params) {

  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npoints);
  double result, error;

  F.function = function;
  F.params = (void *)params;
  /* NOTE: there is an option to use the integrator with an upper limit
   * of infinity, but this is accurate enough for now when setting a
   * high enough maximal frequency. */
  gsl_integration_qags(&F, nu_start, nu_stop, /*espabs=*/0., /*epsrel=*/1e-7,
                       npoints, w, &result, &error);

  /* Clean up after yourself. */
  gsl_integration_workspace_free(w);

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
void rt_cross_sections_init(struct rt_props *restrict rt_props,
                            const struct phys_const *restrict phys_const,
                            const struct unit_system *restrict us) {

  /* Allocate the space to store the (cross section) integrals */
  /* --------------------------------------------------------- */

  double **cse = malloc(RT_NGROUPS * sizeof(double *));
  double **csn = malloc(RT_NGROUPS * sizeof(double *));
  double *av_energy = rt_props->average_photon_energy;
  double *photon_number_integral = rt_props->photon_number_integral;
  for (int group = 0; group < RT_NGROUPS; group++) {
    cse[group] = malloc(rt_ionizing_species_count * sizeof(double));
    csn[group] = malloc(rt_ionizing_species_count * sizeof(double));
    av_energy[group] = 0.;
    photon_number_integral[group] = 0.;
  }

  double integral_E[RT_NGROUPS];
  double integral_N[RT_NGROUPS];
  double integral_sigma_E[RT_NGROUPS][rt_ionizing_species_count];
  double integral_sigma_E_over_hnu[RT_NGROUPS][rt_ionizing_species_count];

  /* Grab constants and conversions in cgs */
  /* ------------------------------------- */
  const double T_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TEMPERATURE);

  const float dimension_kB[5] = {1, 2, -2, 0, -1}; /* [g cm^2 s^-2 K^-1] */
  const double kB_to_cgs =
      units_general_cgs_conversion_factor(us, dimension_kB);
  const double kB_cgs = phys_const->const_boltzmann_k * kB_to_cgs;

  const float dimension_h[5] = {1, 2, -1, 0, 0}; /* [g cm^2 s^-1] */
  const double h_to_cgs = units_general_cgs_conversion_factor(us, dimension_h);
  const double h_planck_cgs = phys_const->const_planck_h * h_to_cgs;

  const double c_cgs = phys_const->const_speed_light_c *
                       units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);

  /* Prepare parameter struct for integration functions */
  /* -------------------------------------------------- */
  const int spectype = rt_props->stellar_spectrum_type;
  const double maxfreq_const_spectrum =
      rt_props->const_stellar_spectrum_max_frequency;
  const double T_bb = rt_props->stellar_spectrum_blackbody_T * T_cgs;
  struct rt_photoion_cs_parameters cs_params = rt_init_photoion_cs_params_cgs();

  struct rt_spectrum_integration_params integration_params = {
      /*species=*/0,
      /*spectrum_type=*/spectype,
      /*freq_max for const spectrum=*/maxfreq_const_spectrum,
      /*T=*/T_bb,
      /*kB=*/kB_cgs,
      /*h_planck=*/h_planck_cgs,
      /*c=*/c_cgs,
      /*cross section params=*/&cs_params};

  /* Set up Integration Limits */
  /* ------------------------- */

  /* Get start and end points of the integrals */
  double nu_start[RT_NGROUPS];
  double nu_stop[RT_NGROUPS];
  for (int group = 0; group < RT_NGROUPS; group++)
    nu_start[group] = rt_props->photon_groups[group];
  for (int group = 0; group < RT_NGROUPS - 1; group++)
    nu_stop[group] = rt_props->photon_groups[group + 1];

  if (RT_NGROUPS == 1) {
    /* If we only have one group, start integrating from the Hydrogen
     * ionization frequency, not from zero. The reasoning here is that
     * typically you define the *ionizing* radiation as stellar emission
     * rates, not the *total* radiation. */
    nu_start[0] = cs_params.E_ion[rt_ionizing_species_HI] / h_planck_cgs;
    if (engine_rank == 0)
      message(
          "WARNING: with only 1 photon group, I'll start integrating"
          " the cross sections at the first ionizing frequency %.3g",
          nu_start[0]);
  } else {
    /* don't start at exactly 0 to avoid unlucky divisions */
    if (nu_start[0] == 0.) nu_start[0] = min(1e-20, 1e-12 * nu_start[1]);
  }

  /* Get frequency at which we stop integrating */
  double nu_stop_final;
  if (rt_props->stellar_spectrum_type == 0) {
    nu_stop_final = rt_props->const_stellar_spectrum_max_frequency;
  } else if (rt_props->stellar_spectrum_type == 1) {
    nu_stop_final = 10. * blackbody_peak_frequency(T_bb, kB_cgs, h_planck_cgs);
  } else {
    nu_stop_final = -1.;
    error("Unknown stellar spectrum type %d", rt_props->stellar_spectrum_type);
  }
  nu_stop[RT_NGROUPS - 1] = nu_stop_final;

  /* Compute Integrals */
  /* ----------------- */
  for (int g = 0; g < RT_NGROUPS; g++) {
    /* This is independent of species. */
    integral_E[g] = rt_cross_sections_integrate_gsl(
        spectrum_integrand, nu_start[g], nu_stop[g], RT_INTEGRAL_NPOINTS,
        &integration_params);
    integral_N[g] = rt_cross_sections_integrate_gsl(
        spectrum_over_hnu_integrand, nu_start[g], nu_stop[g],
        RT_INTEGRAL_NPOINTS, &integration_params);

    for (int s = 0; s < rt_ionizing_species_count; s++) {
      integration_params.species = s;
      integral_sigma_E[g][s] = rt_cross_sections_integrate_gsl(
          spectrum_times_sigma_integrand, nu_start[g], nu_stop[g],
          RT_INTEGRAL_NPOINTS, &integration_params);
      integral_sigma_E_over_hnu[g][s] = rt_cross_sections_integrate_gsl(
          spectrum_times_sigma_over_hnu_integrand, nu_start[g], nu_stop[g],
          RT_INTEGRAL_NPOINTS, &integration_params);
    }
  }

  /* Now compute the actual average cross sections */
  /* --------------------------------------------- */
  for (int g = 0; g < RT_NGROUPS; g++) {
    photon_number_integral[g] = integral_N[g];
    if (integral_N[g] > 0.) {
      av_energy[g] = integral_E[g] / integral_N[g];
    } else {
      av_energy[g] = 0.;
    }
    for (int s = 0; s < rt_ionizing_species_count; s++) {
      if (integral_E[g] > 0.) {
        cse[g][s] = integral_sigma_E[g][s] / integral_E[g];
        csn[g][s] = integral_sigma_E_over_hnu[g][s] / integral_N[g];
      } else {
        /* No radiation = no interaction */
        cse[g][s] = 0.;
        csn[g][s] = 0.;
      }
    }
  }

  /* for (int g = 0; g < RT_NGROUPS; g++) { */
  /*   printf("\nGroup %d\n", g); */
  /*   printf("nu_start:                  %12.6g\n", nu_start[g]); */
  /*   printf("nu_end:                    %12.6g\n", nu_stop[g]); */
  /*   printf("spectrum energy integral:  %12.6g\n", integral_E[g]); */
  /*   printf("spectrum number integral:  %12.6g\n", integral_N[g]); */
  /*   printf("average photon energy:     %12.6g\n", av_energy[g]); */
  /*  */
  /*   printf("species:                   "); */
  /*   for (int s = 0; s < rt_ionizing_species_count; s++) printf("%12d ", s);
   */
  /*   printf("\n"); */
  /*  */
  /*   printf("integral sigma * E:        "); */
  /*   for (int s = 0; s < rt_ionizing_species_count; s++) printf("%12.6g ",  */
  /*       integral_sigma_E[g][s]); */
  /*   printf("\n"); */
  /*  */
  /*   printf("integral sigma * E / h nu: "); */
  /*   for (int s = 0; s < rt_ionizing_species_count; s++) printf("%12.6g ",  */
  /*       integral_sigma_E_over_hnu[g][s]); */
  /*   printf("\n"); */
  /*  */
  /*   printf("energy weighted c.section: "); */
  /*   for (int s = 0; s < rt_ionizing_species_count; s++) printf("%12.6g ",  */
  /*       cse[g][s]); */
  /*   printf("\n"); */
  /*  */
  /*   printf("number weighted c.section: "); */
  /*   for (int s = 0; s < rt_ionizing_species_count; s++) printf("%12.6g ",  */
  /*       csn[g][s]); */
  /*   printf("\n"); */
  /* } */

  /* Store the results */
  /* ----------------- */

  if (rt_props->energy_weighted_cross_sections != NULL) {
    for (int g = 0; g < RT_NGROUPS; g++) {
      if (rt_props->energy_weighted_cross_sections[g] != NULL)
        free(rt_props->energy_weighted_cross_sections[g]);
    }
    free(rt_props->energy_weighted_cross_sections);
  }
  rt_props->energy_weighted_cross_sections = cse;

  if (rt_props->number_weighted_cross_sections != NULL) {
    for (int g = 0; g < RT_NGROUPS; g++) {
      if (rt_props->number_weighted_cross_sections[g] != NULL)
        free(rt_props->number_weighted_cross_sections[g]);
    }
    free(rt_props->number_weighted_cross_sections);
  }
  rt_props->number_weighted_cross_sections = csn;
}
