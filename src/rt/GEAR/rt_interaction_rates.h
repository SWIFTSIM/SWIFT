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
#ifndef SWIFT_RT_GEAR_INTERACION_RATES_H
#define SWIFT_RT_GEAR_INTERACION_RATES_H

#include "rt_blackbody.h"
#include "rt_parameters.h"
#include "rt_properties.h"

#include <gsl/gsl_integration.h>

#define RT_INTEGRAL_NPOINTS 1000

/* HI, HeI, HeII */
#define RT_NIONIZING_SPECIES 3

/*! Struct containing the parametrized cross section parameters
 * for each photoionizing species.
 * Correct usage is to call rt_init_photoion_cs_params_cgs(),
 * which returns a fully initialized struct. */
struct rt_photoion_cs_parameters {
  double E_ion[RT_NIONIZING_SPECIES];
  double E_zero[RT_NIONIZING_SPECIES];
  double sigma_zero[RT_NIONIZING_SPECIES];
  double P[RT_NIONIZING_SPECIES];
  double ya[RT_NIONIZING_SPECIES];
  double yw[RT_NIONIZING_SPECIES];
  double y0[RT_NIONIZING_SPECIES];
  double y1[RT_NIONIZING_SPECIES];
};

/*! Parameters needed to compute the stellar spectra integrals -
 * the data needs to be encapsulated in a single struct to be passed
 * as an argument for GSL integrators. */
struct rt_spectrum_integration_params {
  /* Which species are we dealing with? */
  int species;
  /* Which spectrum type are we dealing with? */
  int spectrum_type;
  /* Max frequency for const spectrum */
  double const_stellar_spectrum_frequency_max;
  /* Temperature of blackbody in correct units */
  double T;
  /* Boltzmann constant in correct units */
  double kB;
  /* Planck constant in correct units */
  double h_planck;
  /* speed of light in correct units */
  double c;
  /* Values for the cross section parametrization */
  struct rt_photoion_cs_parameters *cs_params;
};

/**
 * Initialize the parameters for the cross section computation in cgs,
 * and return a fully and correctly initialized struct.
 * The data is taken from Verner et al. 1996
 * (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V) via Rosdahl et al. 2013
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
 */
static struct rt_photoion_cs_parameters rt_init_photoion_cs_params_cgs(void) {

  struct rt_photoion_cs_parameters photoion_cs_params_cgs = {
      /* E_ion =         {13.60,     24.59,     54.42}          eV */
      /* E_zero =        {0.4298,    0.1361,    1.720},         eV */
      /* E_ion =      */ {2.179e-11, 3.940e-11, 8.719e-11},  /* erg */
      /* E_zero =     */ {6.886e-13, 2.181 - 13, 2.756e-12}, /* erg */
      /* sigma_zero = */ {5.475e-14, 9.492e-16, 1.369e-14},  /* cm^-2 */
      /* P =          */ {2.963, 3.188, 2.963},
      /* ya =         */ {32.88, 1.469, 32.88},
      /* yw =         */ {0., 2.039, 0.},
      /* y0 =         */ {0., 0.4434, 0.},
      /* y1 =         */ {0., 2.136, 0.}};

  return photoion_cs_params_cgs;
}

/**
 * Compute the parametrized cross section for a given energy and species.
 * The parametrization is taken from Verner et al. 1996
 * (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V) via Rosdahl et al. 2013
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
 *
 * @param E energy for which to compute the cross section.
 * @param species index of species, 0 < species < RT_NIONIZING_SPECIES
 * @param params cross section parameters struct
 */
static double photoionization_cross_section(
    const double E, const int species,
    const struct rt_photoion_cs_parameters *params) {

  const double E0 = params->E_zero[species];
  const double E_ion = params->E_ion[species];
  const double y0 = params->y0[species];
  const double y1 = params->y1[species];
  const double yw = params->yw[species];
  const double ya = params->ya[species];
  const double P = params->P[species];
  const double sigma_0 = params->sigma_zero[species];

  if (E < E_ion) return 0.;

  const double x = E / E0 - y0;
  const double y = sqrt(x * x + y1 * y1);
  const double temp1 = pow(y, 0.5 * P - 5.5);
  const double temp2 = pow(1. + sqrt(y / ya), -P);

  return sigma_0 * ((x - 1.) * (x - 1.) + yw * yw) * temp1 * temp2;
}

/**
 * @brief compute the chosen spectrum for a given frequency nu.
 * This function is indended to be used to integrate the spectra
 * over a frequency range in order to obtain averaged ionization
 * cross sections.
 *
 * @param nu frequency at which to evaluate spectrum
 * @param params spectrum integration params. Needs to be of type
 * void* for GSL integrators.
 */
static double rt_interaction_rates_get_spectrum(const double nu, void *params) {

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
static double rt_interaction_rates_integrate_gsl(
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
static void rt_interaction_rates_init(struct rt_props *restrict rt_props,
                                      const struct phys_const *restrict
                                          phys_const,
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

/**
 * @brief compute the heating, ionization, and dissassociation rates
 * for the particle radiation field as needed by grackle, and the
 * net absorption/emission rates for each photon group
 *
 * @param rates (return) Interaction rates for grackle. [0]: heating rate.
 * [1]: HI ionization. [2]: HeI ionization. [3]: HeII ionization.
 * [4]: H2 dissociation.
 * @param heating_rates_by_group (return) net absorption/emission rates of each
 * photon frequency group in internal units.
 * @param p particle to work on
 * @param species_densities the physical densities of all traced species
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param us internal units struct
 * @param cosmo cosmology struct
 **/
__attribute__((always_inline)) INLINE static void
rt_tchem_get_interaction_rates(gr_float rates[5],
                               float heating_rates_by_group[RT_NGROUPS],
                               const struct part *restrict p,
                               gr_float species_densities[6],
                               const struct rt_props *restrict rt_props,
                               const struct phys_const *restrict phys_const,
                               const struct unit_system *restrict us,
                               const struct cosmology *restrict cosmo) {

  rates[0] = 0.; /* Needs to be in [erg / s / cm^3 / nHI] for grackle. */
  rates[1] = 0.; /* [1 / time_units] */
  rates[2] = 0.; /* [1 / time_units] */
  rates[3] = 0.; /* [1 / time_units] */
  rates[4] = 0.; /* [1 / time_units] */
  for (int group = 0; group < RT_NGROUPS; group++) {
    heating_rates_by_group[group] = 0.;
  }

  /* "copy" ionization energies from cross section parameters */
  struct rt_photoion_cs_parameters cs_params_cgs =
      rt_init_photoion_cs_params_cgs();
  const double *E_ion_cgs = cs_params_cgs.E_ion;

  /* Integrate energy spectra and cross sections assuming blackbody spectra
   * to obtain estimate for effective cross sections, then use the actual
   * energies present to get the rates */
  /* TODO: check whether we shouldn't be using actual speed of light here */
  const double c_cgs = rt_params.reduced_speed_of_light *
                       units_cgs_conversion_factor(us, UNIT_CONV_VELOCITY);
  const double to_erg = units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);

  /* First, get species number densities and number densities
   * in units of neutral hydrogen number density. */
  double m_p = phys_const->const_proton_mass;
  double species_number_densities_cgs[RT_NIONIZING_SPECIES]; /* in cm^-3 */
  double species_number_densities_nHI[RT_NIONIZING_SPECIES]; /* in nHI^-1 */
  const double to_inv_volume =
      units_cgs_conversion_factor(us, UNIT_CONV_INV_VOLUME);
  const double mass_to_number_density_cgs = to_inv_volume / m_p;
  /* neutral hydrogen */
  species_number_densities_cgs[0] =
      species_densities[0] * mass_to_number_density_cgs;
  species_number_densities_nHI[0] = 1.;
  /* neutral helium */
  species_number_densities_cgs[1] =
      0.25 * species_densities[2] * mass_to_number_density_cgs;
  species_number_densities_nHI[1] =
      0.25 * species_densities[2] / species_densities[0];
  /* singly ionized helium */
  species_number_densities_cgs[2] =
      0.25 * species_densities[3] * mass_to_number_density_cgs;
  species_number_densities_nHI[2] =
      0.25 * species_densities[3] / species_densities[0];

  const double inv_time_cgs =
      units_cgs_conversion_factor(us, UNIT_CONV_INV_TIME);

  /* For the grackle photoionization, we need to
   * keep track of the rates for each species.
   * For the heating rate, we need to sum up all species.
   * To remove the correct amount of energy from the
   * radiation fields, we additionally need to keep track
   * of rates from each photon group. */

  /* store photoionization rate for each species here */
  double ionization_rates_by_species[RT_NIONIZING_SPECIES];
  for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++)
    ionization_rates_by_species[spec] = 0.;

  for (int group = 0; group < RT_NGROUPS; group++) {

    /* Sum results for this group over all species */
    double heating_rate_group_nHI = 0.;
    double heating_rate_group_cgs = 0.;
    float energy_density_i_cgs =
        p->rt_data.radiation[group].energy_density * to_erg * to_inv_volume;

    for (int spec = 0; spec < RT_NIONIZING_SPECIES; spec++) {
      /* Note: the cross sections are in cgs. */
      const double cse = rt_props->energy_weighted_cross_sections[group][spec];
      const double csn = rt_props->number_weighted_cross_sections[group][spec];

      heating_rate_group_nHI +=
          (cse - E_ion_cgs[spec] * csn) * species_number_densities_nHI[spec];
      heating_rate_group_cgs +=
          (cse - E_ion_cgs[spec] * csn) * species_number_densities_cgs[spec];
      ionization_rates_by_species[spec] +=
          energy_density_i_cgs * cse * species_number_densities_cgs[spec] *
          c_cgs / inv_time_cgs; /* internal units T^-1 */
    }

    /* Store total heating rate for grackle */
    rates[0] += heating_rate_group_nHI * c_cgs * energy_density_i_cgs;
    /* Store rates for each group in internal units WITHOUT THE ENERGY DENSITY
     * TERM */
    heating_rates_by_group[group] +=
        heating_rate_group_cgs * c_cgs / inv_time_cgs;
  }

  /* We're done. Write the results in correct place */
  rates[1] = ionization_rates_by_species[0];
  rates[2] = ionization_rates_by_species[1];
  rates[3] = ionization_rates_by_species[2];
}

#endif /* SWIFT_RT_GEAR_INTERACION_RATES_H */
