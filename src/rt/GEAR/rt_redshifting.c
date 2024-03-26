/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Stan Verhoeve (s06verhoeve@gmail.com)
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

#include "rt_redshifting.h"
#include "rt_blackbody.h"
#include "rt_properties.h"
#include "rt_interaction_cross_sections.h"
#include <gsl/gsl_integration.h>
#include "error.h"


static double rt_redshift_get_spectrum_derivative(const double nu, void *params) {
  rt_spectrum_integration_params *pars = (struct rt_spectrum_integration_params *)params;

  if (pars->spectrum_type ==0) {
    /* Constant spectrum */
    return 0.;
  } else if (pars->spectrum_type == 1) {
    /* Blackbody spectrum */
    const double T = pars->T;
    const double kB = pars->kB;
    const double h_planck - pars->h_planck;
    const double c = pars->c;
    return blackbody_spectrum_intensity_first_derivative(nu, T, kB, h_planck, c);
  } else {
    error("Unknown stellar spectrum type selected: %d", pars->spectrum_type);
    return 0.;
  }
}


static double spectrum_derivative_times_nu_integrand(double nu, void *params) {
  struct rt_spectrum_integration_params *p = 
      (struct rt_spectrum_integration_params *)params;
  const double J = rt_redshift_get_spectrum_derivative(nu, params);

  return nu * J;
}


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


/* Integrate the spectrum */
static double rt_redshift_integrate_gsl( double (*function)(double, void *), double nu_start, double nu_stop,
    int npoints, struct rt_spectrum_integration_params *params) {
  
  gsl_function F;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npoints);
  double result, error;

  F.function = function;
  F.params = (void *)params;

  gsl_integration_qags(&F, nu_start, nu_stop, /*espabs=*/0., /*epsrel=*/1e-7,
		       npoints, w, &result, &error);

  /* Clean up after yourself. */
  gsl_integration_workspace_free(w);

  return result;
}


void rt_redshift_init(struct rt_props *restrict rt_props,
		      const struct phys_const *restrict phys_const,
		      const struct unit_system *restrict us) {

  /* Allocate the space to store the redshift integrals */
  double *av_redshift_energy = rt_props->average_redshift_energy;
  double *av_energy = rt_props->average_photon_energy;
  double integral_nuE[RT_NGROUPS];
  
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
  
  // NOTE: Not needed for the integration, but may be needed for the frequency bounds
  struct rt_photoion_cs_parameters cs_params = rt_init_photoion_cs_params_cgs(); 

  struct rt_spectrum_integration_params integration_params = {
      /*species=*/0,
      /*spectrum_type=*/spectype,
      /*freq_max for const spectrum=*/maxfreq_const_spectrum,
      /*T=*/T_bb,
      /*kB=*/kB_cgs,
      /*h_planck=*/h_planck_cgs,
      /*c=*/c_cgs};

  
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
    
    /* For redshifting, it does not matter where we start integrating
     * (even if we only have one group). */
    nu_start[0] = cs_params.E_ion[rt_ionizing_species_HI] / h_planck_cgs;
    if (engine_rank == 0)
      message(
          "WARNING: with only 1 photon group, I'll start integrating"
          " the redshift energy at the first ionizing frequency %.3g",
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


  /* Compute integrals */
  /* ----------------- */

  for (int g=0; g < RT_NGROUPS; g++) {
    integral_nuE[g] = rt_redshift_integrate_gsl(spectrum_derivative_times_nu_integrand, 
		      nu_start[g], nu_stop[g], RT_INTEGRAL_NPOINTS, 
		      &integration_params);
    
  }
  /* Now compute the actual average redshifted energy */
  /* ------------------------------------------------ */
  for (int g = 0; g < RT_NGROUPS; g++) {
    if (avg_energy[g] > 0.) {
      av_redshift_energy[g] = integral_nuE[g] / avg_energy[g];
    } else {
      av_redshift_energy[g] = 0.;
    }
  }
}
