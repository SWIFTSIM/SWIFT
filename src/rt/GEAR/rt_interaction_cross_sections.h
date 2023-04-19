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
#ifndef SWIFT_RT_GEAR_INTERACION_CROSS_SECTIONS_H
#define SWIFT_RT_GEAR_INTERACION_CROSS_SECTIONS_H

#include "inline.h"
#include "rt_species.h"

#include <math.h>
#include <stdlib.h>

/**
 * @file src/rt/GEAR/rt_interaction_cross_sections.h
 * @brief header file concerning photoionization cross sections
 **/

#define RT_INTEGRAL_NPOINTS 10000

/*! Struct containing the parametrized cross section parameters
 * for each photoionizing species.
 * Correct usage is to call rt_init_photoion_cs_params_cgs(),
 * which returns a fully initialized struct. */
struct rt_photoion_cs_parameters {
  double E_ion[rt_ionizing_species_count];      /* erg */
  double E_zero[rt_ionizing_species_count];     /* erg */
  double sigma_zero[rt_ionizing_species_count]; /* cm^-2 */
  double P[rt_ionizing_species_count];
  double ya[rt_ionizing_species_count];
  double yw[rt_ionizing_species_count];
  double y0[rt_ionizing_species_count];
  double y1[rt_ionizing_species_count];
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
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R) (Table E1)
 */
__attribute__((always_inline)) INLINE static struct rt_photoion_cs_parameters
rt_init_photoion_cs_params_cgs(void) {

  struct rt_photoion_cs_parameters photoion_cs_params_cgs;
  double E_ion[rt_ionizing_species_count];
  rt_species_get_ionizing_energy(E_ion);

  photoion_cs_params_cgs.E_ion[rt_ionizing_species_HI] =
      E_ion[rt_ionizing_species_HI];
  photoion_cs_params_cgs.E_zero[rt_ionizing_species_HI] = 6.886e-13;
  photoion_cs_params_cgs.sigma_zero[rt_ionizing_species_HI] = 5.475e-14;
  photoion_cs_params_cgs.P[rt_ionizing_species_HI] = 2.963;
  photoion_cs_params_cgs.ya[rt_ionizing_species_HI] = 32.88;
  photoion_cs_params_cgs.yw[rt_ionizing_species_HI] = 0.;
  photoion_cs_params_cgs.y0[rt_ionizing_species_HI] = 0.;
  photoion_cs_params_cgs.y1[rt_ionizing_species_HI] = 0.;

  photoion_cs_params_cgs.E_ion[rt_ionizing_species_HeI] =
      E_ion[rt_ionizing_species_HeI];
  /* Note: The value of 0.1361 eV for E_0 of HeI given in table 5.1
   * in Rosdahl et al. is an error. The correct value is 13.61 eV */
  photoion_cs_params_cgs.E_zero[rt_ionizing_species_HeI] = 2.181e-11;
  photoion_cs_params_cgs.sigma_zero[rt_ionizing_species_HeI] = 9.492e-16;
  photoion_cs_params_cgs.P[rt_ionizing_species_HeI] = 3.188;
  photoion_cs_params_cgs.ya[rt_ionizing_species_HeI] = 1.469;
  photoion_cs_params_cgs.yw[rt_ionizing_species_HeI] = 2.039;
  photoion_cs_params_cgs.y0[rt_ionizing_species_HeI] = 0.4434;
  photoion_cs_params_cgs.y1[rt_ionizing_species_HeI] = 2.136;

  photoion_cs_params_cgs.E_ion[rt_ionizing_species_HeII] =
      E_ion[rt_ionizing_species_HeII];
  photoion_cs_params_cgs.E_zero[rt_ionizing_species_HeII] = 2.756e-12;
  photoion_cs_params_cgs.sigma_zero[rt_ionizing_species_HeII] = 1.369e-14;
  photoion_cs_params_cgs.P[rt_ionizing_species_HeII] = 2.963;
  photoion_cs_params_cgs.ya[rt_ionizing_species_HeII] = 32.88;
  photoion_cs_params_cgs.yw[rt_ionizing_species_HeII] = 0.;
  photoion_cs_params_cgs.y0[rt_ionizing_species_HeII] = 0.;
  photoion_cs_params_cgs.y1[rt_ionizing_species_HeII] = 0.;

  return photoion_cs_params_cgs;
}

/**
 * Compute the parametrized cross section for a given energy and species.
 * The parametrization is taken from Verner et al. 1996
 * (ui.adsabs.harvard.edu/abs/1996ApJ...465..487V) via Rosdahl et al. 2013
 * (ui.adsabs.harvard.edu/abs/2013MNRAS.436.2188R)
 *
 * @param E energy for which to compute the cross section, in erg
 * @param species index of species, 0 < species < rt_ionizing_species_count
 * @param params cross section parameters struct
 */
__attribute__((always_inline)) INLINE static double
photoionization_cross_section(const double E, const int species,
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

#endif /* SWIFT_RT_GEAR_INTERACION_CROSS_SECTIONS_H */
