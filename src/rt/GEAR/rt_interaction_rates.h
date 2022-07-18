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

#include "inline.h"

#include <math.h>

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
__attribute__((always_inline)) INLINE static struct rt_photoion_cs_parameters
rt_init_photoion_cs_params_cgs(void) {

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

#endif /* SWIFT_RT_GEAR_INTERACION_RATES_H */
