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

#define RT_INTEGRAL_NPOINTS 10000

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
