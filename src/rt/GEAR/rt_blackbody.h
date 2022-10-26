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
#ifndef SWIFT_RT_GEAR_BLACKBODY_H
#define SWIFT_RT_GEAR_BLACKBODY_H

/**
 * @file src/rt/GEAR/rt_blackbody.h
 * @brief Functions related to the blackbody spectrum
 */

#include "inline.h"

#include <math.h>

/**
 * @brief Return the specific intensity of the blackbody spectrum
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
__attribute__((always_inline)) INLINE double blackbody_spectrum_intensity(
    const double nu, const double T, const double kB, const double h_planck,
    const double c) {

  const double hnu = h_planck * nu;
  const double kT = kB * T;
  const double nu2 = nu * nu;
  double temp;
  if (hnu / kT < 1e-6) {
    /* prevent division by zero, use Taylor approximation */
    temp = kT;
  } else if (hnu / kT > 700.) {
    /* prevent infs */
    temp = 0.;
  } else {
    temp = 1. / (exp(hnu / kT) - 1.);
  }
  return 2. * hnu * nu2 / (c * c) * temp;
}

/**
 * Return the blackbody spectrum energy density
 *
 * @param nu frequency at which to compute specific intensity
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 * @param c speed of light
 */
__attribute__((always_inline)) INLINE double blackbody_spectrum_energy_density(
    const double nu, const double T, const double kB, const double h_planck,
    const double c) {
  return 4. * M_PI / c * blackbody_spectrum_intensity(nu, T, kB, h_planck, c);
}

/**
 * Return the frequency at which the blackbody spectrum at a given tempterature
 * peaks.
 *
 * @param T temperature characterizing the spectrum
 * @param kB Boltzmann constant
 * @param h_planck Planck's constant
 */
__attribute__((always_inline)) INLINE double blackbody_peak_frequency(
    const double T, const double kB, const double h_planck) {
  return 2.82144 * kB * T / h_planck;
}

#endif
