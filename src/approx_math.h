/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_APPROX_MATH_H
#define SWIFT_APPROX_MATH_H

#include "inline.h"

/**
 * @brief Approximate version of the complementay error function erfcf(x).
 *
 * This is based on eq. 7.1.27 of Abramowitz & Stegun, 1972.
 * The absolute error is < 4.7*10^-4 over the range 0 < x < infinity.
 *
 * Returns garbage for x < 0.
 * @param x The number to compute erfc for.
 */
__attribute__((always_inline, const)) INLINE static float approx_erfcf(
    float x) {

  /* 1 + 0.278393*x + 0.230389*x^2 + 0.000972*x^3 + 0.078108*x^4 */
  float arg = 0.078108f;
  arg = x * arg + 0.000972f;
  arg = x * arg + 0.230389f;
  arg = x * arg + 0.278393f;
  arg = x * arg + 1.f;

  /* 1 / arg^4 */
  const float arg2 = arg * arg;
  const float arg4 = arg2 * arg2;
  return 1.f / arg4;
}

/**
 * @brief Approximate version of expf(x) using a 4th order Taylor expansion
 *
 * The absolute error is smaller than 3 * 10^-6 for -0.2 < x < 0.2.
 * The absolute error is smaller than 3 * 10^-7 for -0.1 < x < 0.1.

 * The relative error is smaller than 1 * 10^-6 for -0.2 < x < 0.2.
 * The relative error is smaller than 3 * 10^-7 for -0.1 < x < 0.1.
 *
 * @param x The number to take the exponential of.
 */
__attribute__((always_inline, const)) INLINE static float approx_expf(float x) {
  return 1.f + x * (1.f + x * (0.5f + x * (1.f / 6.f + 1.f / 24.f * x)));
}

/**
 * @brief Approximate version of expf(x) using a 6th order Taylor expansion
 *
 */
__attribute__((always_inline, const)) INLINE static float good_approx_expf(
    float x) {
  return 1.f +
         x * (1.f +
              x * (0.5f + x * ((1.f / 6.f) +
                               x * ((1.f / 24.f) +
                                    x * ((1.f / 120.f) + (1.f / 720.f) * x)))));
}

/**
 * @brief Approximate version of exp(x) using a 6th order Taylor expansion
 */
__attribute__((always_inline, const)) INLINE static double good_approx_exp(
    double x) {
  return 1. +
         x * (1. + x * (0.5 + x * ((1. / 6.) +
                                   x * ((1. / 24.) +
                                        x * ((1. / 120.) + (1. / 720.) * x)))));
}

#endif /* SWIFT_APPROX_MATH_H */
