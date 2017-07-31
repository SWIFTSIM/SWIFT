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
 * @brief Approximate version of expf(x) using a 4th order Taylor expansion
 *
 * The absolute error is smaller than 3 * 10^-6 for -0.2 < x < 0.2.
 * The absolute error is smaller than 2 * 10^-7 for -0.1 < x < 0.1.

 * The relative error is smaller than 1 * 10^-6 for -0.2 < x < 0.2.
 * The relative error is smaller than 4 * 10^-8 for -0.1 < x < 0.1.
 *
 * @param x The number to take the exponential of.
 */
__attribute__((always_inline)) INLINE static float approx_expf(float x) {
  return 1.f + x * (1.f + x * (0.5f + x * (1.f / 6.f + 1.f / 24.f * x)));
}

/**
 * @brief Approximate version of expf(x) using a 6th order Taylor expansion
 *
 */
__attribute__((always_inline)) INLINE static float good_approx_expf(float x) {
  return 1.f +
         x * (1.f +
              x * (0.5f +
                   x * ((1.f / 6.f) +
                        x * ((1.f / 24.f) +
                             x * ((1.f / 120.f) + (1.f / 720.f) * x)))));
}

#endif /* SWIFT_APPROX_MATH_H */
