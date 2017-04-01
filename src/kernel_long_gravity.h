/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_KERNEL_LONG_GRAVITY_H
#define SWIFT_KERNEL_LONG_GRAVITY_H

#include <math.h>

/* Includes. */
#include "const.h"
#include "inline.h"
#include "vector.h"

#define one_over_sqrt_pi ((float)(M_2_SQRTPI * 0.5))

/**
 * @brief Computes the long-range correction term for the FFT calculation.
 *
 * @param u The ratio of the distance to the FFT cell scale $u = x/A$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_eval(
    float u, float *const W) {

  /* const float arg1 = u * 0.5f; */
  /* const float arg2 = u * one_over_sqrt_pi; */
  /* const float arg3 = -arg1 * arg1; */

  /* const float term1 = erfcf(arg1); */
  /* const float term2 = arg2 * expf(arg3); */

  /* *W = term1 + term2; */
  *W = 1.f;
}

#endif  // SWIFT_KERNEL_LONG_GRAVITY_H
