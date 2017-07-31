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

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "approx_math.h"
#include "const.h"
#include "inline.h"

/* Standard headers */
#include <math.h>

/**
 * @brief Computes the long-range correction term for the FFT calculation.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_eval(
    float u, float *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));

  const float arg1 = u * 0.5f;
  const float arg2 = u * one_over_sqrt_pi;
  const float arg3 = -arg1 * arg1;

  const float term1 = erfcf(arg1);
  const float term2 = arg2 * expf(arg3);

  *W = term1 + term2;
#else

  const float arg = 2.f * u;
  const float exp_arg = good_approx_expf(arg);
  const float term = 1.f / (1.f + exp_arg);

  *W = arg * exp_arg * term * term - exp_arg * term + 1.f;
  *W *= 2.f;
#endif
}

/**
 * @brief Returns the long-range truncation of the Poisson potential in Fourier
 * space.
 *
 * @param u2 The square of the Fourier mode times the cell scale
 * \f$u^2 = k^2r_s^2\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void fourier_kernel_long_grav_eval(
    double u2, double *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION
  *W = exp(-u2);
#else
  const double u = sqrt(u2);
  const double arg = M_PI_2 * u;
  *W = arg / sinh(arg);
#endif
}

#endif  // SWIFT_KERNEL_LONG_GRAVITY_H
