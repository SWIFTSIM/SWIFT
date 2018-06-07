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
#include <float.h>
#include <math.h>

//#define GADGET2_LONG_RANGE_CORRECTION

struct truncated_derivatives {

  float chi_0;
  float chi_1;
  float chi_2;
  float chi_3;
  float chi_4;
  float chi_5;
};

__attribute__((always_inline)) INLINE static void kernel_long_grav_derivatives(
    const float r, const float rs_inv,
    struct truncated_derivatives *const derivs) {

  const float constant1 = 2.f * rs_inv;
  const float x = constant1 * r;

  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha_inv = 1.f + exp_x;

  const float alpha1 = 1.f / alpha_inv;
  const float alpha2 = alpha1 * alpha1;
  const float alpha3 = alpha2 * alpha1;
  const float alpha4 = alpha3 * alpha1;
  const float alpha5 = alpha4 * alpha1;
  const float alpha6 = alpha5 * alpha1;

  const float constant2 = constant1 * constant1;
  const float constant3 = constant2 * constant1;
  const float constant4 = constant3 * constant1;
  const float constant5 = constant4 * constant1;

  derivs->chi_0 = 2.f * (1.f - exp_x * alpha1);
  derivs->chi_1 = constant1 * (2.f * alpha2 - 2.f * alpha1);
  derivs->chi_2 = constant2 * (4.f * alpha3 - 6.f * alpha2 + 2.f * alpha1);
  derivs->chi_3 = constant3 * (12.f * alpha4 - 24.f * alpha3 + 14.f * alpha2 -
                               2.f * alpha1);
  derivs->chi_4 = constant4 * (48.f * alpha5 - 120.f * alpha4 + 100.f * alpha3 -
                               30.f * alpha2 + 2.f * alpha1);
  derivs->chi_5 =
      constant5 * (240.f * alpha6 - 720.f * alpha5 + 780.f * alpha4 -
                   360.f * alpha3 + 62.f * alpha2 - 2.f * alpha1);
}

/**
 * @brief Computes the long-range correction term for the potential calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_pot_eval(
    const float u, float *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float arg1 = u * 0.5f;
  const float term1 = erfcf(arg1);

  *W = term1;
#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2 - 2 exp(x) * alpha */
  *W = 1.f - alpha * exp_x;
  *W *= 2.f;
#endif
}

/**
 * @brief Computes the long-range correction term for the force calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline)) INLINE static void kernel_long_grav_force_eval(
    const float u, float *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));

  const float arg1 = u * 0.5f;
  const float arg2 = u * one_over_sqrt_pi;
  const float arg3 = -arg1 * arg1;

  const float term1 = erfcf(arg1);
  const float term2 = arg2 * expf(arg3);

  *W = term1 + term2;
#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2*(x*alpha - x*alpha^2 - exp(x)*alpha + 1) */
  *W = 1.f - alpha;
  *W = *W * x - exp_x;
  *W = *W * alpha + 1.f;
  *W *= 2.f;

/* const float arg = 2.f * u; */
/* const float exp_arg = good_approx_expf(arg); */
/* const float term = 1.f / (1.f + exp_arg); */

/* *W = arg * exp_arg * term * term - exp_arg * term + 1.f; */
/* *W *= 2.f; */
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
    const double u2, double *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION
  *W = exp(-u2);
#else
  const double u = sqrt(u2);
  const double arg = M_PI_2 * u;
  *W = arg / (sinh(arg) + FLT_MIN);
#endif
}

#endif /* SWIFT_KERNEL_LONG_GRAVITY_H */
