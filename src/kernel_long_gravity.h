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
#include "const.h"
#include "exp.h"
#include "inline.h"

/* Standard headers */
#include <float.h>
#include <math.h>

#define GADGET2_LONG_RANGE_CORRECTION

#ifdef GADGET2_LONG_RANGE_CORRECTION
#define kernel_long_gravity_truncation_name "Gadget-like (using erfc())"
#else
#define kernel_long_gravity_truncation_name "Exp-based Sigmoid"
#endif

/**
 * @brief Derivatives of the long-range truncation function \f$\chi(r,r_s)\f$ up
 * to 5th order.
 */
struct chi_derivatives {

  /*! 0th order derivative \f$\chi(r,r_s)\f$ */
  float chi_0;

  /*! 1st order derivative \f$\partial_{r}\chi(r,r_s)\f$ */
  float chi_1;

  /*! 2nd order derivative \f$\partial_{rr}\chi(r,r_s)\f$ */
  float chi_2;

  /*! 3rd order derivative \f$\partial_{rrr}\chi(r,r_s)\f$ */
  float chi_3;

  /*! 4th order derivative \f$\partial_{rrrr}\chi(r,r_s)\f$ */
  float chi_4;

  /*! 5th order derivative \f$\partial_{rrrrr}\chi(r,r_s)\f$ */
  float chi_5;
};

/**
 * @brief Compute the derivatives of the long-range truncation function
 * \f$\chi(r,r_s)\f$ up to 5th order.
 *
 * @param r The distance.
 * @param r_s_inv The inverse of the long-range gravity mesh scale.
 * @param derivs (return) The computed #chi_derivatives.
 */
__attribute__((always_inline, nonnull)) INLINE static void
kernel_long_grav_derivatives(const float r, const float r_s_inv,
                             struct chi_derivatives *const derivs) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  /* Powers of u = (1/2) * (r / r_s) */
  const float u = 0.5f * r * r_s_inv;
  const float u2 = u * u;
  const float u4 = u2 * u2;

  const float exp_u2 = optimized_expf(-u2);

  /* Compute erfcf(u) using eq. 7.1.25 of
   * Abramowitz & Stegun, 1972.
   *
   * This has a *relative* error of less than 4e-3 over
   * the range of interest (0 < u <  5) */

  const float t = 1.f / (1.f + 0.47047f * u);

  /* 0.3480242 * t - 0.0958798 * t^2 + 0.7478556 * t^3 */
  float a = 0.7478556f;
  a = a * t - 0.0958798f;
  a = a * t + 0.3480242f;
  a = a * t;

  const float erfc_u = a * exp_u2;

  /* C = (1/sqrt(pi)) * expf(-u^2) */
  const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));
  const float common_factor = one_over_sqrt_pi * exp_u2;

  /* (1/r_s)^n * C */
  const float r_s_inv_times_C = r_s_inv * common_factor;
  const float r_s_inv2_times_C = r_s_inv_times_C * r_s_inv;
  const float r_s_inv3_times_C = r_s_inv2_times_C * r_s_inv;
  const float r_s_inv4_times_C = r_s_inv3_times_C * r_s_inv;
  const float r_s_inv5_times_C = r_s_inv4_times_C * r_s_inv;

  /* Now, compute the derivatives of \chi */
#ifdef GRAVITY_USE_EXACT_LONG_RANGE_MATH

  /* erfc(u) */
  derivs->chi_0 = erfcf(u);
#else

  /* erfc(u) */
  derivs->chi_0 = erfc_u;
#endif

  /* (-1/r_s) * (1/sqrt(pi)) * expf(-u^2) */
  derivs->chi_1 = -r_s_inv_times_C;

  /* (1/r_s)^2 * u * (1/sqrt(pi)) * expf(-u^2) */
  derivs->chi_2 = r_s_inv2_times_C * u;

  /* (1/r_s)^3 * (1/2 - u^2) * (1/sqrt(pi)) * expf(-u^2) */
  derivs->chi_3 = r_s_inv3_times_C * (0.5f - u2);

  /* (1/r_s)^4 * (u^3 - 3/2 u) * (1/sqrt(pi)) * expf(-u^2) */
  derivs->chi_4 = r_s_inv4_times_C * (u2 - 1.5f) * u;

  /* (1/r_s)^5 * (3/4 - 3u^2 + u^4) * (1/sqrt(pi)) * expf(-u^2) */
  derivs->chi_5 = r_s_inv5_times_C * (0.75f - 3.f * u2 + u4);

#else

  /* Powers of 2/r_s */
  const float c0 = 1.f;
  const float c1 = 2.f * r_s_inv;
  const float c2 = c1 * c1;
  const float c3 = c2 * c1;
  const float c4 = c3 * c1;
  const float c5 = c4 * c1;

  /* 2r / r_s */
  const float x = c1 * r;

  /* e^(2r / r_s) */
  const float exp_x = expf(x);  // good_approx_expf(x);

  /* 1 / alpha(w) */
  const float a_inv = 1.f + exp_x;

  /* Powers of alpha */
  const float a1 = 1.f / a_inv;
  const float a2 = a1 * a1;
  const float a3 = a2 * a1;
  const float a4 = a3 * a1;
  const float a5 = a4 * a1;
  const float a6 = a5 * a1;

  /* Derivatives of \chi */
  derivs->chi_0 = -2.f * exp_x * c0 * a1 + 2.f;
  derivs->chi_1 = -2.f * exp_x * c1 * a2;
  derivs->chi_2 = -2.f * exp_x * c2 * (2.f * a3 - a2);
  derivs->chi_3 = -2.f * exp_x * c3 * (6.f * a4 - 6.f * a3 + a2);
  derivs->chi_4 = -2.f * exp_x * c4 * (24.f * a5 - 36.f * a4 + 14.f * a3 - a2);
  derivs->chi_5 = -2.f * exp_x * c5 *
                  (120.f * a6 - 240.f * a5 + 150.f * a4 - 30.f * a3 + a2);
#endif
}

/**
 * @brief Computes the long-range correction terms for the potential and
 * force calculations due to the mesh truncation.
 *
 * We use an approximation to the erfc() that gives a *relative* accuracy
 * for the potential tem of 3.4e-3 and 2.4e-4 for the force term over the
 * range [0, 5] of r_over_r_s.
 * The accuracy is much better in the range [0, 2] (6e-5 and 2e-5 respectively).
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 */
__attribute__((nonnull)) INLINE static void kernel_long_grav_eval(
    const float r_over_r_s, float *restrict corr_f, float *restrict corr_pot) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float two_over_sqrt_pi = ((float)M_2_SQRTPI);

  const float u = 0.5f * r_over_r_s;
  const float u2 = u * u;
  const float exp_u2 = optimized_expf(-u2);

  /* Compute erfcf(u) using eq. 7.1.25 of
   * Abramowitz & Stegun, 1972.
   *
   * This has a *relative* error of less than 4e-3 over
   * the range of interest (0 < u <  5) */

  const float t = 1.f / (1.f + 0.47047f * u);

  /* 0.3480242 * t - 0.0958798 * t^2 + 0.7478556 * t^3 */
  float a = 0.7478556f;
  a = a * t - 0.0958798f;
  a = a * t + 0.3480242f;
  a = a * t;

  const float erfc_u = a * exp_u2;

  *corr_pot = erfc_u;
  *corr_f = erfc_u + two_over_sqrt_pi * u * exp_u2;

#else
  const float x = 2.f * r_over_r_s;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2 - 2 exp(x) * alpha */
  float W = 1.f - alpha * exp_x;
  W = W * 2.f;

  *corr_pot = W;

  /* We want 2*(x*alpha - x*alpha^2 - exp(x)*alpha + 1) */
  W = 1.f - alpha;
  W = W * x - exp_x;
  W = W * alpha + 1.f;
  W = W * 2.f;

  *corr_f = W;
#endif
}

/**
 * @brief Computes the long-range correction term for the force calculation
 * coming from FFT in double precision.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 * @param W (return) The value of the kernel function.
 */
__attribute__((always_inline, nonnull)) INLINE static void
kernel_long_grav_force_eval_double(const double u, double *const W) {
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
#ifdef GADGET2_LONG_RANGE_CORRECTION

  const double one_over_sqrt_pi = M_2_SQRTPI * 0.5;

  const double arg1 = u * 0.5;
  const double arg2 = -arg1 * arg1;

  const double term1 = erfc(arg1);
  const double term2 = u * one_over_sqrt_pi * exp(arg2);

  *W = term1 + term2;
#else

  const double x = 2. * u;
  const double exp_x = exp(x);
  const double alpha = 1. / (1. + exp_x);

  /* We want 2*(x*alpha - x*alpha^2 - exp(x)*alpha + 1) */
  *W = 1. - alpha;
  *W = *W * x - exp_x;
  *W = *W * alpha + 1.;
  *W *= 2.;
#endif
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
__attribute__((always_inline, nonnull)) INLINE static void
fourier_kernel_long_grav_eval(const double u2, double *const W) {

#ifdef GADGET2_LONG_RANGE_CORRECTION
  *W = exp(-u2);
#else
  const double u = sqrt(u2);
  const double arg = M_PI_2 * u;
  *W = arg / (sinh(arg) + FLT_MIN);
#endif
}

#endif /* SWIFT_KERNEL_LONG_GRAVITY_H */
