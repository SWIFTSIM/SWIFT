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

  /* Powers of u=r/2r_s */
  const float u = 0.5f * r * r_s_inv;
  const float u2 = u * u;
  const float u3 = u2 * u;
  const float u4 = u3 * u;

  /* Powers of (1/r_s) */
  const float r_s_inv2 = r_s_inv * r_s_inv;
  const float r_s_inv3 = r_s_inv2 * r_s_inv;
  const float r_s_inv4 = r_s_inv3 * r_s_inv;
  const float r_s_inv5 = r_s_inv4 * r_s_inv;

  /* Derivatives of \chi */
#ifdef GRAVITY_USE_EXACT_LONG_RANGE_MATH
  derivs->chi_0 = erfcf(u);
#else
  derivs->chi_0 = approx_erfcf(u);
#endif
  derivs->chi_1 = -r_s_inv;
  derivs->chi_2 = r_s_inv2 * u;
  derivs->chi_3 = -r_s_inv3 * (u2 - 0.5f);
  derivs->chi_4 = r_s_inv4 * (u3 - 1.5f * u);
  derivs->chi_5 = -r_s_inv5 * (u4 - 3.f * u2 + 0.75f);

  const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));
  const float common_factor = one_over_sqrt_pi * expf(-u2);

  /* Multiply in the common factors */
  derivs->chi_1 *= common_factor;
  derivs->chi_2 *= common_factor;
  derivs->chi_3 *= common_factor;
  derivs->chi_4 *= common_factor;
  derivs->chi_5 *= common_factor;

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
 * @brief Computes the long-range correction term for the potential calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 */
__attribute__((const)) INLINE static float kernel_long_grav_pot_eval(
    const float u) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  const float arg1 = u * 0.5f;
#ifdef GRAVITY_USE_EXACT_LONG_RANGE_MATH
  return erfcf(arg1);
#else
  return approx_erfcf(arg1);
#endif

#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2 - 2 exp(x) * alpha */
  float W = 1.f - alpha * exp_x;
  W = W * 2.f;

  return W;
#endif
}

/**
 * @brief Computes the long-range correction term for the force calculation
 * coming from FFT.
 *
 * @param u The ratio of the distance to the FFT cell scale \f$u = r/r_s\f$.
 */
__attribute__((const)) INLINE static float kernel_long_grav_force_eval(
    const float u) {

#ifdef GADGET2_LONG_RANGE_CORRECTION

  static const float one_over_sqrt_pi = ((float)(M_2_SQRTPI * 0.5));

  const float arg1 = u * 0.5f;
  const float arg2 = -arg1 * arg1;

#ifdef GRAVITY_USE_EXACT_LONG_RANGE_MATH
  const float term1 = erfcf(arg1);
#else
  const float term1 = approx_erfcf(arg1);
#endif
  const float term2 = u * one_over_sqrt_pi * expf(arg2);

  return term1 + term2;
#else

  const float x = 2.f * u;
  const float exp_x = expf(x);  // good_approx_expf(x);
  const float alpha = 1.f / (1.f + exp_x);

  /* We want 2*(x*alpha - x*alpha^2 - exp(x)*alpha + 1) */
  float W = 1.f - alpha;
  W = W * x - exp_x;
  W = W * alpha + 1.f;
  W = W * 2.f;

  return W;
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
