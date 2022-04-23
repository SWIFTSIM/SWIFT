/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
 *                      Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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
#ifndef SWIFT_ADIABATIC_INDEX_H
#define SWIFT_ADIABATIC_INDEX_H

/**
 * @file adiabatic_index.h
 * @brief Defines the adiabatic index (polytropix index) \f$\gamma\f$ of the
 * problem and (fast) mathematical functions involving it.
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/* Local headers. */
#include "cbrt.h"
#include "error.h"
#include "inline.h"

/* First define some constants */
#if defined(HYDRO_GAMMA_5_3)

#define hydro_gamma 1.66666666666666667f
#define hydro_gamma_minus_one 0.66666666666666667f
#define hydro_gamma_plus_one 2.66666666666666667f
#define hydro_one_over_gamma_minus_one 1.5f
#define hydro_gamma_plus_one_over_two_gamma 0.8f
#define hydro_gamma_minus_one_over_two_gamma 0.2f
#define hydro_gamma_minus_one_over_gamma_plus_one 0.25f
#define hydro_two_over_gamma_plus_one 0.75f
#define hydro_two_over_gamma_minus_one 3.f
#define hydro_gamma_minus_one_over_two 0.33333333333333333f
#define hydro_two_gamma_over_gamma_minus_one 5.f
#define hydro_one_over_gamma 0.6f

#elif defined(HYDRO_GAMMA_7_5)

#define hydro_gamma 1.4f
#define hydro_gamma_minus_one 0.4f
#define hydro_gamma_plus_one 2.4f
#define hydro_one_over_gamma_minus_one 2.5f
#define hydro_gamma_plus_one_over_two_gamma 0.857142857f
#define hydro_gamma_minus_one_over_two_gamma 0.142857143f
#define hydro_gamma_minus_one_over_gamma_plus_one 0.166666667f
#define hydro_two_over_gamma_plus_one 0.833333333
#define hydro_two_over_gamma_minus_one 5.f
#define hydro_gamma_minus_one_over_two 0.2f
#define hydro_two_gamma_over_gamma_minus_one 7.f
#define hydro_one_over_gamma 0.714285714f

#elif defined(HYDRO_GAMMA_4_3)

#define hydro_gamma 1.33333333333333333f
#define hydro_gamma_minus_one 0.33333333333333333f
#define hydro_gamma_plus_one 2.33333333333333333f
#define hydro_one_over_gamma_minus_one 3.f
#define hydro_gamma_plus_one_over_two_gamma 0.875f
#define hydro_gamma_minus_one_over_two_gamma 0.125f
#define hydro_gamma_minus_one_over_gamma_plus_one 0.142857143f
#define hydro_two_over_gamma_plus_one 0.857142857f
#define hydro_two_over_gamma_minus_one 6.f
#define hydro_gamma_minus_one_over_two 0.166666666666666666f
#define hydro_two_gamma_over_gamma_minus_one 8.f
#define hydro_one_over_gamma 0.75f

#elif defined(HYDRO_GAMMA_2_1)

#define hydro_gamma 2.f
#define hydro_gamma_minus_one 1.f
#define hydro_gamma_plus_one 3.f
#define hydro_one_over_gamma_minus_one 1.f
#define hydro_gamma_plus_one_over_two_gamma 0.75f
#define hydro_gamma_minus_one_over_two_gamma 0.25f
#define hydro_gamma_minus_one_over_gamma_plus_one 0.33333333333333333f
#define hydro_two_over_gamma_plus_one 0.66666666666666666f
#define hydro_two_over_gamma_minus_one 2.f
#define hydro_gamma_minus_one_over_two 0.5f
#define hydro_two_gamma_over_gamma_minus_one 4.f
#define hydro_one_over_gamma 0.5f

#else

#error "An adiabatic index needs to be chosen in const.h !"

#endif

/**
 * @brief Returns the argument to the power given by the adiabatic index
 *
 * Computes \f$x^\gamma\f$.
 */
__attribute__((always_inline, const)) INLINE static float pow_gamma(float x) {

#if defined(HYDRO_GAMMA_5_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x); /* x^(-1/3) */
  return icbrt * x * x;          /* x^(5/3) */
#else
  const float cbrt = cbrtf(x);                 /* x^(1/3) */
  return cbrt * cbrt * x;                      /* x^(5/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, 1.4f); /* x^(7/5) */

#elif defined(HYDRO_GAMMA_4_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x);               /* x^(-1/3) */
  return icbrt * icbrt * x * x;                /* x^(4/3) */
#else
  return cbrtf(x) * x;                   /* x^(4/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_2_1)

  return x * x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns the argument to the power given by the adiabatic index minus
 * one
 *
 * Computes \f$x^{(\gamma-1)}\f$.
 */
__attribute__((always_inline, const)) INLINE static float pow_gamma_minus_one(
    float x) {

#if defined(HYDRO_GAMMA_5_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x); /* x^(-1/3) */
  return x * icbrt;              /* x^(2/3) */
#else
  const float cbrt = cbrtf(x);                 /* x^(1/3) */
  return cbrt * cbrt;                          /* x^(2/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, 0.4f); /* x^(2/5) */

#elif defined(HYDRO_GAMMA_4_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x);               /* x^(-1/3) */
  return x * icbrt * icbrt;                    /* x^(1/3) */
#else
  return cbrtf(x);                       /* x^(1/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_2_1)

  return x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns one over the argument to the power given by the adiabatic
 * index minus one
 *
 * Computes \f$x^{-(\gamma-1)}\f$.
 */
__attribute__((always_inline, const)) INLINE static float
pow_minus_gamma_minus_one(float x) {

#if defined(HYDRO_GAMMA_5_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x); /* x^(-1/3) */
  return icbrt * icbrt;          /* x^(-2/3) */
#else
  const float cbrt_inv = 1.f / cbrtf(x);       /* x^(-1/3) */
  return cbrt_inv * cbrt_inv;                  /* x^(-2/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, -0.4f); /* x^(-2/5) */

#elif defined(HYDRO_GAMMA_4_3)

#ifdef WITH_ICBRTF
  return icbrtf(x);                            /* x^(-1/3) */
#else
  return 1.f / cbrtf(x);                 /* x^(-1/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_2_1)

  return 1.f / x;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Returns one over the argument to the power given by the adiabatic
 * index
 *
 * Computes \f$x^{-\gamma}\f$.
 *
 * @param x Argument
 * @return One over the argument to the power given by the adiabatic index
 */
__attribute__((always_inline, const)) INLINE static float pow_minus_gamma(
    float x) {

#if defined(HYDRO_GAMMA_5_3)

#ifdef WITH_ICBRTF
  const float icbrt = icbrtf(x);      /* x^(-1/3) */
  const float icbrt2 = icbrt * icbrt; /* x^(-2/3) */
  return icbrt * icbrt2 * icbrt2;     /* x^(-5/3) */
#else
  const float cbrt_inv = 1.f / cbrtf(x);       /* x^(-1/3) */
  const float cbrt_inv2 = cbrt_inv * cbrt_inv; /* x^(-2/3) */
  return cbrt_inv * cbrt_inv2 * cbrt_inv2;     /* x^(-5/3) */
#endif  // WITH_ICBRTF

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, -1.4f); /* x^(-7/5) */

#elif defined(HYDRO_GAMMA_4_3)

#ifdef WITH_ICBRTF
  const float cbrt_inv = icbrtf(x);            /* x^(-1/3) */
#else
  const float cbrt_inv = 1.f / cbrtf(x); /* x^(-1/3) */
#endif  // WITH_ICBRTF
  const float cbrt_inv2 = cbrt_inv * cbrt_inv; /* x^(-2/3) */
  return cbrt_inv2 * cbrt_inv2;                /* x^(-4/3) */

#elif defined(HYDRO_GAMMA_2_1)

  const float inv = 1.f / x;
  return inv * inv;

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power given by two divided by the adiabatic
 * index minus one
 *
 * Computes \f$x^{\frac{2}{\gamma - 1}}\f$.
 *
 * @param x Argument
 * @return Argument to the power two divided by the adiabatic index minus one
 */
__attribute__((always_inline, const)) INLINE static float
pow_two_over_gamma_minus_one(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return x * x * x; /* x^3 */

#elif defined(HYDRO_GAMMA_7_5)

  const float x2 = x * x;
  const float x3 = x2 * x;
  return x2 * x3;

#elif defined(HYDRO_GAMMA_4_3)

  const float x3 = x * x * x; /* x^3 */
  return x3 * x3;             /* x^6 */

#elif defined(HYDRO_GAMMA_2_1)

  return x * x; /* x^2 */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power given by two times the adiabatic
 * index divided by the adiabatic index minus one
 *
 * Computes \f$x^{\frac{2\gamma}{\gamma - 1}}\f$.
 *
 * @param x Argument
 * @return Argument to the power two times the adiabatic index divided by the
 * adiabatic index minus one
 */
__attribute__((always_inline, const)) INLINE static float
pow_two_gamma_over_gamma_minus_one(float x) {

#if defined(HYDRO_GAMMA_5_3)

  const float x2 = x * x;
  const float x3 = x2 * x;
  return x2 * x3;

#elif defined(HYDRO_GAMMA_7_5)

  const float x2 = x * x;
  const float x4 = x2 * x2;
  return x4 * x2 * x;

#elif defined(HYDRO_GAMMA_4_3)

  const float x2 = x * x;
  const float x4 = x2 * x2;
  return x4 * x4; /* x^8 */

#elif defined(HYDRO_GAMMA_2_1)

  const float x2 = x * x;
  return x2 * x2; /* x^4 */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power given by the adiabatic index minus
 * one  divided by two times the adiabatic index
 *
 * Computes \f$x^{\frac{\gamma - 1}{2\gamma}}\f$.
 *
 * @param x Argument
 * @return Argument to the power the adiabatic index minus one divided by two
 * times the adiabatic index
 */
__attribute__((always_inline, const)) INLINE static float
pow_gamma_minus_one_over_two_gamma(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return powf(x, 0.2f); /* x^0.2 */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, hydro_gamma_minus_one_over_two_gamma);

#elif defined(HYDRO_GAMMA_4_3)

  return powf(x, 0.125f); /* x^0.125 */

#elif defined(HYDRO_GAMMA_2_1)

  return powf(x, 0.25f); /* x^0.25 */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the inverse argument to the power given by the adiabatic index
 * plus one divided by two times the adiabatic index
 *
 * Computes \f$x^{-\frac{\gamma + 1}{2\gamma}}\f$.
 *
 * @param x Argument
 * @return Inverse argument to the power the adiabatic index plus one divided by
 * two times the adiabatic index
 */
__attribute__((always_inline, const)) INLINE static float
pow_minus_gamma_plus_one_over_two_gamma(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return powf(x, -0.8f); /* x^-0.8 */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, -hydro_gamma_plus_one_over_two_gamma);

#elif defined(HYDRO_GAMMA_4_3)

  return powf(x, -0.875f); /* x^-0.875 */

#elif defined(HYDRO_GAMMA_2_1)

  return powf(x, -0.75f); /* x^-0.75 */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power one over the adiabatic index
 *
 * Computes \f$x^{\frac{1}{\gamma}}\f$.
 *
 * @param x Argument
 * @return Argument to the power one over the adiabatic index
 */
__attribute__((always_inline, const)) INLINE static float pow_one_over_gamma(
    float x) {

#if defined(HYDRO_GAMMA_5_3)

  return powf(x, hydro_one_over_gamma); /* x^(3/5) */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, hydro_one_over_gamma);

#elif defined(HYDRO_GAMMA_4_3)

  return powf(x, hydro_one_over_gamma); /* x^(3/4) */

#elif defined(HYDRO_GAMMA_2_1)

  return sqrtf(x); /* x^(1/2) */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power three adiabatic index minus two.
 *
 * Computes \f$x^{3\gamma - 2}\f$.
 *
 * @param x Argument
 */
__attribute__((always_inline, const)) INLINE static float
pow_three_gamma_minus_two(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return x * x * x; /* x^(3) */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, 2.2f); /* x^(11/5) */

#elif defined(HYDRO_GAMMA_4_3)

  return x * x; /* x^(2) */

#elif defined(HYDRO_GAMMA_2_1)

  return x * x * x * x; /* x^(4) */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power three adiabatic index minus five over
 * two.
 *
 * Computes \f$x^{(3\gamma - 5)/2}\f$.
 *
 * @param x Argument
 */
__attribute__((always_inline, const)) INLINE static float
pow_three_gamma_minus_five_over_two(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return 1.f; /* x^(0) */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, -0.4f); /* x^(-2/5) */

#elif defined(HYDRO_GAMMA_4_3)

  return 1.f / sqrtf(x); /* x^(-1/2) */

#elif defined(HYDRO_GAMMA_2_1)

  return sqrtf(x); /* x^(1/2) */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

/**
 * @brief Return the argument to the power three (adiabatic index - 1).
 *
 * Computes \f$x^{3(\gamma - 1)}\f$.
 *
 * @param x Argument
 */
__attribute__((always_inline, const)) INLINE static float
pow_three_gamma_minus_one(float x) {

#if defined(HYDRO_GAMMA_5_3)

  return x * x; /* x^(2) */

#elif defined(HYDRO_GAMMA_7_5)

  return powf(x, 1.2f); /* x^(6/5) */

#elif defined(HYDRO_GAMMA_4_3)

  return x; /* x^(1) */

#elif defined(HYDRO_GAMMA_2_1)

  return x * x * x; /* x^(3) */

#else

  error("The adiabatic index is not defined !");
  return 0.f;

#endif
}

#endif /* SWIFT_ADIABATIC_INDEX_H */
