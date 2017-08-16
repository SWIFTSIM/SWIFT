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
#ifndef SWIFT_GRAVITY_DERIVATIVE_H
#define SWIFT_GRAVITY_DERIVATIVE_H

/**
 * @file gravity_derivatives.h
 * @brief Derivatives (up to 5th order) of the gravitational potential.
 *
 * We use the notation of Dehnen, Computational Astrophysics and Cosmology,
 * 1, 1, pp. 24 (2014), arXiv:1405.2255
 */

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "inline.h"

struct potential_derivatives {

  /* 0th order terms */
  float D_000;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 1st order terms */
  float D_100, D_010, D_001;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 2nd order terms */
  float D_200, D_020, D_002;
  float D_110, D_101, D_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 3rd order terms */
  float D_300, D_030, D_003;
  float D_210, D_201;
  float D_120, D_021;
  float D_102, D_012;
  float D_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 4th order terms */
  float D_400, D_040, D_004;
  float D_310, D_301;
  float D_130, D_031;
  float D_103, D_013;
  float D_220, D_202, D_022;
  float D_211, D_121, D_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  /* 5th order terms */
  float D_005, D_014, D_023;
  float D_032, D_041, D_050;
  float D_104, D_113, D_122;
  float D_131, D_140, D_203;
  float D_212, D_221, D_230;
  float D_302, D_311, D_320;
  float D_401, D_410, D_500;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif
};

__attribute__((always_inline)) INLINE static void compute_potential_derivatives(
    float r_x, float r_y, float r_z, float r_inv,
    struct potential_derivatives *pot) {

/* Compute some (odd) powers of 1/r */
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  const float r_inv2 = r_inv * r_inv;
  const float r_inv3 = -1.f * r_inv * r_inv2; /* -1 / r^3 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  const float r_inv5 = -3.f * r_inv3 * r_inv2; /* 3 / r^5 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  const float r_inv7 = -5.f * r_inv5 * r_inv2; /* -15 / r^7 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  const float r_inv9 = -7.f * r_inv7 * r_inv2; /* 105 / r^9 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  const float r_inv11 = -9.f * r_inv9 * r_inv2; /* -945 / r^11 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

/* Compute some powers of r_x, r_y and r_z */
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  const float r_x2 = r_x * r_x;
  const float r_y2 = r_y * r_y;
  const float r_z2 = r_z * r_z;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  const float r_x3 = r_x2 * r_x;
  const float r_y3 = r_y2 * r_y;
  const float r_z3 = r_z2 * r_z;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  const float r_x4 = r_x3 * r_x;
  const float r_y4 = r_y3 * r_y;
  const float r_z4 = r_z3 * r_z;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  const float r_x5 = r_x4 * r_x;
  const float r_y5 = r_y4 * r_y;
  const float r_z5 = r_z4 * r_z;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

  /* Get the 0th order term */
  pot->D_000 = r_inv;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* 1st order derivatives */
  pot->D_100 = r_x * r_inv3;
  pot->D_010 = r_y * r_inv3;
  pot->D_001 = r_z * r_inv3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* 2nd order derivatives */
  pot->D_200 = r_x2 * r_inv5 + r_inv3;
  pot->D_020 = r_y2 * r_inv5 + r_inv3;
  pot->D_002 = r_z2 * r_inv5 + r_inv3;
  pot->D_110 = r_x * r_y * r_inv5;
  pot->D_101 = r_x * r_z * r_inv5;
  pot->D_011 = r_y * r_z * r_inv5;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* 3rd order derivatives */
  pot->D_300 = r_x3 * r_inv7 + 3.f * r_x * r_inv5;
  pot->D_030 = r_y3 * r_inv7 + 3.f * r_y * r_inv5;
  pot->D_003 = r_z3 * r_inv7 + 3.f * r_z * r_inv5;
  pot->D_210 = r_x2 * r_y * r_inv7 + r_y * r_inv5;
  pot->D_201 = r_x2 * r_z * r_inv7 + r_z * r_inv5;
  pot->D_120 = r_y2 * r_x * r_inv7 + r_x * r_inv5;
  pot->D_021 = r_y2 * r_z * r_inv7 + r_z * r_inv5;
  pot->D_102 = r_z2 * r_x * r_inv7 + r_x * r_inv5;
  pot->D_012 = r_z2 * r_y * r_inv7 + r_y * r_inv5;
  pot->D_111 = r_x * r_y * r_z * r_inv7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* 4th order derivatives */
  pot->D_400 = r_x4 * r_inv9 + 6.f * r_x2 * r_inv7 + 3.f * r_inv5;
  pot->D_040 = r_y4 * r_inv9 + 6.f * r_y2 * r_inv7 + 3.f * r_inv5;
  pot->D_004 = r_z4 * r_inv9 + 6.f * r_z2 * r_inv7 + 3.f * r_inv5;
  pot->D_310 = r_x3 * r_y * r_inv9 + 3.f * r_x * r_y * r_inv7;
  pot->D_301 = r_x3 * r_z * r_inv9 + 3.f * r_x * r_z * r_inv7;
  pot->D_130 = r_y3 * r_x * r_inv9 + 3.f * r_y * r_x * r_inv7;
  pot->D_031 = r_y3 * r_z * r_inv9 + 3.f * r_y * r_z * r_inv7;
  pot->D_103 = r_z3 * r_x * r_inv9 + 3.f * r_z * r_x * r_inv7;
  pot->D_013 = r_z3 * r_y * r_inv9 + 3.f * r_z * r_y * r_inv7;
  pot->D_220 = r_x2 * r_y2 * r_inv9 + r_x2 * r_inv7 + r_y2 * r_inv7 + r_inv5;
  pot->D_202 = r_x2 * r_z2 * r_inv9 + r_x2 * r_inv7 + r_z2 * r_inv7 + r_inv5;
  pot->D_022 = r_y2 * r_z2 * r_inv9 + r_y2 * r_inv7 + r_z2 * r_inv7 + r_inv5;
  pot->D_211 = r_x2 * r_y * r_z * r_inv9 + r_y * r_z * r_inv7;
  pot->D_121 = r_y2 * r_x * r_z * r_inv9 + r_x * r_z * r_inv7;
  pot->D_112 = r_z2 * r_x * r_z * r_inv9 + r_x * r_y * r_inv7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif
}

#if 0

/*************************/
/* 0th order derivatives */
/*************************/

/**
 * @brief \f$ \phi(r_x, r_y, r_z) \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_000(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {

  return r_inv;
}

/*************************/
/* 1st order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_100(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {

  return -r_x * r_inv * r_inv * r_inv;
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_010(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {

  return -r_y * r_inv * r_inv * r_inv;
}

/**
 * @brief \f$ \frac{\partial\phi(r_x, r_y, r_z)}{\partial r_x} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_001(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {

  return -r_z * r_inv * r_inv * r_inv;
}

/*************************/
/* 2nd order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_200(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_x * r_x * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_y^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_020(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_y * r_y * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_z^2} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_002(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv3 = r_inv * r_inv2;
  const double r_inv5 = r_inv3 * r_inv2;
  return 3. * r_z * r_z * r_inv5 - r_inv3;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x\partial r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_110(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_x * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_x\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_101(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_x * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^2\phi(r_x, r_y, r_z)}{\partial r_y\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_011(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  return 3. * r_y * r_z * r_inv5;
}

/*************************/
/* 3rd order derivatives */
/*************************/

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_300(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_x * r_inv7 + 9. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_030(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_y * r_y * r_y * r_inv7 + 9. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_z^3} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_003(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_z * r_z * r_z * r_inv7 + 9. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^2\partial r_y}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_210(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_y * r_inv7 + 3. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x^2\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_201(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_x * r_z * r_inv7 + 3. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x\partial r_y^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_120(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_y * r_y * r_inv7 + 3. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y^2\partial r_z}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_021(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_z * r_y * r_y * r_inv7 + 3. * r_z * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_x\partial r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_102(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_x * r_z * r_z * r_inv7 + 3. * r_x * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_y\partial r_z^2}
 * \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_012(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv2 = r_inv * r_inv;
  const double r_inv5 = r_inv2 * r_inv2 * r_inv;
  const double r_inv7 = r_inv5 * r_inv2;
  return -15. * r_y * r_z * r_z * r_inv7 + 3. * r_y * r_inv5;
}

/**
 * @brief \f$ \frac{\partial^3\phi(r_x, r_y, r_z)}{\partial r_z\partial
 * r_y\partial r_z} \f$.
 *
 * @param r_x x-coordinate of the distance vector (\f$ r_x \f$).
 * @param r_y y-coordinate of the distance vector (\f$ r_y \f$).
 * @param r_z z-coordinate of the distance vector (\f$ r_z \f$).
 * @param r_inv Inverse of the norm of the distance vector (\f$ |r|^{-1} \f$)
 */
__attribute__((always_inline)) INLINE static double D_111(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  const double r_inv3 = r_inv * r_inv * r_inv;
  const double r_inv7 = r_inv3 * r_inv3 * r_inv;
  return -15. * r_x * r_y * r_z * r_inv7;
}

/*********************************/
/* 4th order gravity derivatives */
/*********************************/

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_z^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_004(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^1 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_013(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^2 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_022(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^3 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_031(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_y^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_040(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_y * r_y) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_103(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^1 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_112(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_y);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^2 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_121(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_z);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^1 \partial_y^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_130(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_y);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_202(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_z * r_z) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_y^1 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_211(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_z);
  /* 13 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^2 \partial_y^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_220(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             (r_y * r_y) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv;
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^3 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_301(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_z);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^3 \partial_y^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_310(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x * r_y);
  /* 11 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^4}{ \partial_x^4 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_400(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return +105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 6.0 *
             (r_x * r_x) +
         3. * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0;
  /* 5 zero-valued terms not written out */
}

/*********************************/
/* 5th order gravity derivatives */
/*********************************/

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_z^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_005(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_z * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_z);
  /* 26 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^1 \partial_z^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_014(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^2 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_023(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^3 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_032(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^4 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_041(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_y^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_050(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_y * r_y * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_y);
  /* 26 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_z^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_104(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_z * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^1 \partial_z^3
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_113(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^2 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_122(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_x);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^3 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_131(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^1 \partial_y^4 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_140(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_y * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_z^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_203(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_z * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_z * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^1 \partial_z^2
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_212(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_y);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^2 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_221(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * (r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^2 \partial_y^3 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_230(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_y * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_y * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_z^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_302(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_z * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_z * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_y^1 \partial_z^1
 * }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_311(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_y * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_z);
  /* 48 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^3 \partial_y^2 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_320(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_y * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * (r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 3.0 * (r_x * r_y * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_x);
  /* 44 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^4 \partial_z^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_401(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_z) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_x * r_z) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_z);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^4 \partial_y^1 }\phi(x, y,
 * z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_410(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_y) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 6.0 * (r_x * r_x * r_y) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 3.0 *
             (r_y);
  /* 42 zero-valued terms not written out */
}

/**
 * @brief Compute \f$ \frac{\partial^5}{ \partial_x^5 }\phi(x, y, z} \f$.
 *
 * Note that r_inv = 1./sqrt(r_x^2 + r_y^2 + r_z^2)
 */
__attribute__((always_inline)) INLINE static double D_500(double r_x,
                                                          double r_y,
                                                          double r_z,
                                                          double r_inv) {
  return -945. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * r_inv * r_inv * (r_x * r_x * r_x * r_x * r_x) +
         105. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv *
             r_inv * 10.0 * (r_x * r_x * r_x) -
         15. * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * 15.0 *
             (r_x);
  /* 26 zero-valued terms not written out */
}

#endif

#endif /* SWIFT_GRAVITY_DERIVATIVE_H */
