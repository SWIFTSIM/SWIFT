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
#include "kernel_gravity.h"

/**
 * @brief Structure containing all the derivatives of the potential field
 * required for the M2L kernel
 */
struct potential_derivatives_M2L {

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

/**
 * @brief Structure containing all the derivatives of the potential field
 * required for the M2P kernel
 */
struct potential_derivatives_M2P {

  /* 1st order terms */
  float D_100, D_010, D_001;

  /* 3rd order terms */
  float D_300, D_030, D_003;
  float D_210, D_201;
  float D_120, D_021;
  float D_102, D_012;
  float D_111;
};

/**
 * @brief Compute all the relevent derivatives of the softened and truncated
 * gravitational potential for the M2L kernel.
 *
 * @param r_x x-component of distance vector
 * @param r_y y-component of distance vector
 * @param r_z z-component of distance vector
 * @param r2 Square norm of distance vector
 * @param r_inv Inverse norm of distance vector
 * @param eps Softening length.
 * @param eps_inv Inverse of softening length.
 * @param pot (return) The structure containing all the derivatives.
 */
__attribute__((always_inline)) INLINE static void
compute_potential_derivatives_M2L(float r_x, float r_y, float r_z, float r2,
                                  float r_inv, float eps, float eps_inv,
                                  struct potential_derivatives_M2L *pot) {

  float Dt_1;
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  float Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  float Dt_5;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  float Dt_7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  float Dt_9;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  float Dt_11;
#endif

  /* Un-softened case */
  if (r2 > eps * eps) {

    Dt_1 = r_inv;
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    const float r_inv2 = r_inv * r_inv;
    Dt_3 = -1.f * Dt_1 * r_inv2; /* -1 / r^3 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    Dt_5 = -3.f * Dt_3 * r_inv2; /* 3 / r^5 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    Dt_7 = -5.f * Dt_5 * r_inv2; /* -15 / r^7 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    Dt_9 = -7.f * Dt_7 * r_inv2; /* 105 / r^9 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
    Dt_11 = -9.f * Dt_9 * r_inv2; /* -945 / r^11 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

  } else {
    const float r = r2 * r_inv;
    const float u = r * eps_inv;
    const float u_inv = r_inv * eps;

    Dt_1 = eps_inv * D_soft_1(u, u_inv);
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    const float eps_inv2 = eps_inv * eps_inv;
    const float eps_inv3 = eps_inv * eps_inv2;
    Dt_3 = -eps_inv3 * D_soft_3(u, u_inv);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    const float eps_inv5 = eps_inv3 * eps_inv2;
    Dt_5 = eps_inv5 * D_soft_5(u, u_inv);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    const float eps_inv7 = eps_inv5 * eps_inv2;
    Dt_7 = -eps_inv7 * D_soft_7(u, u_inv);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    const float eps_inv9 = eps_inv7 * eps_inv2;
    Dt_9 = eps_inv9 * D_soft_9(u, u_inv);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
    const float eps_inv11 = eps_inv9 * eps_inv2;
    Dt_11 = -eps_inv11 * D_soft_11(u, u_inv);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif
  }

/* Alright, let's get the full terms */

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
  pot->D_000 = Dt_1;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* 1st order derivatives */
  pot->D_100 = r_x * Dt_3;
  pot->D_010 = r_y * Dt_3;
  pot->D_001 = r_z * Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* 2nd order derivatives */
  pot->D_200 = r_x2 * Dt_5 + Dt_3;
  pot->D_020 = r_y2 * Dt_5 + Dt_3;
  pot->D_002 = r_z2 * Dt_5 + Dt_3;
  pot->D_110 = r_x * r_y * Dt_5;
  pot->D_101 = r_x * r_z * Dt_5;
  pot->D_011 = r_y * r_z * Dt_5;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* 3rd order derivatives */
  pot->D_300 = r_x3 * Dt_7 + 3.f * r_x * Dt_5;
  pot->D_030 = r_y3 * Dt_7 + 3.f * r_y * Dt_5;
  pot->D_003 = r_z3 * Dt_7 + 3.f * r_z * Dt_5;
  pot->D_210 = r_x2 * r_y * Dt_7 + r_y * Dt_5;
  pot->D_201 = r_x2 * r_z * Dt_7 + r_z * Dt_5;
  pot->D_120 = r_y2 * r_x * Dt_7 + r_x * Dt_5;
  pot->D_021 = r_y2 * r_z * Dt_7 + r_z * Dt_5;
  pot->D_102 = r_z2 * r_x * Dt_7 + r_x * Dt_5;
  pot->D_012 = r_z2 * r_y * Dt_7 + r_y * Dt_5;
  pot->D_111 = r_x * r_y * r_z * Dt_7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  /* 4th order derivatives */
  pot->D_400 = r_x4 * Dt_9 + 6.f * r_x2 * Dt_7 + 3.f * Dt_5;
  pot->D_040 = r_y4 * Dt_9 + 6.f * r_y2 * Dt_7 + 3.f * Dt_5;
  pot->D_004 = r_z4 * Dt_9 + 6.f * r_z2 * Dt_7 + 3.f * Dt_5;
  pot->D_310 = r_x3 * r_y * Dt_9 + 3.f * r_x * r_y * Dt_7;
  pot->D_301 = r_x3 * r_z * Dt_9 + 3.f * r_x * r_z * Dt_7;
  pot->D_130 = r_y3 * r_x * Dt_9 + 3.f * r_y * r_x * Dt_7;
  pot->D_031 = r_y3 * r_z * Dt_9 + 3.f * r_y * r_z * Dt_7;
  pot->D_103 = r_z3 * r_x * Dt_9 + 3.f * r_z * r_x * Dt_7;
  pot->D_013 = r_z3 * r_y * Dt_9 + 3.f * r_z * r_y * Dt_7;
  pot->D_220 = r_x2 * r_y2 * Dt_9 + r_x2 * Dt_7 + r_y2 * Dt_7 + Dt_5;
  pot->D_202 = r_x2 * r_z2 * Dt_9 + r_x2 * Dt_7 + r_z2 * Dt_7 + Dt_5;
  pot->D_022 = r_y2 * r_z2 * Dt_9 + r_y2 * Dt_7 + r_z2 * Dt_7 + Dt_5;
  pot->D_211 = r_x2 * r_y * r_z * Dt_9 + r_y * r_z * Dt_7;
  pot->D_121 = r_y2 * r_x * r_z * Dt_9 + r_x * r_z * Dt_7;
  pot->D_112 = r_z2 * r_x * r_y * Dt_9 + r_x * r_y * Dt_7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* 5th order derivatives */
  pot->D_500 = r_x5 * Dt_11 + 10.f * r_x3 * Dt_9 + 15.f * r_x * Dt_7;
  pot->D_050 = r_y5 * Dt_11 + 10.f * r_y3 * Dt_9 + 15.f * r_y * Dt_7;
  pot->D_005 = r_z5 * Dt_11 + 10.f * r_z3 * Dt_9 + 15.f * r_z * Dt_7;
  pot->D_410 = r_x4 * r_y * Dt_11 + 6.f * r_x2 * r_y * Dt_9 + 3.f * r_y * Dt_7;
  pot->D_401 = r_x4 * r_z * Dt_11 + 6.f * r_x2 * r_z * Dt_9 + 3.f * r_z * Dt_7;
  pot->D_140 = r_y4 * r_x * Dt_11 + 6.f * r_y2 * r_x * Dt_9 + 3.f * r_x * Dt_7;
  pot->D_041 = r_y4 * r_z * Dt_11 + 6.f * r_y2 * r_z * Dt_9 + 3.f * r_z * Dt_7;
  pot->D_104 = r_z4 * r_x * Dt_11 + 6.f * r_z2 * r_x * Dt_9 + 3.f * r_x * Dt_7;
  pot->D_014 = r_z4 * r_y * Dt_11 + 6.f * r_z2 * r_y * Dt_9 + 3.f * r_y * Dt_7;
  pot->D_320 = r_x3 * r_y2 * Dt_11 + r_x3 * Dt_9 + 3.f * r_x * r_y2 * Dt_9 +
               3.f * r_x * Dt_7;
  pot->D_302 = r_x3 * r_z2 * Dt_11 + r_x3 * Dt_9 + 3.f * r_x * r_z2 * Dt_9 +
               3.f * r_x * Dt_7;
  pot->D_230 = r_y3 * r_x2 * Dt_11 + r_y3 * Dt_9 + 3.f * r_y * r_x2 * Dt_9 +
               3.f * r_y * Dt_7;
  pot->D_032 = r_y3 * r_z2 * Dt_11 + r_y3 * Dt_9 + 3.f * r_y * r_z2 * Dt_9 +
               3.f * r_y * Dt_7;
  pot->D_203 = r_z3 * r_x2 * Dt_11 + r_z3 * Dt_9 + 3.f * r_z * r_x2 * Dt_9 +
               3.f * r_z * Dt_7;
  pot->D_023 = r_z3 * r_y2 * Dt_11 + r_z3 * Dt_9 + 3.f * r_z * r_y2 * Dt_9 +
               3.f * r_z * Dt_7;
  pot->D_311 = r_x3 * r_y * r_z * Dt_11 + 3.f * r_x * r_y * r_z * Dt_9;
  pot->D_131 = r_y3 * r_x * r_z * Dt_11 + 3.f * r_x * r_y * r_z * Dt_9;
  pot->D_113 = r_z3 * r_x * r_y * Dt_11 + 3.f * r_x * r_y * r_z * Dt_9;
  pot->D_122 = r_x * r_y2 * r_z2 * Dt_11 + r_x * r_y2 * Dt_9 +
               r_x * r_z2 * Dt_9 + r_x * Dt_7;
  pot->D_212 = r_y * r_x2 * r_z2 * Dt_11 + r_y * r_x2 * Dt_9 +
               r_y * r_z2 * Dt_9 + r_y * Dt_7;
  pot->D_221 = r_z * r_x2 * r_y2 * Dt_11 + r_z * r_x2 * Dt_9 +
               r_z * r_y2 * Dt_9 + r_z * Dt_7;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for orders >5"
#endif
}

/**
 * @brief Compute all the relevent derivatives of the softened and truncated
 * gravitational potential for the M2P kernel.
 *
 * @param r_x x-component of distance vector
 * @param r_y y-component of distance vector
 * @param r_z z-component of distance vector
 * @param r2 Square norm of distance vector
 * @param r_inv Inverse norm of distance vector
 * @param eps Softening length.
 * @param eps_inv Inverse of softening length.
 * @param pot (return) The structure containing all the derivatives.
 */
__attribute__((always_inline)) INLINE static void
compute_potential_derivatives_M2P(float r_x, float r_y, float r_z, float r2,
                                  float r_inv, float eps, float eps_inv,
                                  struct potential_derivatives_M2P *pot) {

  float Dt_1;
  float Dt_3;
  float Dt_5;
  float Dt_7;

  /* Un-softened case */
  if (r2 > eps * eps) {

    const float r_inv2 = r_inv * r_inv;

    Dt_1 = r_inv;
    Dt_3 = -1.f * Dt_1 * r_inv2; /* -1 / r^3 */
    Dt_5 = -3.f * Dt_3 * r_inv2; /* 3 / r^5 */
    Dt_7 = -5.f * Dt_5 * r_inv2; /* -15 / r^7 */

  } else {

    const float r = r2 * r_inv;
    const float u = r * eps_inv;
    const float u_inv = r_inv * eps;
    const float eps_inv2 = eps_inv * eps_inv;
    const float eps_inv3 = eps_inv * eps_inv2;
    const float eps_inv5 = eps_inv3 * eps_inv2;
    const float eps_inv7 = eps_inv5 * eps_inv2;

    Dt_1 = eps_inv * D_soft_1(u, u_inv);
    Dt_3 = -eps_inv3 * D_soft_3(u, u_inv);
    Dt_5 = eps_inv5 * D_soft_5(u, u_inv);
    Dt_7 = -eps_inv7 * D_soft_7(u, u_inv);
  }

  /* Compute some powers of r_x, r_y and r_z */
  const float r_x2 = r_x * r_x;
  const float r_y2 = r_y * r_y;
  const float r_z2 = r_z * r_z;
  const float r_x3 = r_x2 * r_x;
  const float r_y3 = r_y2 * r_y;
  const float r_z3 = r_z2 * r_z;

  /* 1st order derivatives */
  pot->D_100 = r_x * Dt_3;
  pot->D_010 = r_y * Dt_3;
  pot->D_001 = r_z * Dt_3;

  /* 3rd order derivatives */
  pot->D_300 = r_x3 * Dt_7 + 3.f * r_x * Dt_5;
  pot->D_030 = r_y3 * Dt_7 + 3.f * r_y * Dt_5;
  pot->D_003 = r_z3 * Dt_7 + 3.f * r_z * Dt_5;
  pot->D_210 = r_x2 * r_y * Dt_7 + r_y * Dt_5;
  pot->D_201 = r_x2 * r_z * Dt_7 + r_z * Dt_5;
  pot->D_120 = r_y2 * r_x * Dt_7 + r_x * Dt_5;
  pot->D_021 = r_y2 * r_z * Dt_7 + r_z * Dt_5;
  pot->D_102 = r_z2 * r_x * Dt_7 + r_x * Dt_5;
  pot->D_012 = r_z2 * r_y * Dt_7 + r_y * Dt_5;
  pot->D_111 = r_x * r_y * r_z * Dt_7;
}

#endif /* SWIFT_GRAVITY_DERIVATIVE_H */
