/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "kernel_long_gravity.h"

/**
 * @brief Structure containing all the derivatives of the potential field
 * required for the M2L kernel
 */
struct potential_derivatives_M2L {

  /* 0th order term */
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

  /* 0th order term */
  float D_000;

  /* 1st order terms */
  float D_100, D_010, D_001;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  /* 2nd order terms */
  float D_200, D_020, D_002;
  float D_110, D_101, D_011;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  /* 3rd order terms */
  float D_300, D_030, D_003;
  float D_210, D_201;
  float D_120, D_021;
  float D_102, D_012;
  float D_111;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  /* 4th order terms */
  float D_400, D_040, D_004;
  float D_310, D_301;
  float D_130, D_031;
  float D_103, D_013;
  float D_220, D_202, D_022;
  float D_211, D_121, D_112;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  /* 5th order terms */
  float D_005, D_014, D_023;
  float D_032, D_041, D_050;
  float D_104, D_113, D_122;
  float D_131, D_140, D_203;
  float D_212, D_221, D_230;
  float D_302, D_311, D_320;
  float D_401, D_410, D_500;
#endif
};

/**
 * @brief Converts the derivatives from a distance vector to its opposite.
 *
 * From a series of tensors D_xxx(r), compute D_xxx(-r).
 * This can be computed efficiently by flipping the sign of all the odd
 * derivative terms.
 *
 * @param pot The derivatives of the potential.
 */
__attribute__((always_inline, nonnull)) INLINE static void
potential_derivatives_flip_signs(struct potential_derivatives_M2L *pot) {

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  /* 1st order terms */
  pot->D_100 = -pot->D_100;
  pot->D_010 = -pot->D_010;
  pot->D_001 = -pot->D_001;
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  /* 3rd order terms */
  pot->D_300 = -pot->D_300;
  pot->D_030 = -pot->D_030;
  pot->D_003 = -pot->D_003;
  pot->D_210 = -pot->D_210;
  pot->D_201 = -pot->D_201;
  pot->D_021 = -pot->D_021;
  pot->D_120 = -pot->D_120;
  pot->D_012 = -pot->D_012;
  pot->D_102 = -pot->D_102;
  pot->D_111 = -pot->D_111;
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  /* 5th order terms */
  pot->D_500 = -pot->D_500;
  pot->D_050 = -pot->D_050;
  pot->D_005 = -pot->D_005;
  pot->D_410 = -pot->D_410;
  pot->D_401 = -pot->D_401;
  pot->D_041 = -pot->D_041;
  pot->D_140 = -pot->D_140;
  pot->D_014 = -pot->D_014;
  pot->D_104 = -pot->D_104;
  pot->D_320 = -pot->D_320;
  pot->D_302 = -pot->D_302;
  pot->D_032 = -pot->D_032;
  pot->D_230 = -pot->D_230;
  pot->D_023 = -pot->D_023;
  pot->D_203 = -pot->D_203;
  pot->D_311 = -pot->D_311;
  pot->D_131 = -pot->D_131;
  pot->D_113 = -pot->D_113;
  pot->D_122 = -pot->D_122;
  pot->D_212 = -pot->D_212;
  pot->D_221 = -pot->D_221;
#endif
}

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
 * @param periodic Is the calculation periodic ?
 * @param r_s_inv Inverse of the long-range gravity mesh smoothing length.
 * @param pot (return) The structure containing all the derivatives.
 */
__attribute__((always_inline, nonnull)) INLINE static void
potential_derivatives_compute_M2L(const float r_x, const float r_y,
                                  const float r_z, const float r2,
                                  const float r_inv, const float eps,
                                  const int periodic, const float r_s_inv,
                                  struct potential_derivatives_M2L *pot) {

  float Dt_1;
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  float Dt_2;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  float Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  float Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  float Dt_5;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  float Dt_6;
#endif

  /* Softened case */
  if (r2 < eps * eps) {

    const float eps_inv = 1.f / eps;
    const float r = r2 * r_inv;
    const float u = r * eps_inv;

    Dt_1 = eps_inv * D_soft_1(u);
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    const float eps_inv2 = eps_inv * eps_inv;
    Dt_2 = eps_inv2 * D_soft_2(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    const float eps_inv3 = eps_inv2 * eps_inv;
    Dt_3 = eps_inv3 * D_soft_3(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    const float eps_inv4 = eps_inv3 * eps_inv;
    Dt_4 = eps_inv4 * D_soft_4(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    const float eps_inv5 = eps_inv4 * eps_inv;
    Dt_5 = eps_inv5 * D_soft_5(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
    const float eps_inv6 = eps_inv5 * eps_inv;
    Dt_6 = eps_inv6 * D_soft_6(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

    /* Un-truncated un-softened case (Newtonian potential) */
  } else if (!periodic) {

    Dt_1 = r_inv; /* 1 / r */
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    Dt_2 = -1.f * Dt_1 * r_inv; /* -1 / r^2 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    Dt_3 = -3.f * Dt_2 * r_inv; /* 3 / r^3 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    Dt_4 = -5.f * Dt_3 * r_inv; /* -15 / r^4 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    Dt_5 = -7.f * Dt_4 * r_inv; /* 105 / r^5 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
    Dt_6 = -9.f * Dt_5 * r_inv; /* -945 / r^6 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif

    /* Truncated case (long-range) */
  } else {

    /* Get the derivatives of the truncated potential */
    const float r = r2 * r_inv;
    struct chi_derivatives derivs;
    kernel_long_grav_derivatives(r, r_s_inv, &derivs);

    Dt_1 = derivs.chi_0 * r_inv;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

    /* -chi^0 r_i^2 + chi^1 r_i^1 */
    Dt_2 = derivs.chi_1 - derivs.chi_0 * r_inv;
    Dt_2 = Dt_2 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

    /* 3chi^0 r_i^3 - 3 chi^1 r_i^2 + chi^2 r_i^1 */
    Dt_3 = derivs.chi_0 * r_inv - derivs.chi_1;
    Dt_3 = Dt_3 * 3.f;
    Dt_3 = Dt_3 * r_inv + derivs.chi_2;
    Dt_3 = Dt_3 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

    /* -15chi^0 r_i^4 + 15 chi^1 r_i^3 - 6 chi^2 r_i^2  + chi^3 r_i^1 */
    Dt_4 = -derivs.chi_0 * r_inv + derivs.chi_1;
    Dt_4 = Dt_4 * 15.f;
    Dt_4 = Dt_4 * r_inv - 6.f * derivs.chi_2;
    Dt_4 = Dt_4 * r_inv + derivs.chi_3;
    Dt_4 = Dt_4 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

    /* 105chi^0 r_i^5 - 105 chi^1 r_i^4 + 45 chi^2 r_i^3 - 10 chi^3 r_i^2 +
     * chi^4 r_i^1 */
    Dt_5 = derivs.chi_0 * r_inv - derivs.chi_1;
    Dt_5 = Dt_5 * 105.f;
    Dt_5 = Dt_5 * r_inv + 45.f * derivs.chi_2;
    Dt_5 = Dt_5 * r_inv - 10.f * derivs.chi_3;
    Dt_5 = Dt_5 * r_inv + derivs.chi_4;
    Dt_5 = Dt_5 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

    /* -945chi^0 r_i^6 + 945 chi^1 r_i^5 - 420 chi^2 r_i^4 + 105 chi^3 r_i^3 -
     * 15 chi^4 r_i^2 + chi^5 r_i^1 */
    Dt_6 = -derivs.chi_0 * r_inv + derivs.chi_1;
    Dt_6 = Dt_6 * 945.f;
    Dt_6 = Dt_6 * r_inv - 420.f * derivs.chi_2;
    Dt_6 = Dt_6 * r_inv + 105.f * derivs.chi_3;
    Dt_6 = Dt_6 * r_inv - 15.f * derivs.chi_4;
    Dt_6 = Dt_6 * r_inv + derivs.chi_5;
    Dt_6 = Dt_6 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for order >5"
#endif
  }

  /* Alright, let's get the full terms */

  /* Compute some powers of (r_x / r), (r_y / r) and (r_z / r) */
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  const float rx_r = r_x * r_inv;
  const float ry_r = r_y * r_inv;
  const float rz_r = r_z * r_inv;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  const float rx_r2 = rx_r * rx_r;
  const float ry_r2 = ry_r * ry_r;
  const float rz_r2 = rz_r * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  const float rx_r3 = rx_r2 * rx_r;
  const float ry_r3 = ry_r2 * ry_r;
  const float rz_r3 = rz_r2 * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  const float rx_r4 = rx_r3 * rx_r;
  const float ry_r4 = ry_r3 * ry_r;
  const float rz_r4 = rz_r3 * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4
  const float rx_r5 = rx_r4 * rx_r;
  const float ry_r5 = ry_r4 * ry_r;
  const float rz_r5 = rz_r4 * rz_r;
#endif

  /* Get the 0th order term */
  pot->D_000 = Dt_1;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  /* 1st order derivatives */
  pot->D_100 = rx_r * Dt_2;
  pot->D_010 = ry_r * Dt_2;
  pot->D_001 = rz_r * Dt_2;
#endif

#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  Dt_2 *= r_inv;

  /* 2nd order derivatives */
  pot->D_200 = rx_r2 * Dt_3 + Dt_2;
  pot->D_020 = ry_r2 * Dt_3 + Dt_2;
  pot->D_002 = rz_r2 * Dt_3 + Dt_2;
  pot->D_110 = rx_r * ry_r * Dt_3;
  pot->D_101 = rx_r * rz_r * Dt_3;
  pot->D_011 = ry_r * rz_r * Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  Dt_3 *= r_inv;

  /* 3rd order derivatives */
  pot->D_300 = rx_r3 * Dt_4 + 3.f * rx_r * Dt_3;
  pot->D_030 = ry_r3 * Dt_4 + 3.f * ry_r * Dt_3;
  pot->D_003 = rz_r3 * Dt_4 + 3.f * rz_r * Dt_3;
  pot->D_210 = rx_r2 * ry_r * Dt_4 + ry_r * Dt_3;
  pot->D_201 = rx_r2 * rz_r * Dt_4 + rz_r * Dt_3;
  pot->D_120 = ry_r2 * rx_r * Dt_4 + rx_r * Dt_3;
  pot->D_021 = ry_r2 * rz_r * Dt_4 + rz_r * Dt_3;
  pot->D_102 = rz_r2 * rx_r * Dt_4 + rx_r * Dt_3;
  pot->D_012 = rz_r2 * ry_r * Dt_4 + ry_r * Dt_3;
  pot->D_111 = rx_r * ry_r * rz_r * Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  Dt_3 *= r_inv;
  Dt_4 *= r_inv;

  /* 4th order derivatives */
  pot->D_400 = rx_r4 * Dt_5 + 6.f * rx_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_040 = ry_r4 * Dt_5 + 6.f * ry_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_004 = rz_r4 * Dt_5 + 6.f * rz_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_310 = rx_r3 * ry_r * Dt_5 + 3.f * rx_r * ry_r * Dt_4;
  pot->D_301 = rx_r3 * rz_r * Dt_5 + 3.f * rx_r * rz_r * Dt_4;
  pot->D_130 = ry_r3 * rx_r * Dt_5 + 3.f * ry_r * rx_r * Dt_4;
  pot->D_031 = ry_r3 * rz_r * Dt_5 + 3.f * ry_r * rz_r * Dt_4;
  pot->D_103 = rz_r3 * rx_r * Dt_5 + 3.f * rz_r * rx_r * Dt_4;
  pot->D_013 = rz_r3 * ry_r * Dt_5 + 3.f * rz_r * ry_r * Dt_4;
  pot->D_220 = rx_r2 * ry_r2 * Dt_5 + rx_r2 * Dt_4 + ry_r2 * Dt_4 + Dt_3;
  pot->D_202 = rx_r2 * rz_r2 * Dt_5 + rx_r2 * Dt_4 + rz_r2 * Dt_4 + Dt_3;
  pot->D_022 = ry_r2 * rz_r2 * Dt_5 + ry_r2 * Dt_4 + rz_r2 * Dt_4 + Dt_3;
  pot->D_211 = rx_r2 * ry_r * rz_r * Dt_5 + ry_r * rz_r * Dt_4;
  pot->D_121 = ry_r2 * rx_r * rz_r * Dt_5 + rx_r * rz_r * Dt_4;
  pot->D_112 = rz_r2 * rx_r * ry_r * Dt_5 + rx_r * ry_r * Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 4

  Dt_4 *= r_inv;
  Dt_5 *= r_inv;

  /* 5th order derivatives */
  pot->D_500 = rx_r5 * Dt_6 + 10.f * rx_r3 * Dt_5 + 15.f * rx_r * Dt_4;
  pot->D_050 = ry_r5 * Dt_6 + 10.f * ry_r3 * Dt_5 + 15.f * ry_r * Dt_4;
  pot->D_005 = rz_r5 * Dt_6 + 10.f * rz_r3 * Dt_5 + 15.f * rz_r * Dt_4;
  pot->D_410 =
      rx_r4 * ry_r * Dt_6 + 6.f * rx_r2 * ry_r * Dt_5 + 3.f * ry_r * Dt_4;
  pot->D_401 =
      rx_r4 * rz_r * Dt_6 + 6.f * rx_r2 * rz_r * Dt_5 + 3.f * rz_r * Dt_4;
  pot->D_140 =
      ry_r4 * rx_r * Dt_6 + 6.f * ry_r2 * rx_r * Dt_5 + 3.f * rx_r * Dt_4;
  pot->D_041 =
      ry_r4 * rz_r * Dt_6 + 6.f * ry_r2 * rz_r * Dt_5 + 3.f * rz_r * Dt_4;
  pot->D_104 =
      rz_r4 * rx_r * Dt_6 + 6.f * rz_r2 * rx_r * Dt_5 + 3.f * rx_r * Dt_4;
  pot->D_014 =
      rz_r4 * ry_r * Dt_6 + 6.f * rz_r2 * ry_r * Dt_5 + 3.f * ry_r * Dt_4;
  pot->D_320 = rx_r3 * ry_r2 * Dt_6 + rx_r3 * Dt_5 + 3.f * rx_r * ry_r2 * Dt_5 +
               3.f * rx_r * Dt_4;
  pot->D_302 = rx_r3 * rz_r2 * Dt_6 + rx_r3 * Dt_5 + 3.f * rx_r * rz_r2 * Dt_5 +
               3.f * rx_r * Dt_4;
  pot->D_230 = ry_r3 * rx_r2 * Dt_6 + ry_r3 * Dt_5 + 3.f * ry_r * rx_r2 * Dt_5 +
               3.f * ry_r * Dt_4;
  pot->D_032 = ry_r3 * rz_r2 * Dt_6 + ry_r3 * Dt_5 + 3.f * ry_r * rz_r2 * Dt_5 +
               3.f * ry_r * Dt_4;
  pot->D_203 = rz_r3 * rx_r2 * Dt_6 + rz_r3 * Dt_5 + 3.f * rz_r * rx_r2 * Dt_5 +
               3.f * rz_r * Dt_4;
  pot->D_023 = rz_r3 * ry_r2 * Dt_6 + rz_r3 * Dt_5 + 3.f * rz_r * ry_r2 * Dt_5 +
               3.f * rz_r * Dt_4;
  pot->D_311 = rx_r3 * ry_r * rz_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_131 = ry_r3 * rx_r * rz_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_113 = rz_r3 * rx_r * ry_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_122 = rx_r * ry_r2 * rz_r2 * Dt_6 + rx_r * ry_r2 * Dt_5 +
               rx_r * rz_r2 * Dt_5 + rx_r * Dt_4;
  pot->D_212 = ry_r * rx_r2 * rz_r2 * Dt_6 + ry_r * rx_r2 * Dt_5 +
               ry_r * rz_r2 * Dt_5 + ry_r * Dt_4;
  pot->D_221 = rz_r * rx_r2 * ry_r2 * Dt_6 + rz_r * rx_r2 * Dt_5 +
               rz_r * ry_r2 * Dt_5 + rz_r * Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 5
#error "Missing implementation for orders >5"
#endif
}

/**
 * @brief Compute all the relevent derivatives of the softened and truncated
 * gravitational potential for the M2P kernel.
 *
 * For M2P, we compute the derivatives to one order higher than
 * SELF_GRAVITY_MULTIPOLE_ORDER, as these are needed for the accelerations.
 *
 * @param r_x x-component of distance vector
 * @param r_y y-component of distance vector
 * @param r_z z-component of distance vector
 * @param r2 Square norm of distance vector
 * @param r_inv Inverse norm of distance vector
 * @param eps Softening length.
 * @param periodic Is the calculation using periodic BCs?
 * @param r_s_inv The inverse of the gravity mesh-smoothing scale.
 * @param pot (return) The structure containing all the derivatives.
 */
__attribute__((always_inline, nonnull)) INLINE static void
potential_derivatives_compute_M2P(const float r_x, const float r_y,
                                  const float r_z, const float r2,
                                  const float r_inv, const float eps,
                                  const int periodic, const float r_s_inv,
                                  struct potential_derivatives_M2P *pot) {

  float Dt_1;
  float Dt_2;
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  float Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  float Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  float Dt_5;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  float Dt_6;
#endif

  /* Softened case */
  if (r2 < eps * eps) {

    const float eps_inv = 1.f / eps;
    const float r = r2 * r_inv;
    const float u = r * eps_inv;

    Dt_1 = eps_inv * D_soft_1(u);

    const float eps_inv2 = eps_inv * eps_inv;
    Dt_2 = eps_inv2 * D_soft_2(u);
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    const float eps_inv3 = eps_inv2 * eps_inv;
    Dt_3 = eps_inv3 * D_soft_3(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    const float eps_inv4 = eps_inv3 * eps_inv;
    Dt_4 = eps_inv4 * D_soft_4(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    const float eps_inv5 = eps_inv4 * eps_inv;
    Dt_5 = eps_inv5 * D_soft_5(u);
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    const float eps_inv6 = eps_inv5 * eps_inv;
    Dt_6 = eps_inv6 * D_soft_6(u);
#endif

    /* Un-truncated un-softened case (Newtonian potential) */
  } else if (!periodic) {

    Dt_1 = r_inv;               /* 1 / r */
    Dt_2 = -1.f * Dt_1 * r_inv; /* -1 / r^2 */
#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
    Dt_3 = -3.f * Dt_2 * r_inv; /* 3 / r^3 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
    Dt_4 = -5.f * Dt_3 * r_inv; /* -15 / r^4 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
    Dt_5 = -7.f * Dt_4 * r_inv; /* 105 / r^5 */
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
    Dt_6 = -9.f * Dt_5 * r_inv; /* -945 / r^6 */
#endif

    /* Truncated case (long-range) */
  } else {

    /* Get the derivatives of the truncated potential */
    const float r = r2 * r_inv;
    struct chi_derivatives derivs;
    kernel_long_grav_derivatives(r, r_s_inv, &derivs);

    Dt_1 = derivs.chi_0 * r_inv;

    /* -chi^0 r_i^2 + chi^1 r_i^1 */
    Dt_2 = derivs.chi_1 - derivs.chi_0 * r_inv;
    Dt_2 = Dt_2 * r_inv;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

    /* 3chi^0 r_i^3 - 3 chi^1 r_i^2 + chi^2 r_i^1 */
    Dt_3 = derivs.chi_0 * r_inv - derivs.chi_1;
    Dt_3 = Dt_3 * 3.f;
    Dt_3 = Dt_3 * r_inv + derivs.chi_2;
    Dt_3 = Dt_3 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

    /* -15chi^0 r_i^4 + 15 chi^1 r_i^3 - 6 chi^2 r_i^2  + chi^3 r_i^1 */
    Dt_4 = -derivs.chi_0 * r_inv + derivs.chi_1;
    Dt_4 = Dt_4 * 15.f;
    Dt_4 = Dt_4 * r_inv - 6.f * derivs.chi_2;
    Dt_4 = Dt_4 * r_inv + derivs.chi_3;
    Dt_4 = Dt_4 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

    /* 105chi^0 r_i^5 - 105 chi^1 r_i^4 + 45 chi^2 r_i^3 - 10 chi^3 r_i^2 +
     * chi^4 r_i^1 */
    Dt_5 = derivs.chi_0 * r_inv - derivs.chi_1;
    Dt_5 = Dt_5 * 105.f;
    Dt_5 = Dt_5 * r_inv + 45.f * derivs.chi_2;
    Dt_5 = Dt_5 * r_inv - 10.f * derivs.chi_3;
    Dt_5 = Dt_5 * r_inv + derivs.chi_4;
    Dt_5 = Dt_5 * r_inv;

#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

    /* -945chi^0 r_i^6 + 945 chi^1 r_i^5 - 420 chi^2 r_i^4 + 105 chi^3 r_i^3 -
     * 15 chi^4 r_i^2 + chi^5 r_i^1 */
    Dt_6 = -derivs.chi_0 * r_inv + derivs.chi_1;
    Dt_6 = Dt_6 * 945.f;
    Dt_6 = Dt_6 * r_inv - 420.f * derivs.chi_2;
    Dt_6 = Dt_6 * r_inv + 105.f * derivs.chi_3;
    Dt_6 = Dt_6 * r_inv - 15.f * derivs.chi_4;
    Dt_6 = Dt_6 * r_inv + derivs.chi_5;
    Dt_6 = Dt_6 * r_inv;

#endif
  }

  /* Alright, let's get the full terms */

  /* Compute some powers of (r_x / r), (r_y / r) and (r_z / r) */
  const float rx_r = r_x * r_inv;
  const float ry_r = r_y * r_inv;
  const float rz_r = r_z * r_inv;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0
  const float rx_r2 = rx_r * rx_r;
  const float ry_r2 = ry_r * ry_r;
  const float rz_r2 = rz_r * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1
  const float rx_r3 = rx_r2 * rx_r;
  const float ry_r3 = ry_r2 * ry_r;
  const float rz_r3 = rz_r2 * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2
  const float rx_r4 = rx_r3 * rx_r;
  const float ry_r4 = ry_r3 * ry_r;
  const float rz_r4 = rz_r3 * rz_r;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3
  const float rx_r5 = rx_r4 * rx_r;
  const float ry_r5 = ry_r4 * ry_r;
  const float rz_r5 = rz_r4 * rz_r;
#endif

  /* Get the 0th order term */
  pot->D_000 = Dt_1;

  /* 1st order derivatives */
  pot->D_100 = rx_r * Dt_2;
  pot->D_010 = ry_r * Dt_2;
  pot->D_001 = rz_r * Dt_2;

#if SELF_GRAVITY_MULTIPOLE_ORDER > 0

  Dt_2 *= r_inv;

  /* 2nd order derivatives */
  pot->D_200 = rx_r2 * Dt_3 + Dt_2;
  pot->D_020 = ry_r2 * Dt_3 + Dt_2;
  pot->D_002 = rz_r2 * Dt_3 + Dt_2;
  pot->D_110 = rx_r * ry_r * Dt_3;
  pot->D_101 = rx_r * rz_r * Dt_3;
  pot->D_011 = ry_r * rz_r * Dt_3;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 1

  Dt_3 *= r_inv;

  /* 3rd order derivatives */
  pot->D_300 = rx_r3 * Dt_4 + 3.f * rx_r * Dt_3;
  pot->D_030 = ry_r3 * Dt_4 + 3.f * ry_r * Dt_3;
  pot->D_003 = rz_r3 * Dt_4 + 3.f * rz_r * Dt_3;
  pot->D_210 = rx_r2 * ry_r * Dt_4 + ry_r * Dt_3;
  pot->D_201 = rx_r2 * rz_r * Dt_4 + rz_r * Dt_3;
  pot->D_120 = ry_r2 * rx_r * Dt_4 + rx_r * Dt_3;
  pot->D_021 = ry_r2 * rz_r * Dt_4 + rz_r * Dt_3;
  pot->D_102 = rz_r2 * rx_r * Dt_4 + rx_r * Dt_3;
  pot->D_012 = rz_r2 * ry_r * Dt_4 + ry_r * Dt_3;
  pot->D_111 = rx_r * ry_r * rz_r * Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 2

  Dt_3 *= r_inv;
  Dt_4 *= r_inv;

  /* 4th order derivatives */
  pot->D_400 = rx_r4 * Dt_5 + 6.f * rx_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_040 = ry_r4 * Dt_5 + 6.f * ry_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_004 = rz_r4 * Dt_5 + 6.f * rz_r2 * Dt_4 + 3.f * Dt_3;
  pot->D_310 = rx_r3 * ry_r * Dt_5 + 3.f * rx_r * ry_r * Dt_4;
  pot->D_301 = rx_r3 * rz_r * Dt_5 + 3.f * rx_r * rz_r * Dt_4;
  pot->D_130 = ry_r3 * rx_r * Dt_5 + 3.f * ry_r * rx_r * Dt_4;
  pot->D_031 = ry_r3 * rz_r * Dt_5 + 3.f * ry_r * rz_r * Dt_4;
  pot->D_103 = rz_r3 * rx_r * Dt_5 + 3.f * rz_r * rx_r * Dt_4;
  pot->D_013 = rz_r3 * ry_r * Dt_5 + 3.f * rz_r * ry_r * Dt_4;
  pot->D_220 = rx_r2 * ry_r2 * Dt_5 + rx_r2 * Dt_4 + ry_r2 * Dt_4 + Dt_3;
  pot->D_202 = rx_r2 * rz_r2 * Dt_5 + rx_r2 * Dt_4 + rz_r2 * Dt_4 + Dt_3;
  pot->D_022 = ry_r2 * rz_r2 * Dt_5 + ry_r2 * Dt_4 + rz_r2 * Dt_4 + Dt_3;
  pot->D_211 = rx_r2 * ry_r * rz_r * Dt_5 + ry_r * rz_r * Dt_4;
  pot->D_121 = ry_r2 * rx_r * rz_r * Dt_5 + rx_r * rz_r * Dt_4;
  pot->D_112 = rz_r2 * rx_r * ry_r * Dt_5 + rx_r * ry_r * Dt_4;
#endif
#if SELF_GRAVITY_MULTIPOLE_ORDER > 3

  Dt_4 *= r_inv;
  Dt_5 *= r_inv;

  /* 5th order derivatives */
  pot->D_500 = rx_r5 * Dt_6 + 10.f * rx_r3 * Dt_5 + 15.f * rx_r * Dt_4;
  pot->D_050 = ry_r5 * Dt_6 + 10.f * ry_r3 * Dt_5 + 15.f * ry_r * Dt_4;
  pot->D_005 = rz_r5 * Dt_6 + 10.f * rz_r3 * Dt_5 + 15.f * rz_r * Dt_4;
  pot->D_410 =
      rx_r4 * ry_r * Dt_6 + 6.f * rx_r2 * ry_r * Dt_5 + 3.f * ry_r * Dt_4;
  pot->D_401 =
      rx_r4 * rz_r * Dt_6 + 6.f * rx_r2 * rz_r * Dt_5 + 3.f * rz_r * Dt_4;
  pot->D_140 =
      ry_r4 * rx_r * Dt_6 + 6.f * ry_r2 * rx_r * Dt_5 + 3.f * rx_r * Dt_4;
  pot->D_041 =
      ry_r4 * rz_r * Dt_6 + 6.f * ry_r2 * rz_r * Dt_5 + 3.f * rz_r * Dt_4;
  pot->D_104 =
      rz_r4 * rx_r * Dt_6 + 6.f * rz_r2 * rx_r * Dt_5 + 3.f * rx_r * Dt_4;
  pot->D_014 =
      rz_r4 * ry_r * Dt_6 + 6.f * rz_r2 * ry_r * Dt_5 + 3.f * ry_r * Dt_4;
  pot->D_320 = rx_r3 * ry_r2 * Dt_6 + rx_r3 * Dt_5 + 3.f * rx_r * ry_r2 * Dt_5 +
               3.f * rx_r * Dt_4;
  pot->D_302 = rx_r3 * rz_r2 * Dt_6 + rx_r3 * Dt_5 + 3.f * rx_r * rz_r2 * Dt_5 +
               3.f * rx_r * Dt_4;
  pot->D_230 = ry_r3 * rx_r2 * Dt_6 + ry_r3 * Dt_5 + 3.f * ry_r * rx_r2 * Dt_5 +
               3.f * ry_r * Dt_4;
  pot->D_032 = ry_r3 * rz_r2 * Dt_6 + ry_r3 * Dt_5 + 3.f * ry_r * rz_r2 * Dt_5 +
               3.f * ry_r * Dt_4;
  pot->D_203 = rz_r3 * rx_r2 * Dt_6 + rz_r3 * Dt_5 + 3.f * rz_r * rx_r2 * Dt_5 +
               3.f * rz_r * Dt_4;
  pot->D_023 = rz_r3 * ry_r2 * Dt_6 + rz_r3 * Dt_5 + 3.f * rz_r * ry_r2 * Dt_5 +
               3.f * rz_r * Dt_4;
  pot->D_311 = rx_r3 * ry_r * rz_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_131 = ry_r3 * rx_r * rz_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_113 = rz_r3 * rx_r * ry_r * Dt_6 + 3.f * rx_r * ry_r * rz_r * Dt_5;
  pot->D_122 = rx_r * ry_r2 * rz_r2 * Dt_6 + rx_r * ry_r2 * Dt_5 +
               rx_r * rz_r2 * Dt_5 + rx_r * Dt_4;
  pot->D_212 = ry_r * rx_r2 * rz_r2 * Dt_6 + ry_r * rx_r2 * Dt_5 +
               ry_r * rz_r2 * Dt_5 + ry_r * Dt_4;
  pot->D_221 = rz_r * rx_r2 * ry_r2 * Dt_6 + rz_r * rx_r2 * Dt_5 +
               rz_r * ry_r2 * Dt_5 + rz_r * Dt_4;
#endif
}

#endif /* SWIFT_GRAVITY_DERIVATIVE_H */
