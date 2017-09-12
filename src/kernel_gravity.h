/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_KERNEL_GRAVITY_H
#define SWIFT_KERNEL_GRAVITY_H

#include <math.h>

/* Includes. */
#include "inline.h"
#include "minmax.h"

/**
 * @brief Computes the gravity softening function.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the softening length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval(
    float u, float *const W) {

  /* W(u) = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  *W = 21.f * u - 90.f;
  *W = *W * u + 140.f;
  *W = *W * u - 84.f;
  *W = *W * u;
  *W = *W * u + 14.f;
}

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

/**
 * @brief Computes the gravity softening function in double precision.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the softening length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval_double(
    double u, double *const W) {

  /* W(u) = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  *W = 21. * u - 90.;
  *W = *W * u + 140.;
  *W = *W * u - 84.;
  *W = *W * u;
  *W = *W * u + 14.;
}
#endif /* SWIFT_GRAVITY_FORCE_CHECKS */

/************************************************/
/* Derivatives of softening kernel used for FMM */
/************************************************/

__attribute__((always_inline)) INLINE static float D_soft_1(float u,
                                                            float u_inv) {

  /* phi(u) = -3u^7 + 15u^6 - 28u^5 + 21u^4 - 7u^2 + 3 */
  float phi = -3.f * u + 15.f;
  phi = phi * u - 28.f;
  phi = phi * u + 21.f;
  phi = phi * u;
  phi = phi * u - 7.f;
  phi = phi * u;
  phi = phi * u + 3.f;

  return phi;
}

__attribute__((always_inline)) INLINE static float D_soft_3(float u,
                                                            float u_inv) {

  /* phi'(u)/u = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  float phi = 21.f * u - 90.f;
  phi = phi * u + 140.f;
  phi = phi * u - 84.f;
  phi = phi * u;
  phi = phi * u + 14.f;

  return phi;
}

__attribute__((always_inline)) INLINE static float D_soft_5(float u,
                                                            float u_inv) {

  /* (phi'(u)/u)'/u = -105u^3 + 360u^2 - 420u + 168 */
  float phi = -105.f * u + 360.f;
  phi = phi * u - 420.f;
  phi = phi * u + 168.f;

  return phi;
}

__attribute__((always_inline)) INLINE static float D_soft_7(float u,
                                                            float u_inv) {

  /* ((phi'(u)/u)'/u)'/u = 315u - 720 + 420u^-1 */
  return 315.f * u - 720.f + 420.f * u_inv;
}

__attribute__((always_inline)) INLINE static float D_soft_9(float u,
                                                            float u_inv) {

  /* (((phi'(u)/u)'/u)'/u)'/u = -315u^-1 + 420u^-3 */
  float phi = 420.f * u_inv;
  phi = phi * u_inv - 315.f;
  phi = phi * u_inv;

  return phi;
}

__attribute__((always_inline)) INLINE static float D_soft_11(float u,
                                                             float u_inv) {

  /* ((((phi'(u)/u)'/u)'/u)'/u)'/u = 315u^-3 - 1260u^-5 */
  float phi = -1260.f * u_inv;
  phi = phi * u_inv + 315.f;
  phi = phi * u_inv;
  phi = phi * u_inv;
  phi = phi * u_inv;

  return phi;
}

#endif /* SWIFT_KERNEL_GRAVITY_H */
