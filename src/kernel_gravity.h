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

__attribute__((always_inline)) INLINE static double D_soft_0(double u) {

  /* phi(u) = -3u^7 + 15u^6 - 28u^5 + 21u^4 - 7u^2 + 3 */
  double phi = -3. * u + 15.;
  phi = phi * u - 28.;
  phi = phi * u + 21.;
  phi = phi * u;
  phi = phi * u - 7.;
  phi = phi * u;
  phi = phi * u + 3.;

  return phi;
}

__attribute__((always_inline)) INLINE static double D_soft_1(double u) {

  /* phi'(u)/u = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  double phi = 21. * u - 90.;
  phi = phi * u + 140.;
  phi = phi * u - 84.;
  phi = phi * u;
  phi = phi * u + 14.;

  return phi;
}

__attribute__((always_inline)) INLINE static double D_soft_2(double u) {

  /* (phi'(u)/u)'/u = -105u^3 + 360u^2 - 420u + 168 */
  double phi = -105. * u + 360.;
  phi = phi * u - 420.;
  phi = phi * u + 168.;

  return phi;
}

__attribute__((always_inline)) INLINE static double D_soft_3(double u) {

  /* ((phi'(u)/u)'/u)'/u = 315u - 720 + 420/u */
  return 315. * u - 720. + 420. / u;
}

#endif /* SWIFT_KERNEL_GRAVITY_H */
