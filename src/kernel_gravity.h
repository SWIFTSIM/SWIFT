/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#ifdef GADGET2_SOFTENING_CORRECTION
/*! Conversion factor between Plummer softening and internal softening */
#define kernel_gravity_softening_plummer_equivalent 2.8
#define kernel_gravity_softening_plummer_equivalent_inv (1. / 2.8)
#define kernel_gravity_softening_name "Gadget-2 (spline kernel)"
#else
/*! Conversion factor between Plummer softening and internal softening */
#define kernel_gravity_softening_plummer_equivalent 3.
#define kernel_gravity_softening_plummer_equivalent_inv (1. / 3.)
#define kernel_gravity_softening_name "Wendland-C2"
#endif /* GADGET2_SOFTENING_CORRECTION */

/**
 * @brief Computes the gravity softening kernel for the potential.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the spline softening length $u = x/H$.
 */
__attribute__((const)) INLINE static float kernel_grav_pot_eval(const float u) {

  float W;
#ifdef GADGET2_SOFTENING_CORRECTION
  if (u < 0.5f)
    W = -2.8f + u * u * (5.333333333333f + u * u * (6.4f * u - 9.6f));
  else
    W = -3.2f + 0.066666666667f / u +
        u * u *
            (10.666666666667f + u * (-16.f + u * (9.6f - 2.133333333333f * u)));
#else

  /* W(u) = 3u^7 - 15u^6 + 28u^5 - 21u^4 + 7u^2 - 3 */
  W = 3.f * u - 15.f;
  W = W * u + 28.f;
  W = W * u - 21.f;
  W = W * u;
  W = W * u + 7.f;
  W = W * u;
  W = W * u - 3.f;
#endif
  return W;
}

/**
 * @brief Computes the gravity softening kernel for the forces.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the spline softening length $u = x/H$.
 */
__attribute__((const)) INLINE static float kernel_grav_force_eval(
    const float u) {

  float W;
#ifdef GADGET2_SOFTENING_CORRECTION
  if (u < 0.5f)
    W = 10.6666667f + u * u * (32.f * u - 38.4f);
  else
    W = 21.3333333f - 48.f * u + 38.4f * u * u - 10.6666667f * u * u * u -
        0.06666667f / (u * u * u);
#else

  /* W(u) = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  W = 21.f * u - 90.f;
  W = W * u + 140.f;
  W = W * u - 84.f;
  W = W * u;
  W = W * u + 14.f;
#endif
  return W;
}

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

/**
 * @brief Computes the gravity softening function for potential in double
 * precision.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the softening length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval_pot_double(
    double u, double *const W) {

#ifdef GADGET2_SOFTENING_CORRECTION
  if (u < 0.5)
    *W = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
  else
    *W = -3.2 + 0.066666666667 / u +
         u * u *
             (10.666666666667 + u * (-16.0 + u * (9.6 - 2.133333333333 * u)));
#else

  /* W(u) = 3u^7 - 15u^6 + 28u^5 - 21u^4 + 7u^2 - 3 */
  *W = 3. * u - 15.;
  *W = *W * u + 28.;
  *W = *W * u - 21.;
  *W = *W * u;
  *W = *W * u + 7.;
  *W = *W * u;
  *W = *W * u - 3;
#endif
}

/**
 * @brief Computes the gravity softening function for forces in double
 * precision.
 *
 * This functions assumes 0 < u < 1.
 *
 * @param u The ratio of the distance to the softening length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval_force_double(
    double u, double *const W) {

#ifdef GADGET2_SOFTENING_CORRECTION
  if (u < 0.5)
    *W = 10.666666666667 + u * u * (32.0 * u - 38.4);
  else
    *W = 21.333333333333 - 48.0 * u + 38.4 * u * u -
         10.666666666667 * u * u * u - 0.066666666667 / (u * u * u);
#else

  /* W(u) = 21u^5 - 90u^4 + 140u^3 - 84u^2 + 14 */
  *W = 21. * u - 90.;
  *W = *W * u + 140.;
  *W = *W * u - 84.;
  *W = *W * u;
  *W = *W * u + 14.;
#endif
}
#endif /* SWIFT_GRAVITY_FORCE_CHECKS */

/************************************************/
/* Derivatives of softening kernel used for FMM */
/************************************************/

__attribute__((const)) INLINE static float D_soft_1(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* -3u^7 + 15u^6 - 28u^5 + 21u^4 - 7u^2 + 3 */
  float phi = -3.f * u + 15.f;
  phi = phi * u - 28.f;
  phi = phi * u + 21.f;
  phi = phi * u;
  phi = phi * u - 7.f;
  phi = phi * u;
  phi = phi * u + 3.f;

  return phi;
}

__attribute__((const)) INLINE static float D_soft_2(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* -21u^6 + 90u^5 - 140u^4 + 84u^3 - 14u */
  float phi = -21.f * u + 90.f;
  phi = phi * u - 140.f;
  phi = phi * u + 84.f;
  phi = phi * u;
  phi = phi * u - 14.f;
  phi = phi * u;

  return phi;
}

__attribute__((const)) INLINE static float D_soft_3(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* -105u^5 + 360u^4 - 420u^3 + 168u^2 */
  float phi = -105.f * u + 360.f;
  phi = phi * u - 420.f;
  phi = phi * u + 168.f;
  phi = phi * u;
  phi = phi * u;

  return phi;
}

__attribute__((const)) INLINE static float D_soft_4(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* -315u^4 + 720u^3 - 420u^2 */
  float phi = -315.f * u + 720.f;
  phi = phi * u - 420.f;
  phi = phi * u;
  phi = phi * u;

  return phi;
}

__attribute__((const)) INLINE static float D_soft_5(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* -315u^3 + 420u */
  float phi = -315.f * u;
  phi = phi * u + 420.f;
  phi = phi * u;

  return phi;
}

__attribute__((const)) INLINE static float D_soft_6(const float u) {

#ifdef GADGET2_SOFTENING_CORRECTION
  error("Invalid choice of softening kernel shape");
#endif

  /* 315u^2 - 1260 */
  float phi = 315 * u;
  phi = phi * u - 1260.f;

  return phi;
}

#endif /* SWIFT_KERNEL_GRAVITY_H */
