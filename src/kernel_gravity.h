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
#include "const.h"
#include "inline.h"
#include "minmax.h"
#include "vector.h"

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

  /* W(u) = 21u^6 - 90u^5 + 140u^4 - 84u^3 + 14u */
  *W = 21.f * u - 90.f;
  *W = *W * u + 140.f;
  *W = *W * u - 84.f;
  *W = *W * u;
  *W = *W * u + 14.f;
  *W = *W * u;
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

  /* W(u) = 21u^6 - 90u^5 + 140u^4 - 84u^3 + 14u */
  *W = 21. * u - 90.;
  *W = *W * u + 140.;
  *W = *W * u - 84.;
  *W = *W * u;
  *W = *W * u + 14.;
  *W = *W * u;
}
#endif /* SWIFT_GRAVITY_FORCE_CHECKS */

#endif /* SWIFT_KERNEL_GRAVITY_H */
