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

/* The gravity kernel is defined as a degree 6 polynomial in the distance
   r. The resulting value should be post-multiplied with r^-3, resulting
   in a polynomial with terms ranging from r^-3 to r^3, which are
   sufficient to model both the direct potential as well as the splines
   near the origin.
   As in the hydro case, the 1/h^3 needs to be multiplied in afterwards */

/* Coefficients for the kernel. */
#define kernel_grav_name "Gadget-2 softening kernel"
#define kernel_grav_degree 6 /* Degree of the polynomial */
#define kernel_grav_ivals 2  /* Number of branches */
static const float
    kernel_grav_coeffs[(kernel_grav_degree + 1) * (kernel_grav_ivals + 1)]
    __attribute__((aligned(16))) = {32.f,
                                    -38.4f,
                                    0.f,
                                    10.66666667f,
                                    0.f,
                                    0.f,
                                    0.f, /* 0 < u < 0.5 */
                                    -10.66666667f,
                                    38.4f,
                                    -48.f,
                                    21.3333333f,
                                    0.f,
                                    0.f,
                                    -0.066666667f, /* 0.5 < u < 1 */
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f}; /* 1 < u */

/**
 * @brief Computes the gravity softening function.
 *
 * @param u The ratio of the distance to the softening length $u = x/h$.
 * @param W (return) The value of the kernel function $W(x,h)$.
 */
__attribute__((always_inline)) INLINE static void kernel_grav_eval(
    float u, float *const W) {

  /* Pick the correct branch of the kernel */
  const int ind = (int)min(u * (float)kernel_grav_ivals, kernel_grav_ivals);
  const float *const coeffs =
      &kernel_grav_coeffs[ind * (kernel_grav_degree + 1)];

  /* First two terms of the polynomial ... */
  float w = coeffs[0] * u + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_grav_degree; k++) w = u * w + coeffs[k];

  /* Return everything */
  *W = w / (u * u * u);
}

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

__attribute__((always_inline)) INLINE static void kernel_grav_eval_double(
    double u, double *const W) {

  static const double kernel_grav_coeffs_double[(kernel_grav_degree + 1) *
                                                (kernel_grav_ivals + 1)]
      __attribute__((aligned(16))) = {32.,
                                      -38.4,
                                      0.,
                                      10.66666667,
                                      0.,
                                      0.,
                                      0., /* 0 < u < 0.5 */
                                      -10.66666667,
                                      38.4,
                                      -48.,
                                      21.3333333,
                                      0.,
                                      0.,
                                      -0.066666667, /* 0.5 < u < 1 */
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.,
                                      0.}; /* 1 < u */

  /* Pick the correct branch of the kernel */
  const int ind = (int)min(u * (double)kernel_grav_ivals, kernel_grav_ivals);
  const double *const coeffs =
      &kernel_grav_coeffs_double[ind * (kernel_grav_degree + 1)];

  /* First two terms of the polynomial ... */
  double w = coeffs[0] * u + coeffs[1];

  /* ... and the rest of them */
  for (int k = 2; k <= kernel_grav_degree; k++) w = u * w + coeffs[k];

  /* Return everything */
  *W = w / (u * u * u);
}
#endif /* SWIFT_GRAVITY_FORCE_CHECKS */

#endif /* SWIFT_KERNEL_GRAVITY_H */
