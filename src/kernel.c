/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

#include <math.h>
#include <stdio.h>

#include "kernel.h"

/**
 * @brief Test the SPH kernel function by dumping it in the interval [0,1].
 *
 * @param N number of intervals in [0,1].
 */
void SPH_kernel_dump(int N) {

  int k;
  float x, w, dw_dx;
  float x4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float w4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  // float dw_dx4[4] __attribute__ ((aligned (16)));

  for (k = 0; k <= N; k++) {
    x = ((float)k) / N;
    x4[3] = x4[2];
    x4[2] = x4[1];
    x4[1] = x4[0];
    x4[0] = x;
    kernel_deval(x, &w, &dw_dx);
    // kernel_deval_vec( (vector *)x4 , (vector *)w4 , (vector *)dw_dx4 );
    printf(" %e %e %e %e %e %e %e\n", x, w, dw_dx, w4[0], w4[1], w4[2], w4[3]);
  }
}

/**
 * @brief The Gadget-2 gravity kernel function
 *
 * @param r The distance between particles
 */
float gadget(float r) {
  float fac, h_inv, u, r2 = r * r;
  if (r >= const_epsilon)
    fac = 1.0f / (r2 * r);
  else {
    h_inv = 1. / const_epsilon;
    u = r * h_inv;
    if (u < 0.5)
      fac = const_iepsilon3 * (10.666666666667 + u * u * (32.0 * u - 38.4));
    else
      fac = const_iepsilon3 *
            (21.333333333333 - 48.0 * u + 38.4 * u * u -
             10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
  }
  return const_G * fac;
}

/**
 * @brief Test the gravity kernel function by dumping it in the interval [0,1].
 *
 * @param r_max The radius up to which the kernel is dumped.
 * @param N number of intervals in [0,1].
 */
void gravity_kernel_dump(float r_max, int N) {

  int k;
  float x, w;
  float x4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float w4[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  // float dw_dx4[4] __attribute__ ((aligned (16)));

  for (k = 1; k <= N; k++) {
    x = (r_max * k) / N;
    x4[3] = x4[2];
    x4[2] = x4[1];
    x4[1] = x4[0];
    x4[0] = x;
    kernel_grav_eval(x, &w);
    w *= const_G / (x * x * x);
    // blender_deval_vec( (vector *)x4 , (vector *)w4 , (vector *)dw_dx4 );
    printf(" %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n", x, w * x, w4[0],
           w4[1], w4[2], w4[3], gadget(x) * x);
  }
}
