/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *                    James Willis (james.s.willis@durham.ac.uk)
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
#include "const.h"
#include "kernel_gravity.h"
#include "kernel_long_gravity.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#define numPoints (1 << 7)

/**
 * @brief The Gadget-2 gravity kernel function
 *
 * @param r The distance between particles
 * @param h The cut-off distance of the kernel
 */
float gadget(float r, float h) {
  float fac;
  const float r2 = r * r;
  if (r >= h)
    fac = 1.0f / (r2 * r);
  else {
    const float h_inv = 1. / h;
    const float h_inv3 = h_inv * h_inv * h_inv;
    const float u = r * h_inv;
    if (u < 0.5)
      fac = h_inv3 * (10.666666666667 + u * u * (32.0 * u - 38.4));
    else
      fac =
          h_inv3 * (21.333333333333 - 48.0 * u + 38.4 * u * u -
                    10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
  }
  return fac;
}

int main() {

  const float h = 3.f;
  const float r_max = 6.f;

  for (int k = 1; k < numPoints; ++k) {

    const float r = (r_max * k) / numPoints;

    const float u = r / h;

    const float gadget_w = gadget(r, h);

    float swift_w;
    if (u < 1.) {
      kernel_grav_eval(u, &swift_w);
      swift_w *= (1 / (h * h * h));
    } else {
      swift_w = 1 / (r * r * r);
    }

    if (fabsf(gadget_w - swift_w) > 2e-7) {

      printf("%2d: r= %f h= %f u= %f Wg(r,h)= %f Ws(r,h)= %f\n", k, r, h, u,
             gadget_w, swift_w);

      printf("Invalid value ! Gadget= %e, SWIFT= %e\n", gadget_w, swift_w);
      return 1;
    }
  }

  printf("\nAll values are consistent\n");

  /* Now test the long range function */
  const float a_smooth = const_gravity_a_smooth;

  for (int k = 1; k < numPoints; ++k) {

    const float r = (r_max * k) / numPoints;

    const float u = r / a_smooth;

    float swift_w;
    kernel_long_grav_eval(u, &swift_w);

    float gadget_w = erfc(u / 2) + u * exp(-u * u / 4) / sqrt(M_PI);

    if (fabsf(gadget_w - swift_w) > 2e-7) {

      printf("%2d: r= %f r_lr= %f u= %f Ws(r)= %f Wg(r)= %f\n", k, r, a_smooth,
             u, swift_w, gadget_w);

      printf("Invalid value ! Gadget= %e, SWIFT= %e\n", gadget_w, swift_w);
      return 1;
    }
  }

  return 0;
}
