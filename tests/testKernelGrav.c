/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
 * Taken from Gadget-2.0.7's forcetree.c lines 2755-2800
 *
 * @param r The distance between particles
 * @param epsilon The cut-off distance of the kernel
 */
float gadget(float r, float epsilon) {

  const float h = epsilon;
  const float h_inv = 1.f / h;

  const float u = r * h_inv;

  if (u >= 1) {
    const float r_inv = 1. / r;

    return r_inv * r_inv * r_inv;
  } else {
    if (u < 0.5)
      return h_inv * h_inv * h_inv *
             (10.666666666667 + u * u * (32.0 * u - 38.4));
    else
      return h_inv * h_inv * h_inv *
             (21.333333333333 - 48.0 * u + 38.4 * u * u -
              10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
  }
}

int main(int argc, char *argv[]) {

  const float h = 3.f;
  const float r_max = 6.f;

  for (int k = 1; k < numPoints; ++k) {

    const float r = (r_max * k) / numPoints;
    const float gadget_w = gadget(r, h);

    const float h_inv = 1.f / h;
    const float h_inv3 = h_inv * h_inv * h_inv;
    const float u = r * h_inv;

    float swift_w;
    if (r >= h) {
      swift_w = 1 / (r * r * r);

    } else {
      kernel_grav_eval(u, &swift_w);
      swift_w *= h_inv3;
    }

    if (fabsf(gadget_w - swift_w) > 1e-5 * fabsf(gadget_w)) {

      printf("%2d: r= %f h= %f u= %f Wg(r,h)= %f Ws(r,h)= %f\n", k, r, h, u,
             gadget_w, swift_w);

      printf("Invalid value ! Gadget= %e, SWIFT= %e\n", gadget_w, swift_w);
      return 1;
    }
  }

  printf("\nAll values are consistent\n");

  /* Now test the long range function */
  /* const float a_smooth = 4.5f; */

  /* for (int k = 1; k < numPoints; ++k) { */

  /*   const float r = (r_max * k) / numPoints; */

  /*   const float u = r / a_smooth; */

  /*   float swift_w; */
  /*   kernel_long_grav_eval(u, &swift_w); */

  /*   float gadget_w = erfcf(u / 2) + u * expf(-u * u / 4) / sqrtf(M_PI); */

  /*   if (fabsf(gadget_w - swift_w) > 1e-4 * fabsf(gadget_w)) { */

  /*     printf("%2d: r= %f r_lr= %f u= %f Ws(r)= %f Wg(r)= %f\n", k, r,
   * a_smooth, */
  /*            u, swift_w, gadget_w); */

  /*     printf("Invalid value ! Gadget= %e, SWIFT= %e\n", gadget_w, swift_w);
   */
  /*     return 1; */
  /*   } */
  /* } */

  return 0;
}
