/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

#include "kernel_hydro.h"

#include <math.h>
#include <stdio.h>

/**
 * @brief Test the SPH kernel function by dumping it in the interval [0,1].
 *
 * @param N number of intervals in [0,1].
 */
void hydro_kernel_dump(int N) {

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
