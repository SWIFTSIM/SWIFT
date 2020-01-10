/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include "../config.h"

#include "approx_math.h"
#include "vector.h"

#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {

  const int numPoints = 60000;

  for (int i = 0; i < numPoints; ++i) {

    const float x = 0.6f * (i / (float)numPoints) - 0.3f;
    const float exp_correct = expf(x);
    const float exp_approx = approx_expf(x);

    const float abs = fabs(exp_correct - exp_approx);
    const float rel =
        0.5f * fabs(exp_correct - exp_approx) / fabs(exp_correct + exp_approx);

    int error = 0;

    if (abs > 3e-6 && fabsf(x) <= 0.2) {
      printf("Absolute difference too large !\n");
      error = 1;
    }
    if (abs > 3e-7 && fabsf(x) <= 0.1) {
      printf("Absolute difference too large !\n");
      error = 2;
    }

    if (rel > 1e-6 && fabsf(x) <= 0.2) {
      printf("Relative difference too large !\n");
      error = 3;
    }
    if (rel > 3e-7 && fabsf(x) <= 0.1) {
      printf("Relative difference too large !\n");
      error = 4;
    }

    if (error > 0) {
      printf("%2d/%d: x= %f exp(x)= %e approx_exp(x)=%e abs=%e rel=%e\n", i, error, x,
             exp_correct, exp_approx, abs, rel);
      return 1;
    }
  }

  printf("\nAll values are consistent\n");

  return 0;
}
