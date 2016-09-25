/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

// Standard includes.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

// Local includes.
#include "cbrt.h"
#include "clocks.h"
#include "error.h"

int main(int argc, char *argv[]) {

  /* Some constants for this test. */
  const int num_vals = 20000000;
  const float range_min = -10.0f;
  const float range_max = 10.0f;
  const float err_rel_tol = 1e-6;
  message("executing %i runs of each command.", num_vals);

  /* Create and fill an array of floats. */
  float *data = (float *)malloc(sizeof(float) * num_vals);
  for (int k = 0; k < num_vals; k++) {
    data[k] = (float)rand() / RAND_MAX;
    data[k] = (1.0f - data[k]) * range_min + data[k] * range_max;
  }

  /* First run just checks for correctnes. */
  for (int k = 0; k < num_vals; k++) {
    const double exact = cbrt(data[k]);  // computed in doule just to be sure.
    const float ours = 1.0f / icbrtf(data[k]);
    const float err_abs = fabsf(exact - ours);
    if (err_abs * fabsf(exact) > err_rel_tol) {
      error("failed for x = %.8e, exact = %.8e, ours = %.8e, err = %.3e.\n",
            data[k], exact, ours, err_abs);
    }
  }

  /* Second run to check the speed of the inverse cube root. */
  float acc_exact = 0.0f;
  ticks tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_exact += 1.0f / cbrtf(data[k]);
  }
  message("1.0f / cbrtf took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  float acc_ours = 0.0f;
  ticks tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_ours += icbrtf(data[k]);
  }
  message("icbrtf took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(),
          acc_ours);

  /* Third run to check the speed of the cube root. */
  acc_exact = 0.0f;
  tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_exact += cbrtf(data[k]);
  }
  message("cbrtf took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  acc_ours = 0.0f;
  tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    const float temp = icbrtf(data[k]);
    acc_ours += data[k] * temp * temp;
  }
  message("x * icbrtf^2 took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(),
          acc_ours);

  /* Fourth run to check the speed of (.)^(2/3). */
  acc_exact = 0.0f;
  tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    const float temp = cbrtf(data[k]);
    acc_exact += temp * temp;
  }
  message("cbrtf^2 took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  acc_ours = 0.0f;
  tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_ours += data[k] * icbrtf(data[k]);
  }
  message("x * icbrtf took %.3f %s (acc = %.8e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(),
          acc_ours);

  return 0;
}
