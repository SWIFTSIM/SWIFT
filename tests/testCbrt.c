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

/* Config parameters. */
#include "../config.h"

/* Standard includes. */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local includes. */
#include "cbrt.h"
#include "clocks.h"
#include "error.h"

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FP-exceptions */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Some constants for this test. */
  const int num_vals = 200000000;
  const float range_min = -1e6f;
  const float range_max = 1e6f;
  const float err_rel_tol = 1e-7;
  message("executing %i runs of each command.", num_vals);

  /* Create and fill an array of floats. */
  float *data = (float *)malloc(sizeof(float) * num_vals);
  for (int k = 0; k < num_vals; k++) {
    data[k] = (float)rand() / RAND_MAX;
    data[k] = (1.0f - data[k]) * range_min + data[k] * range_max;
    if (data[k] == 0.f) k--; /* Skip 0 to avoid spurious mistakes */
  }

  /* First run just checks for correctness. */
  for (int k = 0; k < num_vals; k++) {
    const double exact = cbrt(data[k]);  // computed in double just to be sure.
    const float ours = 1.0f / icbrtf(data[k]);
    const float err_abs = fabs(exact - ours);
    const float err_rel = 0.5f * fabs(exact - ours) / (exact + ours);

    if (err_rel > err_rel_tol && data[k] != 0.f)
      error(
          "failed for x = %.8e, exact = %.8e, ours = %.8e, err = %.3e, rel = "
          "%.3e",
          data[k], exact, ours, err_abs, err_rel);
  }

  /* Second run to check the speed of the inverse cube root. */
  double acc_exact = 0.0f;
  ticks tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_exact += 1.0f / cbrtf(data[k]);
  }
  message("1.0f / cbrtf took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  double acc_ours = 0.0;
  ticks tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_ours += icbrtf(data[k]);
  }
  message("icbrtf       took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(), acc_ours);

  /* Third run to check the speed of the cube root. */
  acc_exact = 0.0f;
  tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_exact += cbrtf(data[k]);
  }
  message("cbrtf        took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  acc_ours = 0.0f;
  tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    const float temp = icbrtf(data[k]);
    acc_ours += data[k] * temp * temp;
  }
  message("x * icbrtf^2 took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(), acc_ours);

  /* Fourth run to check the speed of (.)^(2/3). */
  acc_exact = 0.0f;
  tic_exact = getticks();
  for (int k = 0; k < num_vals; k++) {
    const float temp = cbrtf(data[k]);
    acc_exact += temp * temp;
  }
  message("cbrtf^2      took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_exact), clocks_getunit(),
          acc_exact);

  acc_ours = 0.0f;
  tic_ours = getticks();
  for (int k = 0; k < num_vals; k++) {
    acc_ours += data[k] * icbrtf(data[k]);
  }
  message("x * icbrtf   took %9.3f %s (acc = %18.11e).",
          clocks_from_ticks(getticks() - tic_ours), clocks_getunit(), acc_ours);

  free(data);
  return 0;
}
