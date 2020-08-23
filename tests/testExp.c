/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "swift.h"

/* Standard includes */
#include <fenv.h>
#include <math.h>

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param tol Relative tolerance
 * @param x Value for which we tested the function (for error messages)
 */
void check_value(double a, double b, const double tol, const double x) {
  if (fabs(a - b) / fabs(a + b) > tol)
    error("Values are inconsistent: %12.15e %12.15e rel=%e (for x=%e).", a, b,
          fabs(a - b) / fabs(a + b), x);
}

int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Loop over some values */
  for (float x = 0.; x < 32.; x += 0.000001) {

    const double exact_p = exp(x);
    const double exact_n = exp(-x);
    const double swift_exp_p = optimized_expf(x);
    const double swift_exp_n = optimized_expf(-x);

    check_value(exact_p, swift_exp_p, 2.1e-6, x);
    check_value(exact_n, swift_exp_n, 2.1e-6, x);
  }

  return 0;
}
