/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2021 Willem Elbers (whe@willemelbers.com)
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
#include "log.h"
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

/**
 * @brief Check that a and b are consistent (up to some absolute error)
 *
 * @param a First value
 * @param b Second value
 * @param tol Relative tolerance
 * @param x Value for which we tested the function (for error messages)
 */
void check_value_abs(double a, double b, const double tol, const double x) {
  if (fabs(a - b) > tol)
    error("Values are inconsistent: %12.15e %12.15e abs=%e (for x=%e).", a, b,
          fabs(a - b), x);
}

int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Loop over some values */
  for (float x = 0.; x < 32; x += 0.000001) {

    const double exact_p = x;
    const double exact_n = -x;
    const double swift_log_p = optimized_logf(exp(x));
    const double swift_log_n = optimized_logf(exp(-x));
    const double swift_log10_p = optimized_log10f(exp10(x));
    const double swift_log10_n = optimized_log10f(exp10(-x));

    /* Check the absolute error */
    check_value_abs(exact_p, swift_log_p, 5e-3, x);
    check_value_abs(exact_n, swift_log_n, 5e-3, x);
    check_value_abs(exact_p, swift_log10_p, 5e-3, x);
    check_value_abs(exact_n, swift_log10_n, 5e-3, x);

    /* Check the relative error in the domain where it is small */
    if (x >= 0.3) {
      check_value(exact_p, swift_log_p, 5e-3, x);
      check_value(exact_n, swift_log_n, 5e-3, x);
      check_value(exact_p, swift_log10_p, 5e-3, x);
      check_value(exact_n, swift_log10_n, 5e-3, x);
    }
  }

  message("Success.");

  return 0;
}
