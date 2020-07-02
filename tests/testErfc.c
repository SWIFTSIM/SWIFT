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
 * Compute erfcf(u) using eq. 7.1.26 of
 * Abramowitz & Stegun, 1972.
 *
 * This has a *relative* error of less than 3.4e-3 over
 * the range of interest (0 < x < 5)
 */
float optimized_erfcf(const float x) {

  const float x2 = x * x;
  const float exp_x2 = expf(-x2);

  const float t = 1.f / (1.f + 0.3275911f * x);

  const float a1 = 0.254829592f;
  const float a2 = -0.284496736f;
  const float a3 = 1.421413741f;
  const float a4 = -1.453152027;
  const float a5 = 1.061405429f;

  /* a1 * t + a2 * t^2 + a3 * t^3 + a4 * t^4 + a5 * t^5 */
  float a = a5 * t + a4;
  a = a * t + a3;
  a = a * t + a2;
  a = a * t + a1;
  a = a * t;

  return a * exp_x2;
}

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param rel_tol Relative tolerance
 * @param abs_tol Absolute tolerance
 * @param x Value for which we tested the function (for error messages)
 */
void check_value(const double a, const double b, const double rel_tol,
                 const double abs_tol, const double x) {

  if (fabs(a - b) / fabs(a + b) > rel_tol)
    error("Values are inconsistent: %12.15e %12.15e rel=%e (for x=%e).", a, b,
          fabs(a - b) / fabs(a + b), x);
  if (fabs(a - b) > abs_tol)
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
  for (float x = 0.f; x < 5.f; x += 0.000001f) {

    const double exact = erfc(x);
    const double swift_erfcf = optimized_erfcf(x);

    check_value(exact, swift_erfcf, 3.358e-3, 6.1e-7, x);
  }

  return 0;
}
