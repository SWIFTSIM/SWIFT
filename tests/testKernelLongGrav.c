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

const int num_tests = 1 << 10;

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(double a, double b, const char* s, const double tol,
                 const double r, const double r_s) {
  if (fabs(a - b) / fabs(a + b) > tol)
    error("Values are inconsistent: %12.15e %12.15e rel=%e (%s for r_s=%e r/r_s=%e)!",
          a, b, fabs(a - b) / fabs(a + b), s, r_s, r / r_s);
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

  for (int n = 0; n < num_tests; ++n) {

    const double r_s = exp10(4. * rand() / ((double)RAND_MAX) - 2.);
    const double r_s_inv = 1.f / r_s;

    // message("Testing r_s=%e", r_s);

    /* Loop over some radii */
    for (double i = -4; i < 1; i += 0.001) {

      /* Get a radius in the relevant range */
      const double r = exp10(i) * r_s;

      if (r > 5. * r_s) break;

      /* Compute the SWIFT expressions */
      struct chi_derivatives chi_swift;
      kernel_long_grav_derivatives((float)r, (float)r_s_inv, &chi_swift);

      /* Compute the exact expressions */
      const double one_over_sqrt_pi = M_2_SQRTPI * 0.5;
      const double u = 0.5 * r / r_s;
      const double C = one_over_sqrt_pi * exp(-u * u);

      const double chi_0 = erfc(u);
      const double chi_1 = -C / r_s;
      const double chi_2 = C * 0.5 * r * pow(r_s, -3.);
      const double chi_3 = C * 0.25 * (2. * r_s * r_s - r * r) * pow(r_s, -5.);
      const double chi_4 =
          C * 0.125 * (r * r * r - 6. * r_s * r_s * r) * pow(r_s, -7.);
      const double chi_5 =
          C * 0.0625 *
          (12. * pow(r_s, 4.) - 12. * r_s * r_s * r * r + pow(r, 4.)) *
          pow(r_s, -9.);

      check_value(chi_swift.chi_0, chi_0, "chi_0", 4e-3, r, r_s);
      check_value(chi_swift.chi_1, chi_1, "chi_1", 1e-5, r, r_s);
      check_value(chi_swift.chi_2, chi_2, "chi_2", 1e-5, r, r_s);
      check_value(chi_swift.chi_3, chi_3, "chi_3", 1e-4, r, r_s);
      check_value(chi_swift.chi_4, chi_4, "chi_4", 4e-4, r, r_s);
      check_value(chi_swift.chi_5, chi_5, "chi_5", 4e-4, r, r_s);
    }
  }

  return 0;
}
