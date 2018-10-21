/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swift.h"

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(float a, float b, const char* s) {
  if (fabsf(a - b) / fabsf(a + b) > 1.e-6f)
    error("Values are inconsistent: %12.15e %12.15e rel=%e (%s)!", a, b,
          fabsf(a - b) / fabsf(a + b), s);
}

/**
 * @brief Check that the pre-defined adiabatic index constants contain correct
 * values
 */
void check_constants(void) {
  float val;

  val = 0.5 * (hydro_gamma + 1.0f) / hydro_gamma;
  check_value(val, hydro_gamma_plus_one_over_two_gamma, "(gamma+1)/(2 gamma)");

  val = 0.5 * (hydro_gamma - 1.0f) / hydro_gamma;
  check_value(val, hydro_gamma_minus_one_over_two_gamma, "(gamma-1)/(2 gamma)");

  val = (hydro_gamma - 1.0f) / (hydro_gamma + 1.0f);
  check_value(val, hydro_gamma_minus_one_over_gamma_plus_one,
              "(gamma-1)/(gamma+1)");

  val = 2.0f / (hydro_gamma + 1.0f);
  check_value(val, hydro_two_over_gamma_plus_one, "2/(gamma+1)");

  val = 2.0f / (hydro_gamma - 1.0f);
  check_value(val, hydro_two_over_gamma_minus_one, "2/(gamma-1)");

  val = 0.5f * (hydro_gamma - 1.0f);
  check_value(val, hydro_gamma_minus_one_over_two, "(gamma-1)/2");

  val = 2.0f * hydro_gamma / (hydro_gamma - 1.0f);
  check_value(val, hydro_two_gamma_over_gamma_minus_one, "(2 gamma)/(gamma-1)");

  val = 1.0f / hydro_gamma;
  check_value(val, hydro_one_over_gamma, "1/gamma");
}

/**
 * @brief Check that the adiabatic index power functions return the correct
 * values
 */
void check_functions(float x) {

  float val_a, val_b;
  const double xx = x;

#if defined(HYDRO_GAMMA_5_3)
#define hydro_gamma_d (5. / 3.)
#elif defined(HYDRO_GAMMA_7_5)
#define hydro_gamma_d (7. / 5.)
#elif defined(HYDRO_GAMMA_4_3)
#define hydro_gamma_d (4. / 3.)
#elif defined(HYDRO_GAMMA_2_1)
#define hydro_gamma_d (2. / 1.)
#else
#error "Need to choose an adiabatic index!"
#endif

  val_a = pow(xx, hydro_gamma_d);
  val_b = pow_gamma(x);
  check_value(val_a, val_b, "x^gamma");

  val_a = pow(xx, hydro_gamma_d - 1.0);
  val_b = pow_gamma_minus_one(x);
  check_value(val_a, val_b, "x^(gamma - 1)");

  val_a = pow(xx, -(hydro_gamma_d - 1.0));
  val_b = pow_minus_gamma_minus_one(x);
  check_value(val_a, val_b, "x^(-(gamma - 1))");

  val_a = pow(xx, -hydro_gamma_d);
  val_b = pow_minus_gamma(x);
  check_value(val_a, val_b, "x^(-gamma)");

  val_a = pow(xx, 2.0 / (hydro_gamma_d - 1.0));
  val_b = pow_two_over_gamma_minus_one(x);
  check_value(val_a, val_b, "x^(2/(gamma-1))");

  val_a = pow(xx, 2.0 * hydro_gamma_d / (hydro_gamma_d - 1.0));
  val_b = pow_two_gamma_over_gamma_minus_one(x);
  check_value(val_a, val_b, "x^((2 gamma)/(gamma-1))");

  val_a = pow(xx, (hydro_gamma_d - 1.0) / (2.0 * hydro_gamma_d));
  val_b = pow_gamma_minus_one_over_two_gamma(x);
  check_value(val_a, val_b, "x^((gamma-1)/(2 gamma))");

  val_a = pow(xx, -(hydro_gamma_d + 1.0) / (2.0 * hydro_gamma_d));
  val_b = pow_minus_gamma_plus_one_over_two_gamma(x);
  check_value(val_a, val_b, "x^(-(gamma+1)/(2 gamma))");

  val_a = pow(xx, 1.0 / hydro_gamma_d);
  val_b = pow_one_over_gamma(x);
  check_value(val_a, val_b, "x^(1/gamma)");

  val_a = pow(xx, 3. * hydro_gamma_d - 2.);
  val_b = pow_three_gamma_minus_two(x);
  check_value(val_a, val_b, "x^(3gamma - 2)");

  val_a = pow(xx, (3. * hydro_gamma_d - 5.) / 2.);
  val_b = pow_three_gamma_minus_five_over_two(x);
  check_value(val_a, val_b, "x^((3gamma - 5)/2)");
}

/**
 * @brief Check adiabatic index constants and power functions
 */
int main(int argc, char* argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  message("Testing for gamma=%f", hydro_gamma);

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* check the values of the adiabatic index constants */
  check_constants();

  for (int i = 0; i < 100; ++i) {

    const float x = random_uniform(0., 100.);

    message("Random test %d/100 (x=%e)", i, x);

    /* check the adiabatic index power functions */
    check_functions(x);
  }

  return 0;
}
