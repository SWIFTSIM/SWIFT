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

#include "adiabatic_index.h"
#include "error.h"

/**
 * @brief Check that a and b are consistent (up to some absolute error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 */
void check_value(float a, float b, const char* s) {
  if (fabsf(a - b) > 1.e-5f) {
    error("Values are inconsistent: %g %g (%s)!", a, b, s);
  } else {
    message("Values are consistent: %g %g (%s).", a, b, s);
  }
}

/**
 * @brief Check that the pre-defined adiabatic index constants contain correct
 * values
 */
void check_constants() {
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
void check_functions() {
  float val_a, val_b;
  const float x = 0.4;

  val_a = pow(x, -hydro_gamma);
  val_b = pow_minus_gamma(x);
  check_value(val_a, val_b, "x^(-gamma)");

  val_a = pow(x, 2.0f / (hydro_gamma - 1.0f));
  val_b = pow_two_over_gamma_minus_one(x);
  check_value(val_a, val_b, "x^(2/(gamma-1))");

  val_a = pow(x, 2.0f * hydro_gamma / (hydro_gamma - 1.0f));
  val_b = pow_two_gamma_over_gamma_minus_one(x);
  check_value(val_a, val_b, "x^((2 gamma)/(gamma-1))");

  val_a = pow(x, 0.5f * (hydro_gamma - 1.0f) / hydro_gamma);
  val_b = pow_gamma_minus_one_over_two_gamma(x);
  check_value(val_a, val_b, "x^((gamma-1)/(2 gamma))");

  val_a = pow(x, -0.5f * (hydro_gamma + 1.0f) / hydro_gamma);
  val_b = pow_minus_gamma_plus_one_over_two_gamma(x);
  check_value(val_a, val_b, "x^(-(gamma+1)/(2 gamma))");

  val_a = pow(x, 1.0f / hydro_gamma);
  val_b = pow_one_over_gamma(x);
  check_value(val_a, val_b, "x^(1/gamma)");
}

/**
 * @brief Check adiabatic index constants and power functions
 */
int main() {

  /* check the values of the adiabatic index constants */
  check_constants();

  /* check the adiabatic index power functions */
  check_functions();

  return 0;
}
