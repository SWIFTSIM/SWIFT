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

void check_value(float a, float b, const char* s) {
  if (fabsf(a - b) > 1.e-6f) {
    error("Values are inconsistent: %g %g (%s)!", a, b, s);
  } else {
    message("Values are consistent: %g %g (%s).", a, b, s);
  }
}

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

int main() {

  check_constants();

  return 0;
}
