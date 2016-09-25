/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_CBRT_H
#define SWIFT_CBRT_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <math.h>

/**
 * @brief Compute the inverse cube root of a single-precision floating-point
 * number.
 *
 * @param x_in The input value.
 *
 * @return The inverse cubic root of @c x_in. Note that this function
 *         does not care about non-finite inputs.
 */
__attribute__((always_inline)) inline float icbrtf(float x_in) {

  // Extract the exponent.
  const unsigned int x_as_uint = *((char *)&x_in);
  const int exponent = ((x_as_uint & 0x7f800000) >> 23) - 127;

  // Clear the exponent and sign to get the mantissa.
  const unsigned int x_norm_as_int = (x_as_uint & ~0xff800000) | 0x3f800000;
  const float x_norm = *((char *)&x_norm_as_int);

  // Multiply by sqrt(1/2) and subtract one, should then be in the
  // range [sqrt(1/2) - 1, sqrt(2) - 1).
  const float x = x_norm * (float)M_SQRT1_2 - 1.0f;

  // Compute the polynomial interpolant.
  float res =
      9.99976591940035e-01f +
      x * (-3.32901212909283e-01f +
           x * (2.24361110929912e-01f +
                x * (-1.88913279594895e-01f + x * 1.28384036492344e-01f)));

  // Compute the new exponent and the correction factor.
  int exponent_new = exponent;
  if (exponent_new < 0) exponent_new -= 2;
  exponent_new = -exponent_new / 3;
  const int exponent_rem = exponent + 3 * exponent_new;
  const unsigned int exponent_scale_as_int = (exponent_new + 127) << 23;
  float exponent_scale = *((char *)&exponent_scale_as_int);
  exponent_scale *=
      exponent_rem > 0
          ? (exponent_rem > 1 ? 5.61231024154687e-01f : 7.07106781186548e-01f)
          : 8.90898718140339e-01f;
  res *= exponent_scale;

  // One step of Newton iteration to refine the result.
  res *= (1.0f / 3.0f) * (4.0f - x_in * res * res * res);

  // We're done.
  return copysignf(res, x_in);
}

#endif /* SWIFT_CBRT_H */
