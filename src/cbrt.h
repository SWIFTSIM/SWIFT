/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Local headers. */
#include "inline.h"

/**
 * @brief Compute the inverse cube root of a single-precision floating-point
 * number.
 *
 * This function does not care about non-finite inputs.
 *
 * @warning This function is faster than both gcc and Intel's `cbrtf()`
 * functions on x86 systems. However, Other compilers or other architectures
 * may have faster implementations of the standard function `cbrtf()` that
 * will potentionally outperform this function.
 *
 * @param x_in The input value.
 *
 * @return The inverse cubic root of @c x_in (i.e. \f$x_{in}^{-1/3} \f$) .
 */
__attribute__((always_inline)) INLINE static float icbrtf(float x_in) {

  union {
    float as_float;
    unsigned int as_uint;
    int as_int;
  } cast;

  /* Extract the exponent. */
  cast.as_float = x_in;
  const int exponent = ((cast.as_int & 0x7f800000) >> 23) - 127;

  /* Clear the exponent and sign to get the mantissa. */
  cast.as_uint = (cast.as_uint & ~0xff800000) | 0x3f800000;
  const float x_norm = cast.as_float;

  /* Multiply by sqrt(1/2) and subtract one, should then be in the
     range [sqrt(1/2) - 1, sqrt(2) - 1). */
  const float x = x_norm * (float)M_SQRT1_2 - 1.0f;

  /* Compute the polynomial interpolant. */
  float res =
      9.99976591940035e-01f +
      x * (-3.32901212909283e-01f +
           x * (2.24361110929912e-01f +
                x * (-1.88913279594895e-01f + x * 1.28384036492344e-01f)));

  /* Compute the new exponent and the correction factor. */
  int exponent_new = exponent;
  if (exponent_new < 0) exponent_new -= 2;
  exponent_new = -exponent_new / 3;
  const int exponent_rem = exponent + 3 * exponent_new;
  cast.as_uint = (exponent_new + 127) << 23;
  static const float scale[3] = {8.90898718140339e-01f, 7.07106781186548e-01f,
                                 5.61231024154687e-01f};
  const float exponent_scale = cast.as_float * scale[exponent_rem];

  /* Scale the result and set the correct sign. */
  res = copysignf(res * exponent_scale, x_in);

  /* One step of Newton iteration to refine the result. */
  res *= (1.0f / 3.0f) * (4.0f - x_in * res * res * res);

  /* We're done. */
  return res;
}

#endif /* SWIFT_CBRT_H */
