/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_OPTIMIZED_EXP_H
#define SWIFT_OPTIMIZED_EXP_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "inline.h"

/**
 * @brief Compute the exponential of a number.
 *
 * This function has a relative accuracy of 1.618e-6 over the input
 * range [-32., 32.].
 *
 * @param x The number to take the exponential of.
 */
__attribute__((always_inline, const)) INLINE static float optimized_expf(
    const float x) {

  /* Let's first express e^x as 2^i * e^f with
   * f in the range [-ln(2)/2, ln(2)/2] */
  const float i = rintf(x * ((float)M_LOG2E));
  const float f = x - ((float)M_LN2) * i;

  /* We can now compute exp(f) using a polynomial
   * approximation valid over the range [-ln(2)/2, ln(2)/2].
   * The coefficients come from the Cephes library and
   * have been obtained using a minmax algorithm */
  float exp_f = 0.041944388f;
  exp_f = exp_f * f + 0.168006673f;
  exp_f = exp_f * f + 0.499999940f;
  exp_f = exp_f * f + 0.999956906f;
  exp_f = exp_f * f + 0.999999642f;

  union {
    int i;
    float f;
  } e;

  /* We can now construct the result by taking exp_f
   * as the mantissa of the answer and bit-shifting i
   * into the exponent part of the floating-point
   * number */
  e.f = exp_f;
  e.i += ((int)i) << 23;

  return e.f;
}

#endif /* SWIFT_OPTIMIZED_EXP_H */
