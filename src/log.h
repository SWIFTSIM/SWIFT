/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Willem Elbers (whe@willemelbers.com)
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
#ifndef SWIFT_OPTIMIZED_LOG_H
#define SWIFT_OPTIMIZED_LOG_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "inline.h"

/* Standard headers */
#include <math.h>
#include <stdint.h>

/**
 * @brief Compute the natural logarithm of a number.
 *
 * This function has a maximum absolute error of 5e-3 over the domain
 * [exp(-32.), exp(32.)] and a maximum relative error of 5e-3 over the two
 * intervals [exp(-32.), exp(-0.3)] and [exp(0.3), exp(32.)]
 *
 * @param x The argument of the log
 */
__attribute__((always_inline, const)) INLINE static float optimized_logf(
    float val) {
  union {
    int32_t i;
    float f;
  } e;
  e.f = val;

  /* Isolate the exponent */
  float log_2 = (float)(((e.i >> 23) & 255) - 128);
  e.i &= ~(255 << 23);
  e.i += 127 << 23;
  /* Approximation based on https://stackoverflow.com/a/28730362 comment */
  log_2 += ((-0.34484843f) * e.f + 2.02466578f) * e.f - 0.67487759f;
  /* Return the natural log */
  return log_2 * 0.693147181f;
}

/**
 * @brief Compute the base-10 logarithm of a number.
 *
 * This function has a maximum absolute error of 5e-3 over the domain
 * [exp(-32.), exp(32.)] and a maximum relative error of 5e-3 over the two
 * intervals [exp(-32.), exp(-0.3)] and [exp(0.3), exp(32.)]
 *
 * @param x The argument of the log
 */
__attribute__((always_inline, const)) INLINE static float optimized_log10f(
    float val) {
  /* Compute the natural log */
  float log_e = optimized_logf(val);
  /* Return the base-10 log */
  return log_e * 0.434294482f;
}

#endif /* SWIFT_OPTIMIZED_LOG_H */
