/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_PERIODIC_H
#define SWIFT_PERIODIC_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "inline.h"

/**
 * @brief Limits the value of x to be between a and b
 *
 * Only wraps once. If x > a + 2(b - a), the returned value will be larger than
 * b. Similarly for x < a - (b - a), the value will be smaller than a.
 */
#define box_wrap(x, a, b)                                              \
  ({                                                                   \
    const __typeof__(x) _x = (x);                                      \
    const __typeof__(a) _a = (a);                                      \
    const __typeof__(b) _b = (b);                                      \
    _x < _a ? (_x + (_b - _a)) : ((_x >= _b) ? (_x - (_b - _a)) : _x); \
  })

/**
 * @brief Limits the value of x to be between a and b
 */
__attribute__((always_inline, const)) INLINE static double box_wrap_multiple(
    double x, const double a, const double b) {
  while (x < a) {
    x += (b - a);
  }
  while (x >= b) {
    x -= (b - a);
  }
  return x;
}

/**
 * @brief Find the smallest distance dx along one axis within a box of size
 * box_size
 *
 * This macro evaluates its arguments exactly once.
 *
 * Only wraps once. If dx > 2b, the returned value will be larger than b.
 * Similarly for dx < -b.
 *
 */
__attribute__((always_inline, const)) INLINE static double nearest(
    const double dx, const double box_size) {

  return ((dx > 0.5 * box_size)
              ? (dx - box_size)
              : ((dx < -0.5 * box_size) ? (dx + box_size) : dx));
}

/**
 * @brief Find the smallest distance dx along one axis within a box of size
 * box_size
 *
 * This macro evaluates its arguments exactly once.
 *
 * Only wraps once. If dx > 2b, the returned value will be larger than b.
 * Similarly for dx < -b.
 *
 */
__attribute__((always_inline, const)) INLINE static float nearestf(
    const float dx, const float box_size) {

  return ((dx > 0.5f * box_size)
              ? (dx - box_size)
              : ((dx < -0.5f * box_size) ? (dx + box_size) : dx));
}

#endif /* SWIFT_PERIODIC_H */
