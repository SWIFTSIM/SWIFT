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
#include <config.h>

/* Includes. */
#include "inline.h"

/* Standard headers. */
#include <math.h>

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
 * @brief Keep a wrapped position out of the cell index cdim.
 *
 * box_wrap() can return dim minus one ulp, which times the inverse cell width
 * rounds to cdim. A position at dim goes to 0. A position below dim that
 * scales to cdim goes to 0 in a periodic box, or down one ulp at a time onto
 * the last cell otherwise. A position above dim is left alone, so that the
 * debug checks still see it.
 *
 * @param pos The wrapped position.
 * @param dim The box size along this axis.
 * @param ih The inverse top-level cell width along this axis.
 * @param cdim The number of top-level cells along this axis.
 * @param periodic Whether the box is periodic.
 */
__attribute__((always_inline, const)) INLINE static double box_wrap_cell_edge(
    double pos, const double dim, const double ih, const int cdim,
    const int periodic) {
  if (pos == dim) {
    pos = 0.0;
  } else if (pos < dim && pos * ih >= cdim) {
    if (periodic)
      pos = 0.0;
    else
      while (pos * ih >= cdim) pos = nextafter(pos, 0.0);
  }
  return pos;
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
