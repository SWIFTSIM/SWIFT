/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016   Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_MINMAX_H
#define SWIFT_MINMAX_H

/**
 * @brief Minimum of two numbers
 *
 * This macro evaluates its arguments exactly once.
 */
#define min(a, b)                 \
  ({                              \
    const __typeof__(a) _a = (a); \
    const __typeof__(b) _b = (b); \
    _a < _b ? _a : _b;            \
  })

/**
 * @brief Maximum of two numbers
 *
 * This macro evaluates its arguments exactly once.
 */
#define max(a, b)                 \
  ({                              \
    const __typeof__(a) _a = (a); \
    const __typeof__(b) _b = (b); \
    _a > _b ? _a : _b;            \
  })

/**
 * @brief Minimum of three numbers
 *
 * This macro evaluates its arguments exactly once.
 */
#define min3(x, y, z)                        \
  ({                                         \
    const __typeof__(x) _x = (x);            \
    const __typeof__(y) _y = (y);            \
    const __typeof__(z) _z = (z);            \
    const __typeof__(x) _temp = min(_x, _y); \
    min(_temp, _z);                          \
  })

/**
 * @brief Maximum of three numbers
 *
 * This macro evaluates its arguments exactly once.
 */
#define max3(x, y, z)                        \
  ({                                         \
    const __typeof__(x) _x = (x);            \
    const __typeof__(y) _y = (y);            \
    const __typeof__(z) _z = (z);            \
    const __typeof__(x) _temp = max(_x, _y); \
    max(_temp, _z);                          \
  })

#endif /* SWIFT_MINMAX_H */
