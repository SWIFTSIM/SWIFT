/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_INTEGER_POWER_H
#define SWIFT_INTEGER_POWER_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "error.h"
#include "inline.h"

/* Standard headers */
#include <math.h>

/**
 * @brief Computes the power of x to the n for a (small) positive integer n.
 *
 * Only optimized for values 0 <= n <= 8. Defaults to pow() above.
 */
__attribute__((const)) INLINE static double integer_pow(const double x,
                                                        const unsigned int n) {

  switch (n) {
    case 0:
      return 1.;
    case 1:
      return x;
    case 2:
      return x * x;
    case 3:
      return x * x * x;
    case 4: {
      const double y = x * x;
      return y * y;
    }
    case 5: {
      const double y = x * x;
      return x * y * y;
    }
    case 6: {
      const double y = x * x;
      return y * y * y;
    }
    case 7: {
      const double y = x * x;
      return x * y * y * y;
    }
    case 8: {
      const double y = x * x;
      const double z = y * y;
      return z * z;
    }
    default:
      return pow(x, (double)n);
  }
}

/**
 * @brief Computes the power of x to the n for a (small) positive integer n.
 *
 * Only optimized for values 0 <= n <= 8. Defaults to powf() above.
 */
__attribute__((const)) INLINE static float integer_powf(const float x,
                                                        const unsigned int n) {

  switch (n) {
    case 0:
      return 1.f;
    case 1:
      return x;
    case 2:
      return x * x;
    case 3:
      return x * x * x;
    case 4: {
      const float y = x * x;
      return y * y;
    }
    case 5: {
      const float y = x * x;
      return x * y * y;
    }
    case 6: {
      const float y = x * x;
      return y * y * y;
    }
    case 7: {
      const float y = x * x;
      return x * y * y * y;
    }
    case 8: {
      const float y = x * x;
      const float z = y * y;
      return z * z;
    }
    default:
      return powf(x, (float)n);
  }
}

#endif /* SWIFT_INTEGER_POWER_H */
