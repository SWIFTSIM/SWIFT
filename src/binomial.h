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
#ifndef SWIFT_BINOMIAL_H
#define SWIFT_BINOMIAL_H

#include "error.h"
#include "inline.h"

/**
 * @brief Compute the binomial coefficient (n, k)
 *
 * Only valid for values up to n = 5.
 */
__attribute__((const)) INLINE static int binomial(const int n, const int k) {

#ifdef SWIFT_DEBUG_CHECKS
  assert(n >= 0);
  assert(k >= 0);
  assert(n <= 5);
  assert(k <= n);
#endif

  static const int coeffs[6][6] = {{1, 0, 0, 0, 0, 0}, {1, 1, 0, 0, 0, 0},
                                   {1, 2, 1, 0, 0, 0}, {1, 3, 3, 1, 0, 0},
                                   {1, 4, 6, 4, 1, 0}, {1, 5, 10, 10, 5, 1}};

  return coeffs[n][k];
}

#endif /* SWIFT_BINOMIAL_H */
