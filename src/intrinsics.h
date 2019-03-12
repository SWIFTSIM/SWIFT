/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Matthieu Schaller matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_INTRINSICS_H
#define SWIFT_INTRINSICS_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "inline.h"

/**
 * @brief Returns the number of leading 0-bits in x, starting at the most
 * significant bit position. If x is 0, the result is undefined.
 *
 * This is a wrapper for the GNU intrinsic with an implementation (from
 * Hacker's Delight) if the compiler intrinsics are not available.
 */
__attribute__((always_inline, const)) INLINE static int intrinsics_clz(
    unsigned int x) {

#ifdef __GNUC__
  /* Use GCC intrinsics if possible */
  return __builtin_clz(x);
#else
  int n;

  if (x == 0) return (32);
  n = 0;
  if (x <= 0x0000FFFF) {
    n = n + 16;
    x = x << 16;
  }
  if (x <= 0x00FFFFFF) {
    n = n + 8;
    x = x << 8;
  }
  if (x <= 0x0FFFFFFF) {
    n = n + 4;
    x = x << 4;
  }
  if (x <= 0x3FFFFFFF) {
    n = n + 2;
    x = x << 2;
  }
  if (x <= 0x7FFFFFFF) {
    n = n + 1;
  }
  return n;
#endif
}

/**
 * @brief Returns the number of leading 0-bits in x, starting at the most
 * significant bit position. If x is 0, the result is undefined.
 *
 * This is a wrapper for the GNU intrinsic with a place-holder for a future
 * version in cases where the compiler intrinsic is not available.
 */
__attribute__((always_inline, const)) INLINE static int intrinsics_clzll(
    unsigned long long x) {

#ifdef __GNUC__
  /* Use GCC intrinsics if possible */
  return __builtin_clzll(x);
#else
#error "Missing definition of clz for long long on this platform."
#endif
}

/**
 * @brief Returns the number of 1-bits in x.
 *
 * This is a wrapper for the GNU intrinsic with an implementation (from
 * Hacker's Delight) if the compiler intrinsics are not available.
 */
__attribute__((always_inline, const)) INLINE static int intrinsics_popcount(
    unsigned int x) {

#ifdef __GNUC__
  /* Use GCC intrinsics if possible */
  return __builtin_popcount(x);
#else
  x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
  x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
  x = (x & 0x0000FFFF) + ((x >> 16) & 0x0000FFFF);
  return x;
#endif
}

/**
 * @brief Returns the number of 1-bits in x.
 *
 * This is a wrapper for the GNU intrinsic with an implementation (from
 * Hacker's Delight) if the compiler intrinsics are not available.
 */
__attribute__((always_inline, const)) INLINE static int intrinsics_popcountll(
    unsigned long long x) {

#ifdef __GNUC__
  /* Use GCC intrinsics if possible */
  return __builtin_popcountll(x);
#else
#error "Missing definition of popcount for long long on this platform."
#endif
}

#endif /* SWIFT_INTRINSICS_H */
