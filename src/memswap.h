/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_MEMSWAP_H
#define SWIFT_MEMSWAP_H

/* Config parameters. */
#include "../config.h"

#ifdef HAVE_IMMINTRIN_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <immintrin.h>
#endif

#ifdef HAVE_ALTIVEC_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <altivec.h>
#endif

/* Macro for in-place swap of two values a and b of type t. */
#define swap_loop(t, a, b, c)  \
  while (c >= sizeof(t)) {     \
    register t temp = *(t *)a; \
    *(t *)a = *(t *)b;         \
    *(t *)b = temp;            \
    a += sizeof(t);            \
    b += sizeof(t);            \
    bytes -= sizeof(t);        \
  }

/**
 * @brief Swap the contents of two elements in-place.
 *
 * Keep in mind that this function works best when the underlying data
 * is aligned to the vector length, e.g. with the @c __attribute__((aligned(32)))
 * syntax, and the code is compiled with @c -funroll-loops.
 *
 * @param a Pointer to the first element.
 * @param b Pointer to the second element.
 * @param bytes Size, in bytes, of the data pointed to by @c a and @c b.
 */
__attribute__((always_inline)) inline void memswap(void *void_a, void *void_b,
                                                   size_t bytes) {
  char *a = (char *)void_a, *b = (char *)void_b;
#ifdef __AVX512F__
  swap_loop(__m512i, a, b, bytes)
#endif
#ifdef __AVX__
  swap_loop(__m256i, a, b, bytes)
#endif
#ifdef __SSE2__
  swap_loop(__m128i, a, b, bytes)
#endif
#ifdef __ALTIVEC__
  swap_loop(vector int, a, b, bytes)
#endif
  swap_loop(size_t, a, b, bytes)
  swap_loop(int, a, b, bytes)
  swap_loop(short, a, b, bytes)
  swap_loop(char, a, b, bytes)
}

#endif /* SWIFT_MEMSWAP_H */
