/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2018 STFC (author email aidan.chalk@stfc.ac.uk)
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
#include <stdint.h>
#include "../config.h"

#ifdef HAVE_IMMINTRIN_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <immintrin.h>
#endif

#ifdef HAVE_ALTIVEC_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <altivec.h>
#endif

/* Macro for in-place swap of two values a and b of type t. a and b are
   assumed to be of type uint8_t* so that the pointer arithmetic works. */
#define swap_loop(type, a, b, count) \
  while (count >= sizeof(type)) {    \
    register type temp = *(type *)a; \
    *(type *)a = *(type *)b;         \
    *(type *)b = temp;               \
    a += sizeof(type);               \
    b += sizeof(type);               \
    count -= sizeof(type);           \
  }

/**
 * @brief Swap the contents of two elements in-place.
 *
 * Keep in mind that this function only works when the underlying data
 * is aligned to the vector length, e.g. with the @c
 * __attribute__((aligned(32))) syntax!
 * Furthermore, register re-labeling only seems to work when the code is
 * compiled with @c -funroll-loops.
 *
 * Note that GCC (at least until 7.3) produces incorrect AVX512 code here
 * by automatically assuming alignment.
 *
 * @param void_a Pointer to the first element.
 * @param void_b Pointer to the second element.
 * @param bytes Size, in bytes, of the data pointed to by @c a and @c b.
 */
__attribute__((always_inline)) inline void memswap(void *restrict void_a,
                                                   void *restrict void_b,
                                                   size_t bytes) {
  int8_t *restrict a = (int8_t *)void_a, *restrict b = (int8_t *)void_b;
#if defined(__AVX512F__) && defined(__INTEL_COMPILER)
  swap_loop(__m512i, a, b, bytes);
#endif
#ifdef __AVX__
  swap_loop(__m256i, a, b, bytes);
#endif
#ifdef __SSE2__
  swap_loop(__m128i, a, b, bytes);
#endif
#ifdef __ALTIVEC__
  swap_loop(vector int, a, b, bytes);
#endif
  swap_loop(int_least64_t, a, b, bytes);
  swap_loop(int_least32_t, a, b, bytes);
  swap_loop(int_least16_t, a, b, bytes);
  swap_loop(int_least8_t, a, b, bytes);

  /* This is a known bug for the current version of clang on ARM.
   * We add this synchronization as a temporary bug fix.
   * See https://bugs.llvm.org/show_bug.cgi?id=40051 */
#if defined(__clang__) && defined(__aarch64__)
  __sync_synchronize();
#endif
}

/**
 * @brief Swap the contents of two elements in-place.
 *
 * As opposed to #memswap, this function does not require the parameters
 * to be aligned in any specific way.
 * Furthermore, register re-labeling only seems to work when the code is
 * compiled with @c -funroll-loops.
 *
 * @param void_a Pointer to the first element.
 * @param void_b Pointer to the second element.
 * @param bytes Size, in bytes, of the data pointed to by @c a and @c b.
 */
__attribute__((always_inline)) inline void memswap_unaligned(
    void *restrict void_a, void *restrict void_b, size_t bytes) {
  int8_t *restrict a = (int8_t *)void_a, *restrict b = (int8_t *)void_b;
#ifdef __AVX512F__
  while (bytes >= sizeof(__m512i)) {
    register __m512i temp;
    temp = _mm512_loadu_si512((__m512i *)a);
    _mm512_storeu_si512((__m512i *)a, _mm512_loadu_si512((__m512i *)b));
    _mm512_storeu_si512((__m512i *)b, temp);
    a += sizeof(__m512i);
    b += sizeof(__m512i);
    bytes -= sizeof(__m512i);
  }
#endif
#ifdef __AVX__
  while (bytes >= sizeof(__m256i)) {
    register __m256i temp;
    temp = _mm256_loadu_si256((__m256i *)a);
    _mm256_storeu_si256((__m256i *)a, _mm256_loadu_si256((__m256i *)b));
    _mm256_storeu_si256((__m256i *)b, temp);
    a += sizeof(__m256i);
    b += sizeof(__m256i);
    bytes -= sizeof(__m256i);
  }
#endif
#ifdef __SSE2__
  while (bytes >= sizeof(__m128i)) {
    register __m128i temp;
    temp = _mm_loadu_si128((__m128i *)a);
    _mm_storeu_si128((__m128i *)a, _mm_loadu_si128((__m128i *)b));
    _mm_storeu_si128((__m128i *)b, temp);
    a += sizeof(__m128i);
    b += sizeof(__m128i);
    bytes -= sizeof(__m128i);
  }
#endif
#ifdef __ALTIVEC__
  // Power8 supports unaligned load/stores, but not sure what it will do here.
  swap_loop(vector int, a, b, bytes);
#endif
  swap_loop(int_least64_t, a, b, bytes);
  swap_loop(int_least32_t, a, b, bytes);
  swap_loop(int_least16_t, a, b, bytes);
  swap_loop(int_least8_t, a, b, bytes);

  /* This is a known bug for the current version of clang on ARM.
   * We add this synchronization as a temporary bug fix.
   * See https://bugs.llvm.org/show_bug.cgi?id=40051 */
#if defined(__clang__) && defined(__aarch64__)
  __sync_synchronize();
#endif
}

#endif /* SWIFT_MEMSWAP_H */
