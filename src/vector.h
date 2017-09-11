/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_VECTOR_H
#define SWIFT_VECTOR_H

/* Have I already read this file? */
#ifndef VEC_MACRO

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "inline.h"

#ifdef WITH_VECTORIZATION

/* Need to check whether compiler supports this (IBM does not)
   This will prevent the macros to be defined and switch off
   explicit vectorization if the compiled does not support it */
#ifdef HAVE_IMMINTRIN_H
/* Include the header file with the intrinsics for Intel architecture. */
#include <immintrin.h>
#endif

/* Define the vector macro. */
#define VEC_MACRO(elcount, type) \
  __attribute__((vector_size((elcount) * sizeof(type)))) type

/* So what will the vector size be? */

/* AVX-512 intrinsics*/
#ifdef HAVE_AVX512_F
#define VEC_HAVE_GATHER
#define VEC_SIZE 16
#define VEC_FLOAT __m512
#define VEC_DBL __m512d
#define VEC_INT __m512i
#define KNL_MASK_16 __mmask16
#define vec_load(a) _mm512_load_ps(a)
#define vec_store(a, addr) _mm512_store_ps(addr, a)
#define vec_setzero() _mm512_setzero_ps()
#define vec_setintzero() _mm512_setzero_epi32()
#define vec_set1(a) _mm512_set1_ps(a)
#define vec_setint1(a) _mm512_set1_epi32(a)
#define vec_set(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) \
  _mm512_set_ps(p, o, n, m, l, k, j, i, h, g, f, e, d, c, b, a)
#define vec_dbl_set(a, b, c, d, e, f, g, h) \
  _mm512_set_pd(h, g, f, e, d, c, b, a)
#define vec_add(a, b) _mm512_add_ps(a, b)
#define vec_mask_add(a, b, mask) _mm512_mask_add_ps(a, mask, b, a)
#define vec_sub(a, b) _mm512_sub_ps(a, b)
#define vec_mask_sub(a, b, mask) _mm512_mask_sub_ps(a, mask, a, b)
#define vec_mul(a, b) _mm512_mul_ps(a, b)
#define vec_div(a, b) _mm512_div_ps(a, b)
#define vec_fma(a, b, c) _mm512_fmadd_ps(a, b, c)
#define vec_fnma(a, b, c) _mm512_fnmadd_ps(a, b, c)
#define vec_sqrt(a) _mm512_sqrt_ps(a)
#define vec_rcp(a) _mm512_rcp14_ps(a)
#define vec_rsqrt(a) _mm512_rsqrt14_ps(a)
#define vec_ftoi(a) _mm512_cvttps_epi32(a)
#define vec_fmin(a, b) _mm512_min_ps(a, b)
#define vec_fmax(a, b) _mm512_max_ps(a, b)
#define vec_fabs(a) _mm512_andnot_ps(_mm512_set1_ps(-0.f), a)
#define vec_floor(a) _mm512_floor_ps(a)
#define vec_cmp_gt(a, b) _mm512_cmp_ps_mask(a, b, _CMP_GT_OQ)
#define vec_cmp_lt(a, b) _mm512_cmp_ps_mask(a, b, _CMP_LT_OQ)
#define vec_cmp_lte(a, b) _mm512_cmp_ps_mask(a, b, _CMP_LE_OQ)
#define vec_cmp_gte(a, b) _mm512_cmp_ps_mask(a, b, _CMP_GE_OQ)
#define vec_cmp_result(a) ({ a; })
#define vec_form_int_mask(a) ({ a; })
#define vec_and(a, b) _mm512_and_ps(a, b)
#define vec_mask_and(a, b) _mm512_kand(a, b)
#define vec_and_mask(a, mask) _mm512_maskz_mov_ps(mask, a)
#define vec_init_mask_true(mask) ({ mask = 0xFFFF; })
#define vec_zero_mask(mask) ({ mask = 0; })
#define vec_create_mask(mask, cond) ({ mask = cond; })
#define vec_combine_masks(mask1, mask2) \
  ({ mask1 = vec_mask_and(mask1, mask2); })
#define vec_pad_mask(mask, pad) ({ mask = mask >> (pad); })
#define vec_blend(mask, a, b) _mm512_mask_blend_ps(mask, a, b)
#define vec_todbl_lo(a) _mm512_cvtps_pd(_mm512_extract128_ps(a, 0))
#define vec_todbl_hi(a) _mm512_cvtps_pd(_mm512_extract128_ps(a, 1))
#define vec_dbl_tofloat(a, b) _mm512_insertf128(_mm512_castps128_ps512(a), b, 1)
#define vec_dbl_load(a) _mm512_load_pd(a)
#define vec_dbl_set1(a) _mm512_set1_pd(a)
#define vec_dbl_sqrt(a) _mm512_sqrt_pd(a)
#define vec_dbl_rcp(a) _mm512_rcp_pd(a)
#define vec_dbl_rsqrt(a) _mm512_rsqrt_pd(a)
#define vec_dbl_ftoi(a) _mm512_cvttpd_epi32(a)
#define vec_dbl_fmin(a, b) _mm512_min_pd(a, b)
#define vec_dbl_fmax(a, b) _mm512_max_pd(a, b)
#define vec_getoffsets(ptrs)                                                \
  _mm512_insertf64x4(                                                       \
      _mm512_insertf64x4(_mm512_setzero_pd(),                               \
                         _mm512_cvtepi64_epi32(_mm512_load_epi64(ptrs) -    \
                                               _mm512_set1_epi64(ptrs[0])), \
                         0),                                                \
      _mm512_cvtepi64_epi32(_mm512_load_epi64(&ptrs[4]) -                   \
                            _mm512_set1_epi64(ptrs[0])),                    \
      1)
#define vec_gather(base, offsets) _mm512_i32gather_ps(offsets.m, base, 1)

/* Initialises a vector struct with a default value. */
#define FILL_VEC(a)                                                     \
  {                                                                     \
    .f[0] = a, .f[1] = a, .f[2] = a, .f[3] = a, .f[4] = a, .f[5] = a,   \
    .f[6] = a, .f[7] = a, .f[8] = a, .f[9] = a, .f[10] = a, .f[11] = a, \
    .f[12] = a, .f[13] = a, .f[14] = a, .f[15] = a                      \
  }

/* Performs a horizontal add on the vector and adds the result to a float. */
#ifdef __ICC
#define VEC_HADD(a, b) b += _mm512_reduce_add_ps(a.v)
#else /* _mm512_reduce_add_ps not present in GCC compiler. \
       TODO: Implement intrinsic version.*/
#define VEC_HADD(a, b)                              \
  {                                                 \
    for (int i = 0; i < VEC_SIZE; i++) b += a.f[i]; \
  }
#endif

/* Do nothing in the case of AVX-512 as there are already
 * instructions for left-packing.*/
#define VEC_FORM_PACKED_MASK(mask, packed_mask) packed_mask = mask

/* Finds the horizontal maximum of vector b and returns a float. */
#define VEC_HMAX(a, b) b = _mm512_reduce_max_ps(a.v)

/* Performs a left-pack on a vector based upon a mask and returns the result. */
#define VEC_LEFT_PACK(a, mask, result) \
  _mm512_mask_compressstoreu_ps(result, mask, a)

/* AVX intrinsics */
#elif defined(HAVE_AVX)
#define VEC_SIZE 8
#define VEC_FLOAT __m256
#define VEC_DBL __m256d
#define VEC_INT __m256i
#define vec_load(a) _mm256_load_ps(a)
#define vec_unaligned_load(a) _mm256_loadu_ps(a)
#define vec_store(a, addr) _mm256_store_ps(addr, a)
#define vec_unaligned_store(a, addr) _mm256_storeu_ps(addr, a)
#define vec_setzero() _mm256_setzero_ps()
#define vec_setintzero() _mm256_setzero_si256()
#define vec_set1(a) _mm256_set1_ps(a)
#define vec_setint1(a) _mm256_set1_epi32(a)
#define vec_set(a, b, c, d, e, f, g, h) _mm256_set_ps(h, g, f, e, d, c, b, a)
#define vec_dbl_set(a, b, c, d) _mm256_set_pd(d, c, b, a)
#define vec_add(a, b) _mm256_add_ps(a, b)
#define vec_mask_add(a, b, mask) vec_add(a, vec_and(b, mask.v))
#define vec_sub(a, b) _mm256_sub_ps(a, b)
#define vec_mask_sub(a, b, mask) vec_sub(a, vec_and(b, mask.v))
#define vec_mul(a, b) _mm256_mul_ps(a, b)
#define vec_div(a, b) _mm256_div_ps(a, b)
#define vec_sqrt(a) _mm256_sqrt_ps(a)
#define vec_rcp(a) _mm256_rcp_ps(a)
#define vec_rsqrt(a) _mm256_rsqrt_ps(a)
#define vec_ftoi(a) _mm256_cvttps_epi32(a)
#define vec_fmin(a, b) _mm256_min_ps(a, b)
#define vec_fmax(a, b) _mm256_max_ps(a, b)
#define vec_fabs(a) _mm256_andnot_ps(_mm256_set1_ps(-0.f), a)
#define vec_floor(a) _mm256_floor_ps(a)
#define vec_cmp_lt(a, b) _mm256_cmp_ps(a, b, _CMP_LT_OQ)
#define vec_cmp_gt(a, b) _mm256_cmp_ps(a, b, _CMP_GT_OQ)
#define vec_cmp_lte(a, b) _mm256_cmp_ps(a, b, _CMP_LE_OQ)
#define vec_cmp_gte(a, b) _mm256_cmp_ps(a, b, _CMP_GE_OQ)
#define vec_cmp_result(a) _mm256_movemask_ps(a)
#define vec_form_int_mask(a) _mm256_movemask_ps(a.v)
#define vec_and(a, b) _mm256_and_ps(a, b)
#define vec_mask_and(a, b) _mm256_and_ps(a.v, b.v)
#define vec_and_mask(a, mask) _mm256_and_ps(a, mask.v)
#define vec_init_mask_true(mask) mask.m = vec_setint1(0xFFFFFFFF)
#define vec_create_mask(mask, cond) mask.v = cond
#define vec_combine_masks(mask1, mask2) \
  ({ mask1.v = vec_mask_and(mask1, mask2); })
#define vec_zero_mask(mask) mask.v = vec_setzero()
#define vec_pad_mask(mask, pad) \
  for (int i = VEC_SIZE - (pad); i < VEC_SIZE; i++) mask.i[i] = 0
#define vec_blend(mask, a, b) _mm256_blendv_ps(a, b, mask.v)
#define vec_todbl_lo(a) _mm256_cvtps_pd(_mm256_extract128_ps(a, 0))
#define vec_todbl_hi(a) _mm256_cvtps_pd(_mm256_extract128_ps(a, 1))
#define vec_dbl_tofloat(a, b) _mm256_insertf128(_mm256_castps128_ps256(a), b, 1)
#define vec_dbl_load(a) _mm256_load_pd(a)
#define vec_dbl_set1(a) _mm256_set1_pd(a)
#define vec_dbl_sqrt(a) _mm256_sqrt_pd(a)
#define vec_dbl_rcp(a) _mm256_rcp_pd(a)
#define vec_dbl_rsqrt(a) _mm256_rsqrt_pd(a)
#define vec_dbl_ftoi(a) _mm256_cvttpd_epi32(a)
#define vec_dbl_fmin(a, b) _mm256_min_pd(a, b)
#define vec_dbl_fmax(a, b) _mm256_max_pd(a, b)

/* Initialises a vector struct with a default value. */
#define FILL_VEC(a)                                                   \
  {                                                                   \
    .f[0] = a, .f[1] = a, .f[2] = a, .f[3] = a, .f[4] = a, .f[5] = a, \
    .f[6] = a, .f[7] = a                                              \
  }

/* Performs a horizontal add on the vector and adds the result to a float. */
#define VEC_HADD(a, b)            \
  a.v = _mm256_hadd_ps(a.v, a.v); \
  a.v = _mm256_hadd_ps(a.v, a.v); \
  b += a.f[0] + a.f[4];

/* Performs a horizontal maximum on the vector and takes the maximum of the
 * result with a float, b. */
#define VEC_HMAX(a, b)                                     \
  {                                                        \
    for (int k = 0; k < VEC_SIZE; k++) b = max(b, a.f[k]); \
  }

/* Returns the lower 128-bits of the 256-bit vector. */
#define VEC_GET_LOW(a) _mm256_castps256_ps128(a)

/* Returns the higher 128-bits of the 256-bit vector. */
#define VEC_GET_HIGH(a) _mm256_extractf128_ps(a, 1)

/* Check if we have AVX2 intrinsics alongside AVX */
#ifdef HAVE_AVX2
#define vec_fma(a, b, c) _mm256_fmadd_ps(a, b, c)
#define vec_fnma(a, b, c) _mm256_fnmadd_ps(a, b, c)

/* Used in VEC_FORM_PACKED_MASK */
#define identity_indices 0x0706050403020100
#define VEC_HAVE_GATHER
#define vec_gather(base, offsets) _mm256_i32gather_ps(base, offsets.m, 1)

/* Takes an integer mask and forms a left-packed integer vector
 * containing indices of the set bits in the integer mask.
 * Also returns the total number of bits set in the mask. */
#define VEC_FORM_PACKED_MASK(mask, packed_mask)                                \
  {                                                                            \
    unsigned long expanded_mask = _pdep_u64(mask, 0x0101010101010101);         \
    expanded_mask *= 0xFF;                                                     \
    unsigned long wanted_indices = _pext_u64(identity_indices, expanded_mask); \
    __m128i bytevec = _mm_cvtsi64_si128(wanted_indices);                       \
    packed_mask.m = _mm256_cvtepu8_epi32(bytevec);                             \
  }

/* Performs a left-pack on a vector based upon a mask and returns the result. */
#define VEC_LEFT_PACK(a, mask, result) \
  vec_unaligned_store(_mm256_permutevar8x32_ps(a, mask.m), result)
#endif /* HAVE_AVX2 */

/* Create an FMA using vec_add and vec_mul if AVX2 is not present. */
#ifndef vec_fma
#define vec_fma(a, b, c) vec_add(vec_mul(a, b), c)
#endif

/* Create a negated FMA using vec_sub and vec_mul if AVX2 is not present. */
#ifndef vec_fnma
#define vec_fnma(a, b, c) vec_sub(c, vec_mul(a, b))
#endif

/* Form a packed mask without intrinsics if AVX2 is not present. */
#ifndef VEC_FORM_PACKED_MASK

/* Takes an integer mask and forms a left-packed integer vector
 * containing indices of the set bits in the integer mask.
 * Also returns the total number of bits set in the mask. */
#define VEC_FORM_PACKED_MASK(mask, v_mask, pack)   \
  {                                                \
    for (int i = 0; i < VEC_SIZE; i++)             \
      if ((mask & (1 << i))) v_mask.i[pack++] = i; \
  }

/* Takes two integer masks and forms two left-packed integer vectors
 * containing indices of the set bits in each corresponding integer mask.
 * Also returns the total number of bits set in the mask. */
#define VEC_FORM_PACKED_MASK_2(mask, v_mask, pack, mask2, v_mask2, pack2) \
  {                                                                       \
    for (int i = 0; i < VEC_SIZE; i++) {                                  \
      if ((mask & (1 << i))) v_mask.i[pack++] = i;                        \
      if ((mask2 & (1 << i))) v_mask2.i[pack2++] = i;                     \
    }                                                                     \
  }
#endif

/* Performs a left-pack on a vector based upon a mask and returns the result. */
/* This uses AVX intrinsics, but this is slower than performing the left-pack
 * manually by looping over the vectors. */
#ifndef VEC_LEFT_PACK
#define VEC_LEFT_PACK(a, mask, result)                                     \
  {                                                                        \
    __m256 t1 = _mm256_castps128_ps256(_mm256_extractf128_ps(a, 1));       \
    __m256 t2 = _mm256_insertf128_ps(t1, _mm256_castps256_ps128(a), 1);    \
    __m256 r0 = _mm256_permutevar_ps(a, mask);                             \
    __m256 r1 = _mm256_permutevar_ps(t2, mask);                            \
    __m128i k1 = _mm_slli_epi32(                                           \
        (__m128i)(_mm_xor_si128((__m128i)VEC_GET_HIGH((__m256)mask),       \
                                (__m128i)_mm_set1_epi32(4))),              \
        29);                                                               \
    __m128i k0 = _mm_slli_epi32((__m128i)(VEC_GET_LOW((__m256)mask)), 29); \
    __m256 kk =                                                            \
        _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_castsi128_ps(k0)), \
                             _mm_castsi128_ps(k1), 1);                     \
    *((__m256 *)(result)) = _mm256_blendv_ps(r0, r1, kk);                  \
  }
#endif /* HAVE_AVX2 */

/* SSE intrinsics*/
#elif defined(HAVE_SSE2)
#define VEC_SIZE 4
#define VEC_FLOAT __m128
#define VEC_DBL __m128d
#define VEC_INT __m128i
#define vec_load(a) _mm_load_ps(a)
#define vec_store(a, addr) _mm_store_ps(addr, a)
#define vec_setzero() _mm_setzero_ps()
#define vec_setintzero() _mm_setzero_si256()
#define vec_set1(a) _mm_set1_ps(a)
#define vec_setint1(a) _mm_set1_epi32(a)
#define vec_set(a, b, c, d) _mm_set_ps(d, c, b, a)
#define vec_dbl_set(a, b) _mm_set_pd(b, a)
#define vec_add(a, b) _mm_add_ps(a, b)
#define vec_sub(a, b) _mm_sub_ps(a, b)
#define vec_mul(a, b) _mm_mul_ps(a, b)
#define vec_div(a, b) _mm_div_ps(a, b)
#define vec_sqrt(a) _mm_sqrt_ps(a)
#define vec_rcp(a) _mm_rcp_ps(a)
#define vec_rsqrt(a) _mm_rsqrt_ps(a)
#define vec_ftoi(a) _mm_cvttps_epi32(a)
#define vec_fmin(a, b) _mm_min_ps(a, b)
#define vec_fmax(a, b) _mm_max_ps(a, b)
#define vec_fabs(a) _mm_andnot_ps(_mm_set1_ps(-0.f), a)
#define vec_floor(a) _mm_floor_ps(a)
#define vec_cmp_gt(a, b) _mm_cmpgt_ps(a, b)
#define vec_cmp_lt(a, b) _mm_cmplt_ps(a, b)
#define vec_cmp_lte(a, b) _mm_cmp_ps(a, b, _CMP_LE_OQ)
#define vec_cmp_result(a) _mm_movemask_ps(a)
#define vec_and(a, b) _mm_and_ps(a, b)
#define vec_todbl_lo(a) _mm_cvtps_pd(a)
#define vec_todbl_hi(a) _mm_cvtps_pd(_mm_movehl_ps(a, a))
#define vec_dbl_tofloat(a, b) _mm_movelh_ps(_mm_cvtpd_ps(a), _mm_cvtpd_ps(b))
#define vec_dbl_load(a) _mm_load_pd(a)
#define vec_dbl_set1(a) _mm_set1_pd(a)
#define vec_dbl_sqrt(a) _mm_sqrt_pd(a)
#define vec_dbl_rcp(a) _mm_rcp_pd(a)
#define vec_dbl_rsqrt(a) _mm_rsqrt_pd(a)
#define vec_dbl_ftoi(a) _mm_cvttpd_epi32(a)
#define vec_dbl_fmin(a, b) _mm_min_pd(a, b)
#define vec_dbl_fmax(a, b) _mm_max_pd(a, b)

/* Initialises a vector struct with a default value. */
#define FILL_VEC(a) \
  { .f[0] = a, .f[1] = a, .f[2] = a, .f[3] = a }

/* Performs a horizontal add on the vector and adds the result to a float. */
#define VEC_HADD(a, b)         \
  a.v = _mm_hadd_ps(a.v, a.v); \
  b += a.f[0] + a.f[1];

/* Create an FMA using vec_add and vec_mul if AVX2 is not present. */
#ifndef vec_fma
#define vec_fma(a, b, c) vec_add(vec_mul(a, b), c)
#endif
#else
#define VEC_SIZE 4
#endif /* HAVE_SSE2 */

/* Define the composite types for element access. */
typedef union {
  VEC_FLOAT v;
  VEC_DBL vd;
  VEC_INT m;
  float f[VEC_SIZE];
  double d[VEC_SIZE / 2];
  int i[VEC_SIZE];
} vector;

/* Define the mask type depending on the instruction set used. */
#ifdef HAVE_AVX512_F
typedef __mmask16 mask_t;
#else
typedef vector mask_t;
#endif

/**
 * @brief Calculates the inverse ($1/x$) of a vector using intrinsics and a
 * Newton iteration to obtain the correct level of accuracy.
 *
 * @param x #vector to be inverted.
 * @return x_inv #vector inverted x.
 */
__attribute__((always_inline)) INLINE vector vec_reciprocal(vector x) {

  vector x_inv;

  x_inv.v = vec_rcp(x.v);
  x_inv.v = vec_sub(x_inv.v,
                    vec_mul(x_inv.v, (vec_fma(x.v, x_inv.v, vec_set1(-1.0f)))));

  return x_inv;
}

/**
 * @brief Calculates the inverse and square root (\f$1/\sqrt{x}\f$) of a vector
 * using intrinsics and a Newton iteration to obtain the correct level of
 * accuracy.
 *
 * @param x #vector to be inverted.
 * @return x_inv #vector inverted x.
 */
__attribute__((always_inline)) INLINE vector vec_reciprocal_sqrt(vector x) {

  vector x_inv;

  x_inv.v = vec_rsqrt(x.v);
  x_inv.v = vec_sub(
      x_inv.v,
      vec_mul(vec_mul(vec_set1(0.5f), x_inv.v),
              (vec_fma(x.v, vec_mul(x_inv.v, x_inv.v), vec_set1(-1.0f)))));

  return x_inv;
}

#else
/* Needed for cache alignment. */
#define VEC_SIZE 16
#endif /* WITH_VECTORIZATION */

#endif /* VEC_MACRO */

#endif /* SWIFT_VECTOR_H */
