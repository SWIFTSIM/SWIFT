//
// Created by yuyttenh on 24/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_UTILS_H
#define SWIFTSIM_SHADOWSWIFT_UTILS_H

#include "shadowswift/delaunay.h"

#include <float.h>
#include <math.h>
#ifdef DELAUNAY_3D_HAND_VEC
#include <immintrin.h>
#endif

static inline int sgn(double x) {
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

/**
 * @brief Check whether two doubles are equal up to the given precision.
 *
 * @param double1
 * @param double2
 * @param precision
 * @return 1 for equality 0 else.
 */
inline static int double_cmp(double double1, double double2,
                             unsigned long precision) {
  long long1, long2;
  if (double1 > 0)
    long1 = (long)(double1 * precision + .5);
  else
    long1 = (long)(double1 * precision - .5);
  if (double2 > 0)
    long2 = (long)(double2 * precision + .5);
  else
    long2 = (long)(double2 * precision - .5);
  return (long1 == long2);
}

inline static int approx_equals(double a, double b, double epsilon) {
  if (a == b) return 1;

  double abs_a = fabs(a);
  double abs_b = fabs(b);
  double diff = fabs(a - b);

  if (a == 0. || b == 0. || abs_a + abs_b < DBL_EPSILON) {
    return diff < epsilon * DBL_EPSILON;
  } else {
    return diff / fmin(abs_a + abs_b, DBL_MAX) < epsilon;
  }
}

#ifdef DELAUNAY_3D_HAND_VEC
/*! @brief Compute the absolute values of a simd vector of 4 doubles. */
inline __m256d mm256_abs_pd(__m256d x) {
  static const __m256d sign_mask = {-0., -0., -0., -0.};  // -0. = 1 << 63
  return _mm256_andnot_pd(sign_mask, x);                  // !sign_mask & x
}

/*! @brief Compute the signs (-1, 0, 1) for a simd vector of 4 doubles. */
inline __m256d mm256_sgn_pd(__m256d x) {
  static const __m256d sign_mask = {-0., -0., -0., -0.};  // -0. = 1 << 63
  static const __m256d one_mask = {1., 1., 1., 1.};
  // one_mask | (sign_mask & x): Should be -1 for negative numbers and 1 for
  // positive numbers.
  return _mm256_or_pd(one_mask, _mm256_and_pd(sign_mask, x));
}

/*! @brief Flip the signs of the first 2 elements of a simd vector of 4 doubles.
 **/
inline static __m256d mm256_flip_sign_1100(__m256d x) {
  static const __m256d mask = {-0., -0., 0., 0.};
  return _mm256_xor_pd(x, mask);
}

/**
 * @brief Compute a four by four determinant using simd instructions
 *
 * @param c0, c1, c2, c3 The columns of the matrix
 * @param c_out (Return) the vector of cofactors when developing around c3
 * @param bound (Return) The errorbound for the determinant
 * @return The determinant.
 */
inline static double determinant_4x4_bounded_simd(__m256d c0, __m256d c1,
                                                  __m256d c2, __m256d c3,
                                                  __m256d* c_out,
                                                  double* bound) {
  /* STEP 1 Compute terms with 2 factors */
  /* Contains: ax * by, bx * ay, cx * dy, dx * cy */
  __m256d ab_cd =
      _mm256_mul_pd(c0, _mm256_permute4x64_pd(c1, _MM_SHUFFLE(2, 3, 0, 1)));

  /* Contains: bx * cy, cx * by, dx * ay, ax * dy */
  __m256d bc_da =
      _mm256_mul_pd(_mm256_permute4x64_pd(c0, _MM_SHUFFLE(0, 3, 2, 1)),
                    _mm256_permute4x64_pd(c1, _MM_SHUFFLE(3, 0, 1, 2)));

  /* Contains: ax * cy, cx * ay, bx * dy, dx * by */
  __m256d ac_bd =
      _mm256_mul_pd(_mm256_permute4x64_pd(c0, _MM_SHUFFLE(3, 1, 2, 0)),
                    _mm256_permute4x64_pd(c1, _MM_SHUFFLE(1, 3, 0, 2)));

  /* STEP 2: Compute cofactors */
  /* These contain the corresponding variables from the `geometry3d_in_sphere`
   * function */
  __m256d ab_bc_cd_da = _mm256_hsub_pd(ab_cd, bc_da);
  __m256d ac_ac_bd_bd = _mm256_hsub_pd(ac_bd, ac_bd);

  /* Now compute vector of cofactors:
   * This corresponds to [abc, bcd, cda, dab] for the normal algorithm */
  __m256d cofactors = _mm256_add_pd(
      _mm256_mul_pd(
          c2, _mm256_permute4x64_pd(ab_bc_cd_da, _MM_SHUFFLE(0, 3, 2, 1))),
      _mm256_add_pd(
          _mm256_mul_pd(
              mm256_flip_sign_1100(
                  _mm256_permute4x64_pd(c2, _MM_SHUFFLE(0, 3, 2, 1))),
              _mm256_permute4x64_pd(ac_ac_bd_bd, _MM_SHUFFLE(3, 1, 2, 0))),
          _mm256_mul_pd(_mm256_permute4x64_pd(c2, _MM_SHUFFLE(1, 0, 3, 2)),
                        ab_bc_cd_da)));
  *c_out = cofactors;

  /* STEP 3: Compute the determinant */
  __m256d mul = _mm256_mul_pd(
      _mm256_permute4x64_pd(c3, _MM_SHUFFLE(2, 1, 0, 3)), cofactors);
  mul = _mm256_hsub_pd(mul, mul);
  double det = mul[0] + mul[2];

  /* STEP 4: Compute the errorbound */
  /* Now repeat step 2 and 3 to compute the errorbound */
  /* First take the absolute values */
  ab_cd = mm256_abs_pd(ab_cd);
  bc_da = mm256_abs_pd(bc_da);
  ac_bd = mm256_abs_pd(ac_bd);
  c2 = mm256_abs_pd(c2);
  /* c3 is guaranteed to be positive for our usecase */
  c3 = mm256_abs_pd(c3);
#ifdef SWIFT_DEBUG_CHECKS
  __m256d abs_c3 = mm256_abs_pd(c3);
  for (int i = 0; i < 4; i++) assert(c3[i] == abs_c3[i]);
#endif
  /* Now recompute the cofactors */
  ab_bc_cd_da = _mm256_hadd_pd(ab_cd, bc_da);
  ac_ac_bd_bd = _mm256_hadd_pd(ac_bd, ac_bd);
  cofactors = _mm256_add_pd(
      _mm256_mul_pd(
          c2, _mm256_permute4x64_pd(ab_bc_cd_da, _MM_SHUFFLE(0, 3, 2, 1))),
      _mm256_add_pd(
          _mm256_mul_pd(
              _mm256_permute4x64_pd(c2, _MM_SHUFFLE(0, 3, 2, 1)),
              _mm256_permute4x64_pd(ac_ac_bd_bd, _MM_SHUFFLE(3, 1, 2, 0))),
          _mm256_mul_pd(_mm256_permute4x64_pd(c2, _MM_SHUFFLE(1, 0, 3, 2)),
                        ab_bc_cd_da)));
  /* And finally the errorbound */
  mul = _mm256_mul_pd(_mm256_permute4x64_pd(c3, _MM_SHUFFLE(2, 1, 0, 3)),
                      cofactors);
  mul = _mm256_hadd_pd(mul, mul);
  *bound = 5. * DBL_EPSILON * (mul[0] + mul[2]);

  return det;
}
#endif

#endif  // SWIFTSIM_SHADOWSWIFT_UTILS_H
