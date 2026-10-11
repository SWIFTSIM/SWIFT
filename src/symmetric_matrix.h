/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SYMMETRIC_MATRIX_H
#define SWIFT_SYMMETRIC_MATRIX_H

/* Local includes */
#include "dimension.h"
#include "error.h"

/* Standard headers */
#include <stdint.h>
#include <string.h>

#if defined(HYDRO_DIMENSION_3D)

#define sym_matrix_num_elements 6

#elif defined(HYDRO_DIMENSION_2D)

#define sym_matrix_num_elements 3

#elif defined(HYDRO_DIMENSION_1D)

#define sym_matrix_num_elements 1

#else
#error "A problem dimensionality must be chosen in config.h !"
#endif

/**
 * @brief Symmetric matrix definition in 3D.
 *
 * The matrix elements can be accessed as an array or via their "coordinates".
 * Replicated elements are not stored.
 */
struct sym_matrix {

  union {
    struct {
      float elements[sym_matrix_num_elements];
    };
    struct {
#if defined(HYDRO_DIMENSION_3D)
      float xx;
      float yy;
      float zz;
      float xy;
      float xz;
      float yz;
#elif defined(HYDRO_DIMENSION_2D)
      float xx;
      float yy;
      float xy;
#elif defined(HYDRO_DIMENSION_1D)
      float xx;
#endif
    };
  };
};

/**
 * @brief Zero the matrix
 */
__attribute__((always_inline)) INLINE static void zero_sym_matrix(
    struct sym_matrix *M) {
  for (int i = 0; i < sym_matrix_num_elements; ++i) M->elements[i] = 0.f;
}

/**
 * @brief Construct the identity matrix
 */
__attribute__((always_inline)) INLINE static void sym_matrix_identity(
    struct sym_matrix *M) {
  M->xx = 1.f;
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
  M->yy = 1.f;
  M->xy = 0.f;
#endif
#if defined(HYDRO_DIMENSION_3D)
  M->zz = 1.f;
  M->xz = 0.f;
  M->yz = 0.f;
#endif
}

/**
 * @brief Check whether this is the null matrix
 */
__attribute__((always_inline)) INLINE static int sym_matrix_is_null(
    const struct sym_matrix *M) {

  for (int i = 0; i < sym_matrix_num_elements; ++i) {
    if (M->elements[i] != 0.f) return 0;
  }
  return 1;
}

/**
 * @brief Construct a 3x3 array from a symmetric matrix.
 */
__attribute__((always_inline)) INLINE static void get_matrix_from_sym_matrix(
    float out[hydro_dimension_integer][hydro_dimension_integer],
    const struct sym_matrix *in) {

  out[0][0] = in->xx;
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
  out[1][1] = in->yy;
  out[0][1] = in->xy;
  out[1][0] = in->xy;
#endif
#if defined(HYDRO_DIMENSION_3D)
  out[0][2] = in->xz;
  out[1][2] = in->yz;
  out[2][0] = in->xz;
  out[2][1] = in->yz;
  out[2][2] = in->zz;
#endif
}

/**
 * @brief Construct a symmetric matrix from a 3x3 array.
 *
 * No check is performed to verify the input 3x3 array is indeed symmetric.
 * The upper half of the matrix is simply take as-is.
 */
__attribute__((always_inline)) INLINE static void get_sym_matrix_from_matrix(
    struct sym_matrix *out,
    const float in[hydro_dimension_integer][hydro_dimension_integer]) {

  out->xx = in[0][0];
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
  out->yy = in[1][1];
  out->xy = in[0][1];
#endif
#if defined(HYDRO_DIMENSION_3D)
  out->zz = in[2][2];
  out->xz = in[0][2];
  out->yz = in[1][2];
#endif
}

/**
 * @brief Compute the product of a symmetric matrix and a vector.
 *
 * Performs out = M * v.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_multiply_by_vector(
    float out[hydro_dimension_integer], const struct sym_matrix *M,
    const float v[hydro_dimension_integer]) {
#if defined(HYDRO_DIMENSION_3D)
  out[0] = M->xx * v[0] + M->xy * v[1] + M->xz * v[2];
  out[1] = M->xy * v[0] + M->yy * v[1] + M->yz * v[2];
  out[2] = M->xz * v[0] + M->yz * v[1] + M->zz * v[2];
#elif defined(HYDRO_DIMENSION_2D)
  out[0] = M->xx * v[0] + M->xy * v[1];
  out[1] = M->xy * v[0] + M->yy * v[1];
#elif defined(HYDRO_DIMENSION_1D)
  out[0] = M->xx * v[0];
#else
#error "A problem dimensionality must be chosen in config.h !"
#endif
}

/**
 * @brief Compute the product of a symmetric matrix and a scalar.
 *
 * Performs M *= a.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_multiply_by_scalar(
    struct sym_matrix *M, const float alpha) {
  for (int i = 0; i < sym_matrix_num_elements; ++i) M->elements[i] *= alpha;
}

/**
 * @brief Multiply two symmetric matrices in two operations, ABA.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_multiplication_ABA(
    struct sym_matrix *M_out, const struct sym_matrix *A,
    const struct sym_matrix *B) {
#if defined(HYDRO_DIMENSION_3D)
  float BA_array[3][3] = {0};
  float A_array[3][3], B_array[3][3];
  get_matrix_from_sym_matrix(A_array, A);
  get_matrix_from_sym_matrix(B_array, B);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        BA_array[i][j] += B_array[i][k] * A_array[k][j];
      }
    }
  }

  M_out->xx = (A->xx * BA_array[0][0] + A->xy * BA_array[1][0] +
               A->xz * BA_array[2][0]);
  M_out->yy = (A->xy * BA_array[0][1] + A->yy * BA_array[1][1] +
               A->yz * BA_array[2][1]);
  M_out->zz = (A->xz * BA_array[0][2] + A->yz * BA_array[1][2] +
               A->zz * BA_array[2][2]);
  M_out->xy = (A->xy * BA_array[0][0] + A->yy * BA_array[1][0] +
               A->yz * BA_array[2][0]);
  M_out->xz = (A->xz * BA_array[0][0] + A->yz * BA_array[1][0] +
               A->zz * BA_array[2][0]);
  M_out->yz = (A->xz * BA_array[0][1] + A->yz * BA_array[1][1] +
               A->zz * BA_array[2][1]);
#else
  error("Function only exists in 3D!");
#endif
}

/**
 * @brief Print a symmetric matrix.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_print(
    const struct sym_matrix *M) {
#if defined(HYDRO_DIMENSION_3D)
  message("|%.7f %.7f %.7f|", M->xx, M->xy, M->xz);
  message("|%.7f %.7f %.7f|", M->xy, M->yy, M->yz);
  message("|%.7f %.7f %.7f|", M->xz, M->yz, M->zz);
#elif defined(HYDRO_DIMENSION_2D)
  message("|%.7f %.7f|", M->xx, M->xy);
  message("|%.7f %.7f|", M->xy, M->yy);
#elif defined(HYDRO_DIMENSION_1D)
  message("|%.7f|", M->xx);
#else
#error "A problem dimensionality must be chosen in config.h !"
#endif
}

/**
 * @brief Check that all the elements of a symmetric matrix are finite.
 *
 * Uses the exponent bits, as isfinite()/isnan() and comparisons with NaN are
 * not reliable under -ffast-math (which SWIFT always uses).
 */
__attribute__((always_inline)) INLINE static int sym_matrix_is_finite(
    const struct sym_matrix *M) {
  for (int i = 0; i < sym_matrix_num_elements; ++i) {
    uint32_t bits;
    memcpy(&bits, &M->elements[i], sizeof(bits));
    if ((bits & 0x7f800000u) == 0x7f800000u) return 0;
  }
  return 1;
}

/**
 * @brief Compute the inverse of a symmetric matrix.
 *
 * The inversion is performed in double precision. The singularity and
 * condition number checks are independent of the overall scale of M.
 *
 * @param M_inv (return) The inverse of M.
 * @param M The symmetric matrix to invert.
 * @param max_cond_num Maximal 2-norm condition number to attempt an inversion.
 * Larger values will trigger the singular matrix case.
 * @return 1 if the inversion has failed. The matrix M_inv is then the null
 * matrix. 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int sym_matrix_invert(
    struct sym_matrix *restrict M_inv, const struct sym_matrix *restrict M,
    const double max_cond_num) {

  /* Non-finite input: fail (the checks below cannot be trusted with NaNs) */
  if (!sym_matrix_is_finite(M)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

#if defined(HYDRO_DIMENSION_3D)

  /* Turn the matrix into a (double) 3x3 array */
  float A[3][3];
  get_matrix_from_sym_matrix(A, M);

  double A_d[3][3];
  double M_inv_matrix[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A_d[i][j] = A[i][j];
    }
  }

  /* Abort if the condition number is bad */
  const double cond_number = matrix_3x3_symmetric_2norm_condition_number(A_d);
  if (!(cond_number <= max_cond_num)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  /* Invert */
  if (invert3x3_matrix_LU(A_d, M_inv_matrix)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  /* Save the resulting matrix back into a (float) sym_matrix object */
  M_inv->xx = M_inv_matrix[0][0];
  M_inv->yy = M_inv_matrix[1][1];
  M_inv->zz = M_inv_matrix[2][2];
  M_inv->xy = M_inv_matrix[0][1];
  M_inv->xz = M_inv_matrix[0][2];
  M_inv->yz = M_inv_matrix[1][2];

  /* The (float) result must be finite too */
  if (!sym_matrix_is_finite(M_inv)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  return 0;

#elif defined(HYDRO_DIMENSION_2D)

  const double a = M->xx;
  const double b = M->xy;
  const double c = M->yy;

  /* Eigenvalues of the symmetric matrix. The small one is obtained from the
   * determinant to avoid cancellation. For a symmetric matrix, the singular
   * values are the absolute values of the eigenvalues. */
  const double mean = 0.5 * (a + c);
  const double half_diff = 0.5 * (a - c);
  const double disc = sqrt(half_diff * half_diff + b * b);
  const double ev_big = mean + copysign(disc, mean);
  const double det = a * c - b * b;

  if (ev_big == 0.) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  const double ev_small = det / ev_big;

  /* Abort if the condition number is bad */
  const double cond_number = fabs(ev_small) > 1e-15 * fabs(ev_big)
                                 ? fabs(ev_big / ev_small)
                                 : INFINITY;
  if (!(cond_number <= max_cond_num)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  /* Invert */
  const double det_inv = 1. / det;
  M_inv->xx = c * det_inv;
  M_inv->yy = a * det_inv;
  M_inv->xy = -b * det_inv;

  /* The (float) result must be finite too */
  if (!sym_matrix_is_finite(M_inv)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  return 0;

#elif defined(HYDRO_DIMENSION_1D)

  /* The condition number of a non-zero 1x1 matrix is 1 */
  if (M->xx == 0.f || max_cond_num < 1.) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  M_inv->xx = 1.f / M->xx;

  /* The (float) result must be finite too */
  if (!sym_matrix_is_finite(M_inv)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  return 0;

#else
#error "A problem dimensionality must be chosen in config.h !"
#endif
}

/**
 * @brief Compute the inverse of a symmetric positive semi-definite matrix,
 * regularised in its ill-conditioned directions.
 *
 * With M = sum_k lambda_k e_k e_k^T and the floor l_f = lambda_max /
 * max_cond_num, the result is
 *   M_inv = sum_k e_k e_k^T f(lambda_k),  f(l) = 1/l for l >= l_f,
 *                                          f(l) = max(l, 0) / l_f^2 otherwise,
 * i.e. the exact inverse when the condition number is within the limit, and
 * otherwise a spectrally-filtered inverse: the well-resolved directions keep
 * their exact inverse while the (near-)degenerate ones are damped, down to
 * zero for an exactly singular direction (and for negative eigenvalues from
 * round-off). f is continuous and bounded by max_cond_num / lambda_max. The
 * intended use is the moment matrix of a neighbourhood that is planar or
 * filamentary (rank-deficient in one or two directions) but well sampled in
 * the others.
 *
 * In 1D, there is no direction to regularise: the function reduces to
 * sym_matrix_invert().
 *
 * @param M_inv (return) The (regularised) inverse of M.
 * @param M The symmetric matrix to invert.
 * @param max_cond_num Maximal 2-norm condition number of the result.
 * @param regularised (return) 1 if at least one eigenvalue was clipped.
 * @return 1 if the inversion has failed (non-finite input, no positive
 * eigenvalue). The matrix M_inv is then the null matrix. 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int sym_matrix_invert_regularised(
    struct sym_matrix *restrict M_inv, const struct sym_matrix *restrict M,
    const double max_cond_num, int *restrict regularised) {

  *regularised = 0;

  /* Non-finite input: fail (the checks below cannot be trusted with NaNs) */
  if (!sym_matrix_is_finite(M)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

#if defined(HYDRO_DIMENSION_3D)

  float A[3][3];
  get_matrix_from_sym_matrix(A, M);

  double A_d[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      A_d[i][j] = A[i][j];
    }
  }

  double ev[3], V[3][3];
  matrix_3x3_symmetric_eigendecomposition(A_d, ev, V);

  const double ev_max = fmax(ev[0], fmax(ev[1], ev[2]));
  if (!(ev_max > 0.)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  /* Regularise the spectrum below the floor: 1/lambda above it, a linear
   * ramp lambda / floor^2 below it (continuous at the floor, zero for an
   * exactly degenerate direction). */
  const double ev_floor = ev_max / max_cond_num;
  double ev_inv[3];
  for (int k = 0; k < 3; ++k) {
    if (ev[k] < ev_floor) {
      ev_inv[k] = fmax(ev[k], 0.) / (ev_floor * ev_floor);
      *regularised = 1;
    } else {
      ev_inv[k] = 1. / ev[k];
    }
  }

  /* Well-conditioned case: use the LU inverse, identical to the result of
   * sym_matrix_invert() */
  if (!*regularised) {
    double M_inv_matrix[3][3];
    if (invert3x3_matrix_LU(A_d, M_inv_matrix)) {
      zero_sym_matrix(M_inv);
      return 1;
    }
    M_inv->xx = M_inv_matrix[0][0];
    M_inv->yy = M_inv_matrix[1][1];
    M_inv->zz = M_inv_matrix[2][2];
    M_inv->xy = M_inv_matrix[0][1];
    M_inv->xz = M_inv_matrix[0][2];
    M_inv->yz = M_inv_matrix[1][2];
  } else {

    /* M_inv = V diag(ev_inv) V^T */
    double Minv[3][3];
    for (int i = 0; i < 3; ++i) {
      for (int j = i; j < 3; ++j) {
        double m = 0.;
        for (int k = 0; k < 3; ++k) m += V[i][k] * ev_inv[k] * V[j][k];
        Minv[i][j] = m;
      }
    }
    M_inv->xx = Minv[0][0];
    M_inv->yy = Minv[1][1];
    M_inv->zz = Minv[2][2];
    M_inv->xy = Minv[0][1];
    M_inv->xz = Minv[0][2];
    M_inv->yz = Minv[1][2];
  }

  /* The (float) result must be finite too */
  if (!sym_matrix_is_finite(M_inv)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  return 0;

#elif defined(HYDRO_DIMENSION_2D)

  const double a = M->xx;
  const double b = M->xy;
  const double c = M->yy;

  /* Eigenvalues (the small one from the determinant to avoid cancellation) */
  const double mean = 0.5 * (a + c);
  const double half_diff = 0.5 * (a - c);
  const double disc = sqrt(half_diff * half_diff + b * b);
  const double ev_big = mean + disc;
  const double ev_small = (ev_big > 0.) ? (a * c - b * b) / ev_big : -1.;

  if (!(ev_big > 0.)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  const double ev_floor = ev_big / max_cond_num;

  if (ev_small >= ev_floor) {

    /* Well-conditioned case: direct inverse, identical to sym_matrix_invert() */
    const double det_inv = 1. / (a * c - b * b);
    M_inv->xx = c * det_inv;
    M_inv->yy = a * det_inv;
    M_inv->xy = -b * det_inv;

  } else {

    *regularised = 1;

    /* Eigenvector of ev_big: (M - ev_big I) e = 0. Use the row of
     * M - ev_big I with the larger norm for a well-defined direction. */
    double ex, ey;
    const double r1x = a - ev_big, r1y = b;
    const double r2x = b, r2y = c - ev_big;
    if (r1x * r1x + r1y * r1y >= r2x * r2x + r2y * r2y) {
      ex = -r1y;
      ey = r1x;
    } else {
      ex = -r2y;
      ey = r2x;
    }
    const double norm = sqrt(ex * ex + ey * ey);
    if (norm > 0.) {
      ex /= norm;
      ey /= norm;
    } else {
      /* M is a multiple of the identity: any direction will do */
      ex = 1.;
      ey = 0.;
    }

    /* M_inv = e e^T / ev_big + e_perp e_perp^T * (ev_small / ev_floor^2)
     * with e_perp = (-ey, ex): linear ramp below the floor */
    const double inv_big = 1. / ev_big;
    const double inv_floor = fmax(ev_small, 0.) / (ev_floor * ev_floor);
    M_inv->xx = ex * ex * inv_big + ey * ey * inv_floor;
    M_inv->yy = ey * ey * inv_big + ex * ex * inv_floor;
    M_inv->xy = ex * ey * (inv_big - inv_floor);
  }

  /* The (float) result must be finite too */
  if (!sym_matrix_is_finite(M_inv)) {
    zero_sym_matrix(M_inv);
    return 1;
  }

  return 0;

#elif defined(HYDRO_DIMENSION_1D)

  return sym_matrix_invert(M_inv, M, max_cond_num);

#else
#error "A problem dimensionality must be chosen in config.h !"
#endif
}

#endif /* SWIFT_SYMMETRIC_MATRIX_H */
