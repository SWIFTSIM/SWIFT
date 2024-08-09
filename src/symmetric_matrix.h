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

/**
 * @brief Symmetric matrix definition in 3D.
 *
 * The matrix elements can be accessed as an array or via their "coordinates".
 * Replicated elements are not stored.
 */
struct sym_matrix {

  union {
    struct {
      float elements[6];
    };
    struct {
      float xx;
      float yy;
      float zz;
      float xy;
      float xz;
      float yz;
    };
  };
};

/**
 * @brief Zero the matrix
 */
__attribute__((always_inline)) INLINE static void zero_sym_matrix(
    struct sym_matrix *M) {
  for (int i = 0; i < 6; ++i) M->elements[i] = 0.f;
}

/**
 * @brief Construct a 3x3 array from a symmetric matrix.
 */
__attribute__((always_inline)) INLINE static void get_matrix_from_sym_matrix(
    float out[3][3], const struct sym_matrix *in) {

  out[0][0] = in->xx;
  out[0][1] = in->xy;
  out[0][2] = in->xz;
  out[1][0] = in->xy;
  out[1][1] = in->yy;
  out[1][2] = in->yz;
  out[2][0] = in->xz;
  out[2][1] = in->yz;
  out[2][2] = in->zz;
}

/**
 * @brief Construct a symmetric matrix from a 3x3 array.
 *
 * No check is performed to verify the input 3x3 array is indeed symmetric.
 */
__attribute__((always_inline)) INLINE static void get_sym_matrix_from_matrix(
    struct sym_matrix *out, const float in[3][3]) {
  out->xx = in[0][0];
  out->yy = in[1][1];
  out->zz = in[2][2];
  out->xy = in[0][1];
  out->xz = in[0][2];
  out->yz = in[1][2];
}

/**
 * @brief Compute the product of a symmetric matrix and a vector.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_multiply_by_vector(
    float out[3], const struct sym_matrix *M, const float v[3]) {

  out[0] = M->xx * v[0] + M->xy * v[1] + M->xz * v[2];
  out[1] = M->xy * v[0] + M->yy * v[1] + M->yz * v[2];
  out[2] = M->xz * v[0] + M->yz * v[1] + M->zz * v[2];
}

/**
 * @brief Multiply two symmetric matrices in two operations, ABA.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_multiplication_ABA(
    struct sym_matrix *M_out, const struct sym_matrix *A,
    const struct sym_matrix *B) {

  M_out->xx = (A->xx * (A->xx * B->xx + A->xy * B->xy + A->xz * B->xz) +
               A->xy * (A->xx * B->xy + A->xy * B->yy + A->xz * B->yz) +
               A->xz * (A->xx * B->xz + A->xy * B->yz + A->xz * B->zz));
  M_out->yy = (A->xy * (A->xy * B->xx + A->yy * B->xy + A->yz * B->xz) +
               A->yy * (A->xy * B->xy + A->yy * B->yy + A->yz * B->yz) +
               A->yz * (A->xy * B->xz + A->yy * B->yz + A->yz * B->zz));
  M_out->zz = (A->xz * (A->xz * B->xx + A->yz * B->xy + A->zz * B->xz) +
               A->yz * (A->xz * B->xy + A->yz * B->yy + A->zz * B->yz) +
               A->zz * (A->xz * B->xz + A->yz * B->yz + A->zz * B->zz));
  M_out->xy = (A->xy * (A->xx * B->xx + A->xy * B->xy + A->xz * B->xz) +
               A->yy * (A->xx * B->xy + A->xy * B->yy + A->xz * B->yz) +
               A->yz * (A->xx * B->xz + A->xy * B->yz + A->xz * B->zz));
  M_out->xz = (A->xz * (A->xx * B->xx + A->xy * B->xy + A->xz * B->xz) +
               A->yz * (A->xx * B->xy + A->xy * B->yy + A->xz * B->yz) +
               A->zz * (A->xx * B->xz + A->xy * B->yz + A->xz * B->zz));
  M_out->yz = (A->xz * (A->xy * B->xx + A->yy * B->xy + A->yz * B->xz) +
               A->yz * (A->xy * B->xy + A->yy * B->yy + A->yz * B->yz) +
               A->zz * (A->xy * B->xz + A->yy * B->yz + A->yz * B->zz));
}

/**
 * @brief Print a symmetric matrix.
 */
__attribute__((always_inline)) INLINE static void sym_matrix_print(
    const struct sym_matrix *M) {
  message("|%.4f %.4f %.4f|", M->xx, M->xy, M->xz);
  message("|%.4f %.4f %.4f|", M->xy, M->yy, M->yz);
  message("|%.4f %.4f %.4f|", M->xz, M->yz, M->zz);
}

/**
 * @brief Compute the inverse of a symmetric matrix and a vector.
 *
 * Returned as a symmetric matrix
 */
__attribute__((always_inline)) INLINE static void sym_matrix_invert(
    struct sym_matrix *M_inv, const struct sym_matrix *M) {

  // Can probably do this by just typing out all the elements,
  // so we don't have to go back and forth between arrays and sym_matrix.
  // But, for now use invert_dimension_by_dimension_matrix function
  float M_inv_matrix[3][3];
  get_matrix_from_sym_matrix(M_inv_matrix, M);
  int res = invert_dimension_by_dimension_matrix(M_inv_matrix);
  if (res) {
    sym_matrix_print(M);
    error("Error inverting matrix");
  }

  get_sym_matrix_from_matrix(M_inv, M_inv_matrix);
}

/**
 * @brief Calculate eigenvalues of a sym_matrix
 *
 * ###Temp note: see calculate_all_eigenvalues in
 * https://github.com/christophmschaefer/miluphcuda/blob/main/linalg.cu
 *
 */
__attribute__((always_inline)) INLINE static void
sym_matrix_compute_eigenvalues(float eigenvalues[3],
                               const struct sym_matrix M_sym_matrix) {

  float M[3][3];
  get_matrix_from_sym_matrix(M, &M_sym_matrix);

  int i, j;
  // Current iteration towards diagonalisation
  float M_iter[3][3];
  // The largest (absolute value) off-diagonal element ...
  float max_offdiag;
  // ... and its indices
  int e, f;
  // trig functions corresponding to angle of rotation
  float sin_theta, cos_theta, tan_theta, one_over_tan_2theta;
  // rotated M_iter
  float M_iter_rotated[3][3];

  // init M_iter
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      M_iter[i][j] = M[i][j];
    }
  }

  // Note Schafer has limit of 5
  int max_iter = 5;
  float tol = 1e-10;
  // iterate until diagonalised within tolerance or until max_iter is reached
  for (int iter = 0; iter < max_iter; iter++) {

    max_offdiag = 0.f;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        // init M_iter_rotated
        M_iter_rotated[i][j] = M_iter[i][j];
        // Find max off-diagonal element and its indices
        if (i != j) {
          if (fabs(M_iter[i][j]) >= max_offdiag) {
            max_offdiag = fabs(M_iter[i][j]);
            e = i;
            f = j;
          }
        }
      }
    }

    // Has matrix has been diagonalised?
    if (max_offdiag < tol) break;

    // Get sin_theta and cos_theta
    one_over_tan_2theta = (M_iter[f][f] - M_iter[e][e]) / (2 * M_iter[e][f]);
    // check this
    if (one_over_tan_2theta < 0)
      tan_theta =
          -1.f / (fabs(one_over_tan_2theta) +
                  sqrtf(one_over_tan_2theta * one_over_tan_2theta + 1.f));
    else
      tan_theta =
          1.f / (fabs(one_over_tan_2theta) +
                 sqrtf(one_over_tan_2theta * one_over_tan_2theta + 1.f));

    cos_theta = 1.f / (sqrtf(tan_theta * tan_theta + 1));
    sin_theta = tan_theta * cos_theta;

    // Get M_iter_rotated by rotating M_iter
    M_iter_rotated[e][e] = cos_theta * cos_theta * M_iter[e][e] +
                           sin_theta * sin_theta * M_iter[f][f] -
                           2.f * sin_theta * cos_theta * M_iter[e][f];
    M_iter_rotated[f][f] = cos_theta * cos_theta * M_iter[f][f] +
                           sin_theta * sin_theta * M_iter[e][e] +
                           2.f * sin_theta * cos_theta * M_iter[e][f];
    M_iter_rotated[e][f] =
        (cos_theta * cos_theta - sin_theta * sin_theta) * M_iter[e][f] +
        sin_theta * cos_theta * (M_iter[e][e] - M_iter[f][f]);
    M_iter_rotated[f][e] = M_iter_rotated[e][f];

    /* the other element in column e and row f*/
    for (i = 0; i < 3; i++) {
      if (i != f && i != e) {
        M_iter_rotated[e][i] =
            cos_theta * M_iter[i][e] - sin_theta * M_iter[i][f];
        M_iter_rotated[i][e] = M_iter_rotated[e][i];
        M_iter_rotated[f][i] =
            cos_theta * M_iter[i][f] + sin_theta * M_iter[i][e];
        M_iter_rotated[i][f] = M_iter_rotated[f][i];
      }
    }

    /* update  M_iter to M_iter_rotated for next iter or output */
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
        M_iter[i][j] = M_iter_rotated[i][j];
      }
    }
  }

  eigenvalues[0] = M_iter[0][0];
  eigenvalues[1] = M_iter[1][1];
  eigenvalues[2] = M_iter[2][2];
}

#endif /* SWIFT_SYMMETRIC_MATRIX_H */
