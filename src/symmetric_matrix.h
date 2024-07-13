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
    
  M_inv->xx = M_inv_matrix[0][0];
  M_inv->yy = M_inv_matrix[1][1];
  M_inv->zz = M_inv_matrix[2][2];
  M_inv->xy = M_inv_matrix[0][1];
  M_inv->xz = M_inv_matrix[0][2];
  M_inv->yz = M_inv_matrix[1][2];
}

#endif /* SWIFT_SYMMETRIC_MATRIX_H */
