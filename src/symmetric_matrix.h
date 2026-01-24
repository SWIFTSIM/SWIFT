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
 * @brief Construc the identity matrix
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
  M->xz = 0.f;
#endif
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
  ut[0] = M->xx * v[0];
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
  message("|%.4f %.4f %.4f|", M->xx, M->xy, M->xz);
  message("|%.4f %.4f %.4f|", M->xy, M->yy, M->yz);
  message("|%.4f %.4f %.4f|", M->xz, M->yz, M->zz);
#elif defined(HYDRO_DIMENSION_2D)
  message("|%.4f %.4f|", M->xx, M->xy);
  message("|%.4f %.4f|", M->xy, M->yy);
#elif defined(HYDRO_DIMENSION_1D)
  message("|%.4f|", M->xx);
#else
#error "A problem dimensionality must be chosen in config.h !"
#endif
}

/**
 * @brief Compute the inverse of a symmetric matrix.
 *
 * Returns as a symmetric matrix
 */
__attribute__((always_inline)) INLINE static void sym_matrix_invert(
    struct sym_matrix *restrict M_inv, const struct sym_matrix *restrict M) {

  float M_inv_matrix[hydro_dimension_integer][hydro_dimension_integer];
  get_matrix_from_sym_matrix(M_inv_matrix, M);
  const int res = invert_dimension_by_dimension_matrix(M_inv_matrix);
  if (res) {
    sym_matrix_print(M);
    error("Error inverting matrix");
  }

  M_inv->xx = M_inv_matrix[0][0];
#if defined(HYDRO_DIMENSION_2D) || defined(HYDRO_DIMENSION_3D)
  M_inv->yy = M_inv_matrix[1][1];
  M_inv->xy = M_inv_matrix[0][1];
#endif
#if defined(HYDRO_DIMENSION_3D)
  M_inv->zz = M_inv_matrix[2][2];
  M_inv->xz = M_inv_matrix[0][2];
  M_inv->yz = M_inv_matrix[1][2];
#endif
}

#endif /* SWIFT_SYMMETRIC_MATRIX_H */
