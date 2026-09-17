/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_STRENGTH_UTILITIES_H
#define SWIFT_STRENGTH_UTILITIES_H

/**
 * @file strength/strength_utilities.h
 * @brief Utilities used throughout the material strength scheme.
 */

#include "math.h"
#include <string.h>

/**
 * @brief Computes the J_2 invariant of a deviatoric symmetric tensor.
 *
 * @param M The deviatoric symmetric matrix.
 */
__attribute__((always_inline)) INLINE static float strength_compute_deviatoric_sym_matrix_J_2(
    const struct sym_matrix M) {

  // ### Does j_2 need to be decreased by a factor of (1 - damage)^2 for B&A?

  /* Compute J_2 invariant. */
  return 0.5f * M.xx * M.xx + 0.5f * M.yy * M.yy + 0.5f * M.zz * M.zz +
         M.xy * M.xy + M.xz * M.xz + M.yz * M.yz;
}

/**
 * @brief Computes the J_2 invariant of a symmetric tensor.
 *
 * @param M The symmetric matrix.
 */
__attribute__((always_inline)) INLINE static float strength_compute_sym_matrix_J_2(
    const struct sym_matrix M) {

  /* Calculate deviatoric tensor. */
  struct sym_matrix M_dev = M;
  M_dev.xx -= (M.xx + M.yy + M.zz) / 3.f;
  M_dev.yy -= (M.xx + M.yy + M.zz) / 3.f;
  M_dev.zz -= (M.xx + M.yy + M.zz) / 3.f;

  /* Compute J_2 invariant. */
  return strength_compute_deviatoric_sym_matrix_J_2(M_dev);
}

/**
 * @brief Computes the strain rate tensor.
 *
 * @param strain_rate_tensor The strain rate tensor to be computed.
 * @param dv The velocity gradient dv/dr.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_strain_rate_tensor(float strain_rate_tensor[3][3], const float dv[3][3]) {

  /* Compute strain rate tensor elements. */
  strain_rate_tensor[0][0] = dv[0][0];
  strain_rate_tensor[1][1] = dv[1][1];
  strain_rate_tensor[2][2] = dv[2][2];
  strain_rate_tensor[0][1] = 0.5f * (dv[0][1] + dv[1][0]);
  strain_rate_tensor[0][2] = 0.5f * (dv[0][2] + dv[2][0]);
  strain_rate_tensor[1][0] = 0.5f * (dv[1][0] + dv[0][1]);
  strain_rate_tensor[1][2] = 0.5f * (dv[1][2] + dv[2][1]);
  strain_rate_tensor[2][0] = 0.5f * (dv[2][0] + dv[0][2]);
  strain_rate_tensor[2][1] = 0.5f * (dv[2][1] + dv[1][2]);
}

/**
 * @brief Computes the rotation rate tensor.
 *
 * @param rotation_rate_tensor The rotation rate tensor to be computed.
 * @param dv The velocity gradient dv/dr.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_rotation_rate_tensor(float rotation_rate_tensor[3][3], const float dv[3][3]) {

  /* Compute rotation rate tensor elements. */
  rotation_rate_tensor[0][0] = 0.f;
  rotation_rate_tensor[1][1] = 0.f;
  rotation_rate_tensor[2][2] = 0.f;
  rotation_rate_tensor[0][1] = 0.5f * (dv[1][0] - dv[0][1]);
  rotation_rate_tensor[0][2] = 0.5f * (dv[2][0] - dv[0][2]);
  rotation_rate_tensor[1][0] = 0.5f * (dv[0][1] - dv[1][0]);
  rotation_rate_tensor[1][2] = 0.5f * (dv[2][1] - dv[1][2]);
  rotation_rate_tensor[2][0] = 0.5f * (dv[0][2] - dv[2][0]);
  rotation_rate_tensor[2][1] = 0.5f * (dv[1][2] - dv[2][1]);
}

/**
 * @brief Computes the rotation contribution for transforming a tensor into the co-rotating frame.
 *
 * This function calculates the term M*R - R*M, where R is the rotation rate tensor,
 * and M is the tensor being rotated. This term accounts for the apparent change
 * of the tensor due to rotation of the reference frame.
 *
 * Note: Papers often make errors in the signs in this equation. For the correct
 *       equation, see Dienes 1979 for a detailed derivation, which leads to
 *       the final expression in Eqn. 4.8.
 *
 * @param rotation_term The rotation term to be computed.
 * @param rotation_rate_tensor The rotation rate tensor.
 * @param M The tensor to rotate.
 */
__attribute__((always_inline)) INLINE static void strength_compute_rotation_term(float rotation_term[3][3],
const float rotation_rate_tensor[3][3], const float M[3][3]) {

  /* Compute rotation term elements. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;
      for (int k = 0; k < 3; k++) {
        rotation_term[i][j] += M[i][k] * rotation_rate_tensor[k][j] -
                               rotation_rate_tensor[i][k] * M[k][j];
      }
    }
  }
}

#endif /* SWIFT_STRENGTH_UTILITIES_H */