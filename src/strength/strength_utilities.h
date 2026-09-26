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
 * @brief Adds the contribution of a neighbour to a particle's velocity gradient.
 *
 * Used in the force loop to build dv_force_loop, with dv[i][j] = dv_j/dx_i.
 *
 * @param dv The velocity gradient contribution to add to dv/dr.
 * @param vi Velocity of the particle.
 * @param vj Velocity of the neighbour.
 * @param G Kernel gradient for the particle pair.
 * @param volume_j Volume of the neighbour.
 */
__attribute__((always_inline)) INLINE static void
strength_add_velocity_gradient_contribution(float dv[3][3], const float vi[3],
                                    const float vj[3], const float G[3],
                                    const float volume_j) {

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      dv[i][j] += (vj[j] - vi[j]) * G[i] * volume_j;
    }
  }
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
 * The rotation rate tensor is R_ij = 0.5 * (dv_i/dx_j - dv_j/dx_i).
 *
 * @param rotation_rate_tensor The rotation rate tensor to be computed.
 * @param dv The velocity gradient dv/dr, with dv[i][j] = dv_j/dx_i
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
 * This function calculates the term R*M - M*R, where R is the rotation rate tensor,
 * R_ij = 0.5 * (dv_i/dx_j - dv_j/dx_i), and M is the tensor being rotated. This
 * term accounts for the apparent change of the tensor due to rotation of the
 * reference frame.
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
        rotation_term[i][j] += rotation_rate_tensor[i][k] * M[k][j] -
                               M[i][k] * rotation_rate_tensor[k][j];
      }
    }
  }
}

#endif /* SWIFT_STRENGTH_UTILITIES_H */