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
#include "symmetric_matrix.h"

/**
 * @brief Computes the J_2 invariant of the deviatoric stress tensor.
 *
 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 */
__attribute__((always_inline)) INLINE static float strength_compute_stress_tensor_J_2(
    const struct sym_matrix deviatoric_stress_tensor) {

  // ### Does j_2 need to be decreased by a factor of (1 - damage)^2 for B&A?

  /* Compute J_2 invariant. */
  return 0.5f * deviatoric_stress_tensor.xx * deviatoric_stress_tensor.xx +
         0.5f * deviatoric_stress_tensor.yy * deviatoric_stress_tensor.yy +
         0.5f * deviatoric_stress_tensor.zz * deviatoric_stress_tensor.zz +
         deviatoric_stress_tensor.xy * deviatoric_stress_tensor.xy +
         deviatoric_stress_tensor.xz * deviatoric_stress_tensor.xz +
         deviatoric_stress_tensor.yz * deviatoric_stress_tensor.yz;
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
  strain_rate_tensor[0][0] = 0.5f * (dv[0][0] + dv[0][0]);
  strain_rate_tensor[0][1] = 0.5f * (dv[0][1] + dv[1][0]);
  strain_rate_tensor[0][2] = 0.5f * (dv[0][2] + dv[2][0]);
  strain_rate_tensor[1][0] = 0.5f * (dv[1][0] + dv[0][1]);
  strain_rate_tensor[1][1] = 0.5f * (dv[1][1] + dv[1][1]);
  strain_rate_tensor[1][2] = 0.5f * (dv[1][2] + dv[2][1]);
  strain_rate_tensor[2][0] = 0.5f * (dv[2][0] + dv[0][2]);
  strain_rate_tensor[2][1] = 0.5f * (dv[2][1] + dv[1][2]);
  strain_rate_tensor[2][2] = 0.5f * (dv[2][2] + dv[2][2]);
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
  rotation_rate_tensor[0][0] = 0.5f * (dv[0][0] - dv[0][0]);
  rotation_rate_tensor[0][1] = 0.5f * (dv[1][0] - dv[0][1]);
  rotation_rate_tensor[0][2] = 0.5f * (dv[2][0] - dv[0][2]);
  rotation_rate_tensor[1][0] = 0.5f * (dv[0][1] - dv[1][0]);
  rotation_rate_tensor[1][1] = 0.5f * (dv[1][1] - dv[1][1]);
  rotation_rate_tensor[1][2] = 0.5f * (dv[2][1] - dv[1][2]);
  rotation_rate_tensor[2][0] = 0.5f * (dv[0][2] - dv[2][0]);
  rotation_rate_tensor[2][1] = 0.5f * (dv[1][2] - dv[2][1]);
  rotation_rate_tensor[2][2] = 0.5f * (dv[2][2] - dv[2][2]);
}

/**
 * @brief Computes the rotation term to transform the deviatoric stress tensor into the co-rotating frame.
 *
 * Note: Papers often make errors in the signs in this equation. For the correct
 *       equation, see Dienes 1979 for a detailed derivation, which leads to
 *       the final expression in Eqn. 4.8.
 *
 * @param rotation_term The rotation term to be computed.
 * @param rotation_rate_tensor The rotation rate tensor.
 * @param sym_matrix_deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 */
__attribute__((always_inline)) INLINE static void strength_compute_rotation_term(float rotation_term[3][3],
const float rotation_rate_tensor[3][3], const struct sym_matrix sym_matrix_deviatoric_stress_tensor) {

  /* Convert deviatoric stress to 3x3 float for matrix multiplication to compute
   * the rotation term . */
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &sym_matrix_deviatoric_stress_tensor);

  /* Compute rotation term elements. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;
      for (int k = 0; k < 3; k++) {
        rotation_term[i][j] +=
            deviatoric_stress_tensor[i][k] * rotation_rate_tensor[k][j] -
            rotation_rate_tensor[i][k] * deviatoric_stress_tensor[k][j];
      }
    }
  }
}

#endif /* SWIFT_STRENGTH_UTILITIES_H */
