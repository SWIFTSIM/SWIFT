/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_STRENGTH_UTILS_H
#define SWIFT_PLANETARY_STRENGTH_UTILS_H
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/strength_utilities.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Compute the J_2 invariant of the deviatoric stress tensor.
 */
__attribute__((always_inline)) INLINE static float J_2_from_stress_tensor(
    struct sym_matrix *sym_matrix_deviatoric_stress_tensor) {

  // ### Does j_2 need to be decreased by a factor of (1 - damage)^2 for B&A?
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             sym_matrix_deviatoric_stress_tensor);
  float J_2 = 0.f;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      J_2 += 0.5f * deviatoric_stress_tensor[i][j] *
              deviatoric_stress_tensor[j][i];
    }
  }

  return J_2;
}

/**
 * @brief Calculate the strain rate tensor.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_strain_rate_tensor(
    struct part *restrict p, float strain_rate_tensor[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_rate_tensor[i][j] =
          0.5f * (p->dv_force_loop[i][j] + p->dv_force_loop[j][i]);
    }
  }
}

/**
 * @brief Calculate the rotation rate tensor.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
calculate_rotation_rate_tensor(struct part *restrict p,
                               float rotation_rate_tensor[3][3]) {

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_rate_tensor[i][j] =
          0.5f * (p->dv_force_loop[j][i] - p->dv_force_loop[i][j]);
    }
  }
}

/**
 * @brief Calculate the rotation term to transform into the co-rotating frame.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_rotation_term(
    float rotation_term[3][3], const float rotation_rate_tensor[3][3],
    const float deviatoric_stress_tensor[3][3]) {

  // Set rotation to transform into the corotating frame
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;

      for (int k = 0; k < 3; k++) {
        // See Dienes 1978 (eqn 4.8)
        rotation_term[i][j] +=
            deviatoric_stress_tensor[i][k] * rotation_rate_tensor[k][j] -
            rotation_rate_tensor[i][k] * deviatoric_stress_tensor[k][j];
      }
    }
  }
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRENGTH_UTILS_H */
