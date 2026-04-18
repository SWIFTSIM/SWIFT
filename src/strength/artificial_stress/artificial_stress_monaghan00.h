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
#ifndef SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H
#define SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H

/**
 * @file strength/artificial_stress/artificial_stress_monaghan00.h
 * @brief Monaghan (2000) artificial stress method.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Apply artificial stress to pairwise stress tensors.
 *
 * Adds a negative factor to positive elements of the stress tensor. These
 * factors consist of the corresponding element of the stress tensor multiplied
 * by a factor that increases with reduced particle separation, scaling by the
 * kernel function to a given power. This method is presented by Monaghan2000
 * and is used by e.g. Schäfer+2016. Note that the expression for delta_p differs
 * from the one used in other papers, so that it can be computed more easily,
 * so direct comparisons can't be made.
 *
 * Method parameters needed in material parameter file:
 * ArtificialStress:
 *     n: Kernel exponent in particle separation factor.
 *     epsilon: Artificial stress multiplication factor.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interaction with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interaction with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
 */
__attribute__((always_inline)) INLINE static void artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r) {

  /* Method parameters. */
  const float artif_stress_n = method_artif_stress_n();
  const float artif_stress_epsilon = method_artif_stress_epsilon();

  /* Calculate the expected separation of closest neighbour, delta_p.
   * Note that the expression for delta_p differs from the one used in other
   * papers, so direct comparisons can't be made. */
  const float max_h = fmaxf(pi->h, pj->h);
  const float delta_p = max_h / 1.487; // ### hardcoded kernel eta for now

  /* Calculate separation factor, artif_stress_f. */
  float wij_delta_p, wij_r;
  kernel_eval(delta_p / max_h, &wij_delta_p);
  kernel_eval(r / max_h, &wij_r);
  const float artif_stress_f = powf(wij_r / wij_delta_p, artif_stress_n);

  /* Apply artificial stress. */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (pairwise_stress_tensor_i[i][j] > 0.f) {
        pairwise_stress_tensor_i[i][j] -= artif_stress_f * artif_stress_epsilon * pairwise_stress_tensor_i[i][j];
      }
      if (pairwise_stress_tensor_j[i][j] > 0.f) {
        pairwise_stress_tensor_j[i][j] -= artif_stress_f * artif_stress_epsilon * pairwise_stress_tensor_j[i][j];
      }
    }
  }
}

#endif /* SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H */
