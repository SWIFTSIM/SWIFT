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
#ifndef SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H
#define SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H

/**
 * @file strength/artificial_stress/artificial_stress_monaghan00.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Update the pairwise stress tensors with artificial stress.
 */
__attribute__((always_inline)) INLINE static void strength_add_artif_stress(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r) {
    
  // Artificial stress (Monaghan, 2000)
  const float max_h = max(pi->h, pj->h);
  const float delta_p = max_h / 1.487; // ### hardcoded for now

  float wij_delta_p;
  kernel_eval(delta_p / max_h, &wij_delta_p);

  float wij_r;
  kernel_eval(r / max_h, &wij_r);

  // This factor should be set in extra parameter file
  const float artif_stress_n = method_artif_stress_n();
  const float artif_stress_f = powf(wij_r / wij_delta_p, artif_stress_n);

  // This factor should be set in extra parameter file
  const float artif_stress_epsilon = method_artif_stress_epsilon();

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
    if (pairwise_stress_tensor_i[i][j] > 0.f)
      pairwise_stress_tensor_i[i][j] -= artif_stress_f *
                  artif_stress_epsilon * pairwise_stress_tensor_i[i][j]; 
    if (pairwise_stress_tensor_j[i][j] > 0.f) 
      pairwise_stress_tensor_j[i][j] -= artif_stress_f * 
                  artif_stress_epsilon * pairwise_stress_tensor_j[i][j];
    }
  }
}

#endif /* SWIFT_ARTIFICIAL_STRESS_MONAGHAN00_H */