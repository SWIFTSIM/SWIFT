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
#ifndef SWIFT_PLANETARY_STRENGTH_IACT_DEFAULT_H
#define SWIFT_PLANETARY_STRENGTH_IACT_DEFAULT_H

/**
 * @file Planetary/strength/default/hydro_iact_strength.h
 * @brief Planetary implementation of SPH with default material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength.h"

/**
 * @brief Calculates the stress tensor with strength for the force interaction.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_set_pairwise_stress_tensors_strength(float pairwise_stress_tensor_i[3][3],
                                           float pairwise_stress_tensor_j[3][3],
                                           const struct part *restrict pi,
                                           const struct part *restrict pj,
                                           const float r) {

  strength_set_pairwise_stress_tensors(pairwise_stress_tensor_i, pairwise_stress_tensor_j, pi, pj, r);
}


/**
 * @brief Extra strength density interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_extra_strength(struct part *restrict pi,
                                         struct part *restrict pj,
                                         const float dx[3], const float wi,
                                         const float wj, const float wi_dx,
                                         const float wj_dx) {}

/**
 * @brief Extra strength density interaction between two particles
 * (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_extra_strength(struct part *restrict pi,
                                                const struct part *restrict pj,
                                                const float dx[3],
                                                const float wi,
                                                const float wi_dx) {}

/**
 * @brief Extra strength force interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_force_extra_strength(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float Gi[3],
                                       const float Gj[3], const float dv_dot_G_i, 
                                       const float dv_dot_G_j) {
    
  pi->strength_data.drho_dt += mj * dv_dot_G_i;
  pj->strength_data.drho_dt += mi * dv_dot_G_j;

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->strength_data.dv_force_loop[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->strength_data.rho_evol);
        pj->strength_data.dv_force_loop[i][j] +=
            (pi->v[j] - pj->v[j]) * Gj[i] * (pi->mass / pi->strength_data.rho_evol);
      }
    }
  }
}

/**
 * @brief Extra strength force interaction between two particles (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_force_extra_strength(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3],
                                              const float Gi[3], const float dv_dot_G_i) {
  
  pi->strength_data.drho_dt += mj * dv_dot_G_i;

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->strength_data.dv_force_loop[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->strength_data.rho_evol);
      }
    }
  }
}

#endif /* SWIFT_PLANETARY_STRENGTH_IACT_DEFAULT_H */
