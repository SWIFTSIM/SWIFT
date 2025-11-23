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
 * @brief Calculates the pairwise stress tensors for the force interaction.
 *
 * The stress tensors used for the force interaction between a specific pair of
 * particles. These differ from the particle's own stress tensor, since they
 * factor in the phases of the two particles as well as the contribution of
 * artificial stress for the pairwie interaction.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interactiion with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interactiion with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
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
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param wi The value of the unmodified kernel function W(r, hi) * hi^d.
 * @param wj The value of the unmodified kernel function W(r, hj) * hj^d.
 * @param wi_dx The norm of the gradient of wi: dW(r, hi)/dr * hi^(d+1).
 * @param wj_dx The norm of the gradient of wj: dW(r, hj)/dr * hj^(d+1).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_density_strength(struct part *restrict pi,
                                         struct part *restrict pj,
                                         const float dx[3], const float wi,
                                         const float wj, const float wi_dx,
                                         const float wj_dx) {}

/**
 * @brief Extra strength density interaction between two particles
 * (non-symmetric)
 *
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param wi The value of the unmodified kernel function W(r, hi) * hi^d.
 * @param wi_dx The norm of the gradient of wi: dW(r, hi)/dr * hi^(d+1).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_density_strength(struct part *restrict pi,
                                                const struct part *restrict pj,
                                                const float dx[3],
                                                const float wi,
                                                const float wi_dx) {}

/**
 * @brief Extra strength force interaction between two particles
 *
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param Gi Kernel gradient for first particle.
 * @param Gj Kernel gradient for second particle.
 * @param dv_dot_G_i Dot product of kernel gradient for first particle and (vi - vj).
 * @param dv_dot_G_j Dot product of kernel gradient for second particle and (vi - vj).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_force_strength(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float Gi[3],
                                       const float Gj[3], const float dv_dot_G_i,
                                       const float dv_dot_G_j) {

  /* Add contribution to drho/dt. */
  pi->strength_data.drho_dt += pj->mass * dv_dot_G_i;
  pj->strength_data.drho_dt += pi->mass * dv_dot_G_j;

  /* Add constribution to dv/dr if both particles are solid. */
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
 * @param pi First particle.
 * @param pj Second particle.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param Gi Kernel gradient for first particle.
 * @param Gj Kernel gradient for second particle.
 * @param dv_dot_G_i Dot product of kernel gradient for first particle and (vi - vj).
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_force_strength(struct part *restrict pi,
                                              const struct part *restrict pj,
                                              const float dx[3],
                                              const float Gi[3], const float dv_dot_G_i) {

  /* Add contribution to drho/dt. */
  pi->strength_data.drho_dt += pj->mass * dv_dot_G_i;

  /* Add constribution to dv/dr if both particles are solid. */
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
