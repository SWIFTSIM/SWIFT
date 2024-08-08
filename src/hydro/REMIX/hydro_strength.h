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
#ifndef SWIFT_PLANETARY_STRENGTH_H
#define SWIFT_PLANETARY_STRENGTH_H

/**
 * @file Planetary/hydro_strength.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the pressure.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_stress_tensor(
    struct part *restrict p, const float pressure) {

    for (int i = 0; i < 6; i++) {
        p->stress_tensor.elements[i] = p->deviatoric_stress_tensor.elements[i];
    }

    float pressure_hat = pressure;

    /*
    // If strength, if damage
    // See schafer 2016 for this...
    if (pressure < 0.f){
        pressure_hat *= (1 - p->damage_D);
    }
    */

    p->stress_tensor.xx -= pressure_hat;
    p->stress_tensor.yy -= pressure_hat;
    p->stress_tensor.zz -= pressure_hat;
}

/**
 * @brief Prepares extra strength parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_extra_strength(struct part *restrict p) {}

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
 * @brief Finishes extra strength parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_extra_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_extra_strength(struct part *restrict p) {}

/**
 * @brief Extra strength gradient interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_gradient_extra_strength(struct part *restrict pi,
                                          struct part *restrict pj,
                                          const float dx[3], const float wi,
                                          const float wj, const float wi_dx,
                                          const float wj_dx) {}

/**
 * @brief Extra strength gradient interaction between two particles
 * (non-symmetric)
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_nonsym_gradient_extra_strength(struct part *restrict pi,
                                                 const struct part *restrict pj,
                                                 const float dx[3],
                                                 const float wi,
                                                 const float wi_dx) {}

/**
 * @brief Finishes extra strength parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_extra_strength(struct part *restrict p) {}

/**
 * @brief Calculates the stress tensor with strength for the force interaction.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_pairwise_stress_tensors_strength(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj) {
    // Use the full stress tensor for solid particles
    if ((pi->phase_state == eos_phase_state_solid) &&
        (pj->phase_state == eos_phase_state_solid)) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                pairwise_stress_tensor_i[i][j] = pi->stress_tensor[i][j];
                pairwise_stress_tensor_j[i][j] = pj->stress_tensor[i][j];
            }
        }

        // Artificial stress etc...
    }
}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_strength(struct part *restrict p) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->dv_lin_repr_kernel[i][j] = 0.f;
    }
  }

  hydro_set_stress_tensor(p, pressure);
}

/**
 * @brief Extra strength force interaction between two particles
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_runner_iact_force_extra_strength(struct part *restrict pi,
                                       struct part *restrict pj,
                                       const float dx[3], const float Gi[3],
                                       const float Gj[3]) {

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == eos_phase_state_solid) &&
      (pj->phase_state == eos_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_lin_repr_kernel[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
        pj->dv_lin_repr_kernel[i][j] +=
            (pi->v[j] - pj->v[j]) * (-Gj[i]) * (pi->mass / pi->rho_evol);
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
                                              const float Gi[3]) {

  // Compute velocity gradient if both particles are solid
  if ((pi->phase_state == eos_phase_state_solid) &&
      (pj->phase_state == eos_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_lin_repr_kernel[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
      }
    }
  }
}

/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_extra_strength(struct part *restrict p) {

  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3], rotation_term[3][3];

  // Set the strain and rotation rates
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      strain_rate_tensor[i][j] =
          0.5f * (p->dv_lin_repr_kernel[i][j] + p->dv_lin_repr_kernel[j][i]);
      rotation_rate_tensor[i][j] =
          0.5f * (p->dv_lin_repr_kernel[j][i] - p->dv_lin_repr_kernel[i][j]);
    }
  }

  // Set rotation to transform into the solid body rotation frame
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      rotation_term[i][j] = 0.f;

      for (int k = 0; k < 3; k++) {
        // See Dienes 1978 (eqn 4.8)
        rotation_term[i][j] += p->deviatoric_stress_tensor[i][k] *
                                   rotation_rate_tensor[k][j] -
                               rotation_rate_tensor[i][k] *
                                   p->deviatoric_stress_tensor[k][j];
      }
    }
  }

  // Compute time derivative of the deviatoric stress tensor (Hooke's law)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      p->dS_dt[i][j] = 2.0f * p->shear_modulus_mu * strain_rate_tensor[i][j] + rotation_term[i][j];
    }

    p->dS_dt[i][i] -= 2.0f * p->shear_modulus_mu *
        (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] + strain_rate_tensor[2][2]) / 3.f;
  }
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_deviatoric_stress(
    struct part *restrict p, const float dt_therm) {

  if (p->phase_state == eos_phase_state_fluid) {
    // No stress for fluids
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        p->deviatoric_stress_tensor[i][j] = 0.f;
      }
    }
  } else {
    // Update solid stress
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        p->deviatoric_stress_tensor[i][j] += p->dS_dt[i][j] * dt_therm;
      }
    }
  }
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength(
    struct part *restrict p, const float dt_therm) {

    evolve_deviatoric_stress(p, dt_therm);
}

#endif /* SWIFT_PLANETARY_STRENGTH_H */
