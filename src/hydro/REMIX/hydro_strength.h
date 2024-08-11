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
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/hydro_strength.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_utilities.h"
#include "strength_damage.h"
#include "strength_stress.h"
#include "strength_yield.h"

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
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_extra_strength(struct part *restrict p,
                                   const float pressure) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      p->dv_force_loop[i][j] = 0.f;
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
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_force_loop[i][j] +=
            (pj->v[j] - pi->v[j]) * Gi[i] * (pj->mass / pj->rho_evol);
        pj->dv_force_loop[i][j] +=
            (pi->v[j] - pj->v[j]) * Gj[i] * (pi->mass / pi->rho_evol);
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
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        pi->dv_force_loop[i][j] +=
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

  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3],
      rotation_term[3][3];
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &p->deviatoric_stress_tensor);

  // Set the strain and rotation rates
  calculate_strain_rate_tensor(p, strain_rate_tensor);
  calculate_rotation_rate_tensor(p, rotation_rate_tensor);
  calculate_rotation_term(rotation_term, rotation_rate_tensor,
                          deviatoric_stress_tensor);

  // Compute time derivative of the deviatoric stress tensor (Hooke's law)
  const float shear_mod = material_shear_mod(p->mat_id);
  float dS_dt[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dS_dt[i][j] =
          2.0f * shear_mod * strain_rate_tensor[i][j] + rotation_term[i][j];
    }

    dS_dt[i][i] -= 2.0f * shear_mod *
                   (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                    strain_rate_tensor[2][2]) /
                   3.f;
  }

  get_sym_matrix_from_matrix(&p->dS_dt, dS_dt);
}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_extra_strength(
    struct part *restrict p, const float dt_therm, const float density,
    const float pressure) {

  evolve_damage(p, pressure, dt_therm);

  evolve_deviatoric_stress(p, dt_therm, density, pressure);
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRENGTH_H */
