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
#ifndef SWIFT_STRENGTH_STRESS_TENSOR_H
#define SWIFT_STRENGTH_STRESS_TENSOR_H

/**
 * @file strength/strength_stress_tensor.h
 * @brief Hooke's Law model for elastc stress 
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_utilities.h"

__attribute__((always_inline)) INLINE static void strength_compute_timestep_stress_tensor(
    const struct part *restrict p, const struct hydro_props *restrict hydro_properties, float *dt_cfl) {

  const float elastic_timestep_factor = hydro_properties->CFL_condition; // ### Set as same as CFL factor for now. Treat this similarly to CFL
  const float norm_dS_dt = norm_sym_matrix(&p->strength_data.dS_dt);
  const float shear_mod = material_shear_mod(p->mat_id);
    
  if (norm_dS_dt * *dt_cfl > elastic_timestep_factor * shear_mod) {
    *dt_cfl = elastic_timestep_factor * shear_mod / norm_dS_dt;
  }
}

__attribute__((always_inline)) INLINE static void
strength_compute_max_wave_speed_stress_tensor(const struct part *restrict p, const float soundspeed, const float density, float *wave_speed) {
  if (p->phase_state == mat_phase_state_solid) {
    const float shear_mod = material_shear_mod(p->mat_id);
    
    // Speed of longitudinal elastic wave
    *wave_speed = sqrtf(soundspeed * soundspeed + (4.f / 3.f) * shear_mod / density);
  }
}

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void strength_compute_stress_tensor(
    struct part *restrict p, const float pressure) {

  p->strength_data.stress_tensor = p->strength_data.deviatoric_stress_tensor;
  p->strength_data.stress_tensor.xx -= pressure;
  p->strength_data.stress_tensor.yy -= pressure;
  p->strength_data.stress_tensor.zz -= pressure;

  const float damage = strength_get_damage(p);
  const struct sym_matrix damaged_deviatoric_stress_tensor = yield_apply_damage_to_deviatoric_stress_tensor(p->strength_data.deviatoric_stress_tensor, damage);
  damage_compute_stress_tensor(&p->strength_data.stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);
}

/**
 * @brief Calculates the stress tensor with strength for the force interaction.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
strength_set_pairwise_stress_tensors(float pairwise_stress_tensor_i[3][3],
                                     float pairwise_stress_tensor_j[3][3],
                                     const struct part *restrict pi,
                                     const struct part *restrict pj,
                                     const float r) {

  // Use the full stress tensor for solid particles
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {

    get_matrix_from_sym_matrix(pairwise_stress_tensor_i, &pi->strength_data.stress_tensor);
    get_matrix_from_sym_matrix(pairwise_stress_tensor_j, &pj->strength_data.stress_tensor);

    artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i,
                                                               pairwise_stress_tensor_j, pi, pj, r);
  }
}

__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_stress_tensor(
    struct part *restrict p, const struct xpart *restrict xp) {
    
  p->strength_data.deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;
}

/**
 * @brief Calculate time derivative of the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stress_tensor_compute_dS_dt(struct part *restrict p, const float dv[3][3]) {

  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3],
      rotation_term[3][3];
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &p->strength_data.deviatoric_stress_tensor);

  // Set the strain and rotation rates
  strength_compute_strain_rate_tensor(dv, strain_rate_tensor);
  strength_compute_rotation_rate_tensor(dv, rotation_rate_tensor);
  strength_compute_rotation_term(rotation_term, rotation_rate_tensor,
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

  get_sym_matrix_from_matrix(&p->strength_data.dS_dt, dS_dt);
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void stress_tensor_evolve_deviatoric_stress_tensor(
    struct part *restrict p, struct sym_matrix *deviatoric_stress_tensor, const int phase_state, 
    float dt_therm) {

  if (phase_state == mat_phase_state_fluid) {
    // No stress for fluids
    zero_sym_matrix(deviatoric_stress_tensor);
  } else {
    // Update solid stress
    for (int i = 0; i < 6; i++) {
      deviatoric_stress_tensor->elements[i] +=
          p->strength_data.dS_dt.elements[i] * dt_therm;
    }
  }    
}

__attribute__((always_inline)) INLINE static void strength_first_init_part_stress_tensor(
    struct part *restrict p, struct xpart *restrict xp) {

  zero_sym_matrix(&p->strength_data.deviatoric_stress_tensor);
  zero_sym_matrix(&xp->strength_data.deviatoric_stress_tensor_full);
}

#endif /* SWIFT_STRENGTH_STRESS_TENSOR_H */
