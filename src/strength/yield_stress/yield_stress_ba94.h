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
#ifndef SWIFT_YIELD_STRESS_BA94_H
#define SWIFT_YIELD_STRESS_BA94_H

/**
 * @file strength/yield_stress/yield_stress_ba94.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Yield stress method dependent way of adding contribution of damage to deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_apply_damage_to_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, const float damage) {

  // See BA94 Eqn. 14    
  deviatoric_stress_tensor.xx *= (1.f - damage);
  deviatoric_stress_tensor.yy *= (1.f - damage);
  deviatoric_stress_tensor.zz *= (1.f - damage);
  deviatoric_stress_tensor.xy *= (1.f - damage);
  deviatoric_stress_tensor.xz *= (1.f - damage);
  deviatoric_stress_tensor.yz *= (1.f - damage);

  return deviatoric_stress_tensor;
}

/**
 * @brief Adjust the yield stress depending on the damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
yield_apply_damage_to_yield_stress(const float yield_stress_intact, const float yield_stress_fully_damaged, const float damage) {

  return yield_stress_intact;
}

/**
 * @brief Compute the intact yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_intact(
    const int mat_id, const int phase_state, const float pressure) {

  float yield_stress_intact = 0.f;

  if (phase_state != mat_phase_state_fluid) {
    // Constant yield stress
    yield_stress_intact = material_Y_0(mat_id);
  }

  return yield_stress_intact;
}

/**
 * @brief Compute the fully damaged yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const int mat_id, const int phase_state, const float pressure,
    const float yield_stress_intact) {

  // Constant yield stress. 
  // Instead damage appears in effective_deviatoric_stress_from_damage()
  return yield_compute_yield_stress_intact(mat_id, phase_state, pressure);
}

/**
 * @brief Calculates the yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress(
    const int mat_id, const int phase_state, const float density, const float u, const float damage) {

  float yield_stress = 0.f;

  if (phase_state != mat_phase_state_fluid) {
    const float pressure =
      gas_pressure_from_internal_energy(density, u, mat_id);    
      
    // Constant yield stress
    yield_stress = yield_compute_yield_stress_intact(mat_id, phase_state, pressure);

    yield_stress =
        yield_softening_apply_density_to_yield_stress(yield_stress, mat_id, density);

    yield_stress =
        yield_softening_apply_temperature_to_yield_stress(yield_stress, mat_id, density, u);
  }

  return yield_stress;
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_deviatoric_stress_tensor(
    struct sym_matrix *deviatoric_stress_tensor,
    const float yield_stress, const float density, const float u) {

  float J_2 = strength_compute_stress_tensor_J_2(deviatoric_stress_tensor);

  // ...
  float f = fminf((yield_stress * yield_stress) / (3.f * J_2), 1.f);

  //## should have some dt dependence?
  for (int i = 0; i < 6; i++) deviatoric_stress_tensor->elements[i] *= f;
}

#endif /* SWIFT_YIELD_STRESS_BA94_H */
