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
#ifndef SWIFT_YIELD_STRESS_COLLINS04_H
#define SWIFT_YIELD_STRESS_COLLINS04_H

/**
 * @file strength/yield_stress/yield_stress_collins04.h
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

  // This method does not apply damage to deviatoric_stress_tensor 
  return deviatoric_stress_tensor;
}

/**
 * @brief Adjust the yield stress depending on the damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_intact, const float yield_stress_fully_damaged, const float damage) {

  return (1.f - damage) * yield_stress_intact + damage * yield_stress_fully_damaged;
}

/**
 * @brief Compute the intact yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_intact(
    const int mat_id, const int phase_state, const float pressure) {

  if (phase_state == mat_phase_state_fluid) {
    return 0.f;
  }

  const float mu_i = material_mu_i(mat_id);
  const float Y_0 = material_Y_0(mat_id);
  const float Y_M = material_Y_M(mat_id);

  // Should be able to decrease if negative pressures until yield_stress=0?
  // miluphcuda does this so maybe not wrong?
  if (pressure > 0.f && Y_M != Y_0) {
    return Y_0 + mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
  } else {
    return Y_0;
  }
}

/**
 * @brief Compute the fully damaged yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const int mat_id, const int phase_state, const float pressure,
    const float yield_stress_intact) {

  if (phase_state == mat_phase_state_fluid) {
    return 0.f;
  }
      
  if (pressure > 0.f) {    
    // Maybe yield_stress_damaged also needs have a max value indep of
    // yield_stress_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018
    return fminf(material_mu_d(mat_id), yield_stress_intact);
  } else {
    return 0.f;
  }
}

/**
 * @brief Calculates the yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress(
    const int mat_id, const int phase_state, const float density, const float u, const float damage) {

  if (phase_state == mat_phase_state_fluid) {
    return 0.f;
  }

  const float pressure =
    gas_pressure_from_internal_energy(density, u, mat_id);    
      
  const float yield_stress_intact = yield_compute_yield_stress_intact(mat_id, phase_state, pressure);
  const float yield_stress_fully_damaged =
      yield_compute_yield_stress_fully_damaged(mat_id, phase_state, pressure, yield_stress_intact);

  // ...
  float yield_stress = 
      yield_compute_damaged_yield_stress(yield_stress_intact, yield_stress_fully_damaged, damage);

  // ### This was previously only for intact.
  yield_softening_apply_density_to_yield_stress(&yield_stress, mat_id, density);
  yield_softening_apply_temperature_to_yield_stress(&yield_stress, mat_id, density, u);

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
  float f = fminf(yield_stress / sqrtf(J_2), 1.f);

  // ## should have some dt dependence?
  deviatoric_stress_tensor->xx *= f;
  deviatoric_stress_tensor->yy *= f;
  deviatoric_stress_tensor->zz *= f;
  deviatoric_stress_tensor->xy *= f;
  deviatoric_stress_tensor->xz *= f;
  deviatoric_stress_tensor->yz *= f;
}

#endif /* SWIFT_YIELD_STRESS_COLLINS04_H */
