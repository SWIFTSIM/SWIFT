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
__attribute__((always_inline)) INLINE static struct sym_matrix yield_model_adjust_deviatoric_stress_tensor_by_damage(
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
adjust_yield_stress_by_damage(struct part *restrict p,
  const float yield_stress_intact, const float yield_stress_fully_damaged, const float damage) {

  return (1.f - damage) * yield_stress_intact + damage * yield_stress_fully_damaged;
}

/**
 * @brief Compute the intact yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_intact(
    struct part *restrict p, const int phase_state, const float pressure) {

  float yield_stress_intact = 0.f;

  if (phase_state != mat_phase_state_fluid) {

    const float mu_i = material_mu_i(p->mat_id);
    const float Y_0 = material_Y_0(p->mat_id);
    const float Y_M = material_Y_M(p->mat_id);
      
    yield_stress_intact = Y_0;

    // Should be able to decrease if negative pressures until yield_stress=0?
    // miluphcuda does this so maybe not wrong?
    if (pressure > 0.f) {
      if (Y_M != Y_0) {
        yield_stress_intact +=
            mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
      }
    }
  }

  return yield_stress_intact;
}

/**
 * @brief Compute the fully damaged yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_fully_damaged(
    struct part *restrict p, const int phase_state, const float pressure,
    const float yield_stress_intact) {

  float yield_stress_damaged = 0.f;

  if (phase_state != mat_phase_state_fluid) {

    const float mu_d = material_mu_d(p->mat_id);

    if (pressure > 0.f) {
      yield_stress_damaged = mu_d * pressure;
    }

    // Maybe yield_stress_damaged also needs have a max value indep of
    // yield_stress_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018

    yield_stress_damaged = min(yield_stress_damaged, yield_stress_intact);
  }

  return yield_stress_damaged;
}

/**
 * @brief Calculates the yield stress
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress(
    struct part *restrict p, const int phase_state, const float density, const float u, const float damage) {

  float yield_stress = 0.f;

  if (phase_state != mat_phase_state_fluid) {

    const float pressure =
      gas_pressure_from_internal_energy(density, u, p->mat_id);    
      
    float yield_stress_intact = compute_yield_stress_intact(p, phase_state, pressure);
    float yield_stress_fully_damaged =
        compute_yield_stress_fully_damaged(p, phase_state, pressure, yield_stress_intact);

    // ...
    yield_stress = 
        adjust_yield_stress_by_damage(p, yield_stress_intact, yield_stress_fully_damaged, damage);

    // ### This was previously only for intact.
    yield_stress =
        adjust_yield_stress_by_density(p, yield_stress, density);

    yield_stress =
        adjust_yield_stress_by_temperature(p, yield_stress, density, u);
  }

  return yield_stress;
}

/**
 * @brief Calculates the yield stress using p->strength_data.damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_damage(
    struct part *restrict p, const int phase_state, const float density, const float u) {

  const float damage = strength_get_damage(p);

  return compute_yield_stress(p, phase_state, density, u, damage);
}

/**
 * @brief Calculates the yield stress using xp->strength_data.damage_full
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float compute_yield_stress_damage_full(
    struct part *restrict p, struct xpart *restrict xp, const int phase_state, const float density, const float u) {

  const float damage = strength_get_damage_full(xp);
    
  return compute_yield_stress(p, phase_state, density, u, damage);
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
adjust_deviatoric_stress_tensor_by_yield_stress(
    struct part *restrict p, struct sym_matrix *deviatoric_stress_tensor,
    const float yield_stress, const float density, const float u) {

  float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);

  // ...
  float f = min(yield_stress / sqrtf(J_2), 1.f);

  // ## should have some dt dependence?
  for (int i = 0; i < 6; i++) deviatoric_stress_tensor->elements[i] *= f;
}

#endif /* SWIFT_YIELD_STRESS_COLLINS04_H */
