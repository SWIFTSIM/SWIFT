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
#ifndef SWIFT_YIELD_STRESS_BA95_H
#define SWIFT_YIELD_STRESS_BA95_H

/**
 * @file strength/yield_stress/yield_stress_ba95.h
 * Benz&Asphaug95 constant yield stress.
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Compute damaged deviatoric stress tensor.
 *
 * Yield stress method dependent way of adding contribution of damage to
 * deviatoric stress tensor. Note that this does not "apply" damage to the
 * deviatoric stress tensor, since the deviatoric stress tensor that is evolved
 * in time doesn't get modified by the damage. Instead this function returns a
 * new tensor, which gets used, for example, to construct the pairwise stress
 * tensors for the interaction of a particle pair.
 *
 * With the Benz&Asphaug95 yield stress method, damage acts to reduce all
 * elements of the deviatoric stress tensor, meaning that a fully damaged
 * material acts as afluid. See Benz&Asphaug95 Eqn. 13.

 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, const float damage) {

  /* Compute elements of damaged deviatoric stress tensor. */
  struct sym_matrix damaged_deviatoric_stress_tensor = deviatoric_stress_tensor;
  damaged_deviatoric_stress_tensor.xx *= (1.f - damage);
  damaged_deviatoric_stress_tensor.yy *= (1.f - damage);
  damaged_deviatoric_stress_tensor.zz *= (1.f - damage);
  damaged_deviatoric_stress_tensor.xy *= (1.f - damage);
  damaged_deviatoric_stress_tensor.xz *= (1.f - damage);
  damaged_deviatoric_stress_tensor.yz *= (1.f - damage);

  return damaged_deviatoric_stress_tensor;
}

/**
 * @brief Compute damaged yield stress.
 *
 * With the Benz&Asphaug95 yield stress method, damage does not affect the yield
 * stress. Damage is instead applied to the deviatoric stress tensor.
 *
 * @param yield_stress_fully_intact The yield stress for fully intact material.
 * @param yield_stress_fully_damaged The yield stress for fully damaged material.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact, const float yield_stress_fully_damaged, const float damage) {

  return yield_stress_fully_intact;
}

/**
 * @brief Compute the yield stress of fully intact material.
 *
 * Benz&Asphaug95 constant yield stress of solids.
 *
 * Method parameters needed in material parameter file:
 * YieldStress:
 *    Y_0: Constant yield stress (Pa).
 *
 * @param mat_id The material ID.
 * @param phase_state The phase ID.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_intact(
    const int mat_id, const int phase_state, const float pressure) {

  /* Return 0.f if the material is not solid. */
  if (phase_state != mat_phase_state_solid) {
    return 0.f;
  }

  /* Return constant yield stress, Y_0. */
  return material_Y_0(mat_id);
}

/**
 * @brief Compute the yield stress of fully damaged material.
 *
 * With the Benz&Asphaug95 yield stress method, damage does not affect the yield
 * stress. Damage is instead applied to the deviatoric stress tensor.
 *
 * @param mat_id The material ID.
 * @param phase_state The phase ID.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const int mat_id, const int phase_state, const float pressure) {

  /* Damage does not affect the yield stress. */
  return yield_compute_yield_stress_fully_intact(mat_id, phase_state, pressure);
}

/**
 * @brief Compute the yield stress.
 *
 * Calculates the yield stress that combines all relevant methods e.g. damage
 * and weakening.
 *
 * @param mat_id The material ID.
 * @param phase_state The phase ID.
 * @param density The density.
 * @param u The specific internal energy.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress(
    const int mat_id, const int phase_state, const float density, const float u, const float damage) {

  /* Return 0.f if the material is not solid. */
  if (phase_state != mat_phase_state_solid) {
    return 0.f;
  }

  /* Calculate pressure. */
  const float pressure =
    gas_pressure_from_internal_energy(density, u, mat_id);

  /* Get constant yield stress. */
  float yield_stress = yield_compute_yield_stress_fully_intact(mat_id, phase_state, pressure);

  /* Apply weakening to yield stress. */
  yield_weakening_apply_density_to_yield_stress(&yield_stress, mat_id, density);
  yield_weakening_apply_temperature_to_yield_stress(&yield_stress, mat_id, density, u);

  return yield_stress;
}

/**
 * @brief Apply the yield stress to the deviatoric stress tensor.
 *
 * ### Something about yield criterion
 *
 * @param deviatoric_stress_tensor The deviatoric stress tensor to be modified.
 * @param yield_stress The yield stress.
 * @param density The density.
 * @param u The specific internal energy.
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_deviatoric_stress_tensor(
    struct sym_matrix *deviatoric_stress_tensor,
    const float yield_stress, const float density, const float u) {

  /* Calculate the J_2 invariant of the deviatoric stress tensor. */
  float J_2 = strength_compute_stress_tensor_J_2(*deviatoric_stress_tensor);

  /* ### Comment with name of yield criterion */
  float f = fminf((yield_stress * yield_stress) / (3.f * J_2), 1.f);

  /* Reduce elements of deviatoric stress tensor based on yield stress */
  deviatoric_stress_tensor->xx *= f;
  deviatoric_stress_tensor->yy *= f;
  deviatoric_stress_tensor->zz *= f;
  deviatoric_stress_tensor->xy *= f;
  deviatoric_stress_tensor->xz *= f;
  deviatoric_stress_tensor->yz *= f;
}

#endif /* SWIFT_YIELD_STRESS_BA95_H */
