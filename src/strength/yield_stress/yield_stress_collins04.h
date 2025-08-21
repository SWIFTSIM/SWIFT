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
#ifndef SWIFT_YIELD_STRESS_COLLINS04_H
#define SWIFT_YIELD_STRESS_COLLINS04_H

/**
 * @file strength/yield_stress/yield_stress_collins04.h
 * Collins+2004 pressure-dependent yield stress.
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
 * With the Collins+2004 yield stress method, damage does not act on the
 * deviatoric stress tensor. Damage is instead applied to the yield stress.
 *
 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, const float damage) {

  return deviatoric_stress_tensor;
}

/**
 * @brief Compute damaged yield stress.
 *
 * With the Collins+2004 yield stress method, damage acts to linearly transition
 * the yield stress between the fully intact and fully damaged yield stresses.
 *
 * @param yield_stress_fully_intact The yield stress for fully intact material.
 * @param yield_stress_fully_damaged The yield stress for fully damaged material.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact, const float yield_stress_fully_damaged, const float damage) {

  /* Compute damaged yield stress, Collins+2004 Eqn. ### */
  return (1.f - damage) * yield_stress_fully_intact + damage * yield_stress_fully_damaged;
}

/**
 * @brief Compute the yield stress of fully intact material.
 *
 * Collins+2004 pressure-dependent yield stress of fully intact solids.
 *
 * Method parameters needed in material parameter file:
 * YieldStress:
 *    mu_i: ###.
 *    Y_0: ### (Pa).
 *    Y_M: ### (Pa).
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

  /* Method parameters. */
  const float mu_i = material_mu_i(mat_id);
  const float Y_0 = material_Y_0(mat_id);
  const float Y_M = material_Y_M(mat_id);

  /* Compute yield stress of fully intact material, Collins+2004 Eqn. ### */
  // ### Should be able to decrease if negative pressures until yield_stress=0? miluphcuda does this so maybe not wrong?
  if (pressure > 0.f && Y_M != Y_0) {
    return Y_0 + mu_i * pressure / (1.f + (mu_i * pressure) / (Y_M - Y_0));
  } else {
    return Y_0;
  }
}

/**
 * @brief Compute the yield stress of fully damaged material.
 *
 * Collins+2004 pressure-dependent yield stress of fully damaged solids.
 *
 * Method parameters needed in material parameter file:
 * YieldStress:
 *    mu_i: ###.
 *    Y_0: ### (Pa).
 *    Y_M: ### (Pa).
 *
 * @param mat_id The material ID.
 * @param phase_state The phase ID.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const int mat_id, const int phase_state, const float pressure) {

  /* Return 0.f if the material is not solid. */
  if (phase_state != mat_phase_state_solid) {
    return 0.f;
  }

  /* Compute yield stress of fully damaged material, Collins+2004 Eqn. ### */
  // ### Maybe yield_stress_damaged also needs have a max value indep of
  // ### yield_stress_fully_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018
  if (pressure > 0.f) {
    const float mu_d = material_mu_d(mat_id);
    const float yield_stress_fully_intact = yield_compute_yield_stress_fully_intact(mat_id, phase_state, pressure);
    return fminf(mu_d, yield_stress_fully_intact);
  } else {
    return 0.f;
  }
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

  /* Calculate yield stresses of fully intact and fully damaged material. */
  const float yield_stress_fully_intact =
      yield_compute_yield_stress_fully_intact(mat_id, phase_state, pressure);
  const float yield_stress_fully_damaged =
      yield_compute_yield_stress_fully_damaged(mat_id, phase_state, pressure);

  /* Combine yield stresses based on how damaged the material is. */
  float yield_stress =
      yield_compute_damaged_yield_stress(yield_stress_fully_intact, yield_stress_fully_damaged, damage);

  /* Apply weakening to yield stress. */
  yield_weakening_apply_density_to_yield_stress(&yield_stress, mat_id, density);
  yield_weakening_apply_temperature_to_yield_stress(&yield_stress, mat_id, density, u);

  return yield_stress;
}

/**
 * @brief Apply the yield stress to a symmetric matrix.
 *
 * ### Something about yield criterion
 *
 * @param M The symmetric matrix to be modified.
 * @param deviatoric_stress_tensor The deviatoric stress tensor.
 * @param density The density.
 * @param u The specific internal energy.
 * @param yield_stress The yield stress.
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, struct sym_matrix deviatoric_stress_tensor,
    const float density, const float u, const float yield_stress) {

  /* Calculate the J_2 invariant of the deviatoric stress tensor. */
  float J_2 = strength_compute_deviatoric_sym_matrix_J_2(deviatoric_stress_tensor);

  /* ### Comment with name of yield criterion */
  float f = fminf(yield_stress / sqrtf(J_2), 1.f);

  /* Reduce elements of deviatoric stress tensor based on yield stress */
  M->xx *= f;
  M->yy *= f;
  M->zz *= f;
  M->xy *= f;
  M->xz *= f;
  M->yz *= f;
}

#endif /* SWIFT_YIELD_STRESS_COLLINS04_H */
