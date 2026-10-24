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
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, const float damage,
    const struct mat_params *mp) {

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
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact,
                                   const float yield_stress_fully_damaged,
                                   const float damage, const struct mat_params *mp) {

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
 * @param phase The phase ID.
 * @param pressure The pressure.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_intact(
    const int phase, const float pressure, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (phase != mat_phase_solid) {
    return 0.f;
  }

  /* Method parameters. */
  const float mu_i = mp->mu_i;
  const float Y_0 = mp->Y_0;
  const float Y_M = mp->Y_M;

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
 * @param phase The phase ID.
 * @param pressure The pressure.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const int phase, const float pressure, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (phase != mat_phase_solid) {
    return 0.f;
  }

  /* Compute yield stress of fully damaged material, Collins+2004 Eqn. ### */
  // ### Maybe yield_stress_damaged also needs have a max value indep of
  // ### yield_stress_fully_intact? See e.g. Güldemeister et al. 2015; Winkler et al. 2018
  if (pressure > 0.f) {
    const float mu_d = mp->mu_d;
    const float yield_stress_fully_intact = yield_compute_yield_stress_fully_intact(phase, pressure, mp);
    return fminf(mu_d * pressure, yield_stress_fully_intact);
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
 * @param phase The phase ID.
 * @param density The density.
 * @param pressure The pressure.
 * @param u The specific internal energy.
 * @param damage The damage.
 * @param mat_id The material ID, for the temperature from the EoS.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress(
    const int phase, const float density, const float pressure, const float u,
    const float damage, const int mat_id, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (phase != mat_phase_solid) {
    return 0.f;
  }

  /* Calculate yield stresses of fully intact and fully damaged material. */
  const float yield_stress_fully_intact =
      yield_compute_yield_stress_fully_intact(phase, pressure, mp);
  const float yield_stress_fully_damaged =
      yield_compute_yield_stress_fully_damaged(phase, pressure, mp);

  /* Combine yield stresses based on how damaged the material is. */
  float yield_stress =
      yield_compute_damaged_yield_stress(yield_stress_fully_intact, yield_stress_fully_damaged, damage, mp);

  /* Apply weakening to yield stress. */
  yield_weakening_apply_density_to_yield_stress(&yield_stress, density, mp);
  yield_weakening_apply_temperature_to_yield_stress(&yield_stress, density, u, mat_id, mp);

  return yield_stress;
}

/**
 * @brief Check the yield criterion.
 *
 * Yield sqrt(J_2) >= Y, evaluated as J_2 >= Y^2 (Y is positive).
 * Returns 1 if the stress is on or beyond the yield surface, 0 if it is inside.
 *
 * @param deviatoric_stress_tensor The deviatoric stress tensor.
 * @param yield_stress The yield stress.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static int yield_check_yield_criterion(
    struct sym_matrix deviatoric_stress_tensor, const float yield_stress,
    const struct mat_params *mp) {

  const float J_2 =
      strength_compute_deviatoric_sym_matrix_J_2(deviatoric_stress_tensor);

  return J_2 >= yield_stress * yield_stress;
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
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, struct sym_matrix deviatoric_stress_tensor,
    const float density, const float u, const float yield_stress,
    const struct mat_params *mp) {

  /* Calculate the J_2 invariant of the deviatoric stress tensor. */
  float J_2 = strength_compute_deviatoric_sym_matrix_J_2(deviatoric_stress_tensor);

  /* ### Comment with name of yield criterion */
  float f = 1.f;
  if (J_2 > 0.f) {
      f = fminf(yield_stress / sqrtf(J_2), 1.f);
  }

  /* Reduce elements of deviatoric stress tensor based on yield stress */
  M->xx *= f;
  M->yy *= f;
  M->zz *= f;
  M->xy *= f;
  M->xz *= f;
  M->yz *= f;
}

#endif /* SWIFT_YIELD_STRESS_COLLINS04_H */
