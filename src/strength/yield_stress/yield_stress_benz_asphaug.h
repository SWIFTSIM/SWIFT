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
#ifndef SWIFT_YIELD_STRESS_BENZ_ASPHAUG_H
#define SWIFT_YIELD_STRESS_BENZ_ASPHAUG_H

/**
 * @file strength/yield_stress/yield_stress_benz_asphaug.h
 * Benz&Asphaug95 constant yield stress.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Sets the yield stress properties of the strength state from the particle.
 *
 * Nothing extra is needed in the state for this method.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void yield_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

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
 * material acts as a fluid. See Benz&Asphaug95 Eqn. 13.

 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, struct strength_state *s,
    const struct mat_params *mp) {

  /* Compute elements of damaged deviatoric stress tensor. */
  struct sym_matrix damaged_deviatoric_stress_tensor = deviatoric_stress_tensor;
  damaged_deviatoric_stress_tensor.xx *= (1.f - s->damage.total);
  damaged_deviatoric_stress_tensor.yy *= (1.f - s->damage.total);
  damaged_deviatoric_stress_tensor.zz *= (1.f - s->damage.total);
  damaged_deviatoric_stress_tensor.xy *= (1.f - s->damage.total);
  damaged_deviatoric_stress_tensor.xz *= (1.f - s->damage.total);
  damaged_deviatoric_stress_tensor.yz *= (1.f - s->damage.total);

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
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact,
                                   const float yield_stress_fully_damaged,
                                   struct strength_state *s,
                                   const struct mat_params *mp) {

  return yield_stress_fully_intact;
}

/**
 * @brief Compute the yield stress of fully intact material.
 *
 * Benz&Asphaug95 constant yield stress of solids.
 *
 * Method parameters needed in material parameter file:
 * YieldStressBenzAsphaug:
 *    Y_0: Constant yield stress (Pa).
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_intact(
    struct strength_state *s, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (s->phase != mat_phase_solid) {
    return 0.f;
  }

  /* Return constant yield stress, Y_0. */
  return mp->yield_stress_benz_asphaug_Y_0;
}

/**
 * @brief Compute the yield stress of fully damaged material.
 *
 * With the Benz&Asphaug95 yield stress method, damage does not affect the yield
 * stress. Damage is instead applied to the deviatoric stress tensor.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    struct strength_state *s, const struct mat_params *mp) {

  /* Damage does not affect the yield stress. */
  return yield_compute_yield_stress_fully_intact(s, mp);
}

/**
 * @brief Compute the yield stress.
 *
 * Calculates the yield stress that combines all relevant methods e.g. damage
 * and weakening.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress(
    struct strength_state *s, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (s->phase != mat_phase_solid) {
    return 0.f;
  }

  /* Get constant yield stress. */
  float yield_stress = yield_compute_yield_stress_fully_intact(s, mp);

  /* Apply weakening to yield stress. */
  yield_weakening_apply_to_yield_stress(&yield_stress, s, mp);

  return yield_stress;
}

/**
 * @brief Check the yield criterion.
 *
 * Yield criterion sqrt(3 J_2) >= Y, evaluated as 3 J_2 >= Y^2 (Y is positive).
 * Returns 1 if the stress is on or beyond the yield surface, 0 if it is inside.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static int yield_check_yield_criterion(
    struct strength_state *s, const struct mat_params *mp) {

  const float yield_stress = strength_state_get_yield_stress(s, mp);

  return 3.f * strength_state_get_J_2(s, mp) >= yield_stress * yield_stress;
}

/**
 * @brief Apply the yield stress to a symmetric matrix.
 *
 * ### Something about yield criterion
 *
 * @param M The symmetric matrix to be modified.
 * @param deviatoric_stress_tensor The deviatoric stress tensor to check
 * against the yield stress.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, const struct sym_matrix deviatoric_stress_tensor,
    struct strength_state *s, const struct mat_params *mp) {

  const float yield_stress = strength_state_get_yield_stress(s, mp);

  /* ### Comment with name of yield criterion */
  float f = 1.f;
  const float J_2 = strength_compute_deviatoric_sym_matrix_J_2(deviatoric_stress_tensor);
  if (J_2 > 0.f) {
      f = fminf((yield_stress * yield_stress) / (3.f * J_2), 1.f);
  }

  /* Reduce elements of symmetric matrix based on yield stress */
  M->xx *= f;
  M->yy *= f;
  M->zz *= f;
  M->xy *= f;
  M->xz *= f;
  M->yz *= f;
}

#endif /* SWIFT_YIELD_STRESS_BENZ_ASPHAUG_H */
