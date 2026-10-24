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
#ifndef SWIFT_YIELD_STRESS_NONE_H
#define SWIFT_YIELD_STRESS_NONE_H

/**
 * @file strength/yield_stress/yield_stress_none.h
 * @brief No yield stress method. If a yield stress is required for other
 * methods, functions return a yield stress of FLT_MAX.
 */

#include <float.h>
#include "const.h"
#include "math.h"

/**
 * @brief Sets the yield stress properties of the strength state from the particle.
 *
 * Nothing extra is needed in the state for this method.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void yield_set_state(
    struct strength_state *s, const struct part *restrict p) {}

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
 * With no yield stress method, damage acts to reduce all elements of the
 * deviatoric stress tensor, meaning that a fully damaged material acts as a
 * fluid. See Benz&Asphaug95 Eqn. 13.
 *
 * @param deviatoric_stress_tensor (sym_matrix) The deviatoric stress tensor.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static struct sym_matrix yield_compute_damaged_deviatoric_stress_tensor(
    struct sym_matrix deviatoric_stress_tensor, const struct strength_state *s,
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
 * With no yield stress method, damage does not affect the yield stress. Damage
 * is instead applied to the deviatoric stress tensor.
 *
 * @param yield_stress_fully_intact The yield stress for fully intact material.
 * @param yield_stress_fully_damaged The yield stress for fully damaged material.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact,
                                   const float yield_stress_fully_damaged,
                                   const struct strength_state *s,
                                   const struct mat_params *mp) {

  return yield_stress_fully_intact;
}

/**
 * @brief Compute the yield stress of fully intact material.
 *
 * With no yield stress method, the yield stress of a fully intact solid is set
 * to FLT_MAX.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_intact(
    const struct strength_state *s, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (s->phase != mat_phase_solid) {
    return 0.f;
  }

  /* Return constant yield stress, FLT_MAX. */
  return FLT_MAX;
}

/**
 * @brief Compute the yield stress of fully damaged material.
 *
 * With no yield stress method, damage does not affect the yield stress. Damage
 * is instead applied to the deviatoric stress tensor.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_compute_yield_stress_fully_damaged(
    const struct strength_state *s, const struct mat_params *mp) {

  /* Damage does not affect the yield stress i.e. it is still FLT_MAX. */
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
    const struct strength_state *s, const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (s->phase != mat_phase_solid) {
    return 0.f;
  }

  /* Get constant yield stress, FLT_MAX. */
  float yield_stress = yield_compute_yield_stress_fully_intact(s, mp);

  return yield_stress;
}

/**
 * @brief Check the yield criterion.
 *
 * With no yield stress method, there is no yield surface.
 *
 * @param yield_stress The yield stress.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static int yield_check_yield_criterion(
    const float yield_stress, const struct strength_state *s,
    const struct mat_params *mp) {

  return 0;
}

/**
 * @brief Apply the yield stress to a symmetric matrix.
 *
 * Empty function when configuring without yield stress method.
 *
 * @param M The symmetric matrix to be modified.
 * @param yield_stress The yield stress.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, const float yield_stress,
    const struct strength_state *s, const struct mat_params *mp) {}

#endif /* SWIFT_YIELD_STRESS_NONE_H */
