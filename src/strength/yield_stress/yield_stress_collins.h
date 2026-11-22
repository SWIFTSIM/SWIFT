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
#ifndef SWIFT_YIELD_STRESS_COLLINS_H
#define SWIFT_YIELD_STRESS_COLLINS_H

/**
 * @file strength/yield_stress/yield_stress_collins.h
 * Collins+2004 pressure-dependent yield stress.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Sets the yield stress properties of the strength state from the
 * particle.
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
 * @brief Returns the damage factor by which damage scales the stress and
 * pressure.
 *
 * Damage does not scale the stress. It is instead applied to the yield stress,
 * and tension is limited by the minimum pressure.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_damage_factor(
    struct strength_state *s, const struct mat_params *mp) {

  return 1.f;
}

/**
 * @brief Compute damaged yield stress.
 *
 * With the Collins+2004 yield stress method, damage acts to linearly transition
 * the yield stress between the fully intact and fully damaged yield stresses.
 *
 * @param yield_stress_fully_intact The yield stress for fully intact material.
 * @param yield_stress_fully_damaged The yield stress for fully damaged
 * material.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_damaged_yield_stress(const float yield_stress_fully_intact,
                                   const float yield_stress_fully_damaged,
                                   struct strength_state *s,
                                   const struct mat_params *mp) {

  return (1.f - s->damage.total) * yield_stress_fully_intact +
         s->damage.total * yield_stress_fully_damaged;
}

/**
 * @brief Compute the yield stress of fully intact material.
 *
 * Pressure-dependent yield stress, iSALE-Dellen manual Eqn. 4.68. This also
 * applies in tension, down to Y=0.
 *
 * Method parameters needed in material parameter file:
 * YieldStressCollins:
 *    Y_0i: Cohesion (yield stress at zero pressure) of intact material.
 *    Y_Mi: Limiting yield stress at high pressure of intact material.
 *    mu_i: Coefficient of internal friction of intact material.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_yield_stress_fully_intact(struct strength_state *s,
                                        const struct mat_params *mp) {

  /* Return 0.f if the material is not solid. */
  if (strength_state_get_phase(s, mp) != mat_phase_solid) {
    return 0.f;
  }

  /* Method parameters. */
  const float mu_i = mp->yield_stress_collins_mu_i;
  const float Y_0i = mp->yield_stress_collins_Y_0i;
  const float Y_Mi = mp->yield_stress_collins_Y_Mi;
  const float pressure = strength_state_get_pressure(s, mp);

  /* Constant yield stress case. */
  if ((Y_Mi <= Y_0i) || (mu_i == 0.f)) {
    return Y_0i;
  }

  /* The yield stress reaches zero at pressure_zero in tension, and is zero
   * below it. */
  const float pressure_zero = -Y_0i * (Y_Mi - Y_0i) / (mu_i * Y_Mi);
  if (pressure <= pressure_zero) {
    return 0.f;
  }

  return Y_0i + mu_i * pressure / (1.f + mu_i * pressure / (Y_Mi - Y_0i));
}

/**
 * @brief Compute the yield stress of fully damaged material.
 *
 * Pressure-dependent yield stress, iSALE-Dellen manual Eqn. 4.67, capped at
 * the intact yield stress. It is left negative in tension, so that the
 * combined yield stress interpolates smoothly between the intact and damaged
 * profiles. The combined value is then limited to zero or above.
 *
 * Method parameters needed in material parameter file:
 * YieldStressCollins:
 *    Y_0d: Cohesion (yield stress at zero pressure) of damaged material.
 *    Y_Md: Limiting yield stress at high pressure of damaged material.
 *    mu_d: Coefficient of internal friction of damaged material.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
yield_compute_yield_stress_fully_damaged(struct strength_state *s,
                                         const struct mat_params *mp) {

  /* Method parameters. */
  const float mu_d = mp->yield_stress_collins_mu_d;
  const float Y_0d = mp->yield_stress_collins_Y_0d;
  const float Y_Md = mp->yield_stress_collins_Y_Md;
  const float pressure = strength_state_get_pressure(s, mp);

  /* Damaged material cannot be stronger than intact material. */
  const float yield_stress_fully_intact =
      yield_compute_yield_stress_fully_intact(s, mp);

  return fminf(fminf(Y_0d + mu_d * pressure, Y_Md), yield_stress_fully_intact);
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
  if (strength_state_get_phase(s, mp) != mat_phase_solid) {
    return 0.f;
  }

  /* Calculate yield stresses of fully intact and fully damaged material. */
  const float yield_stress_fully_intact =
      yield_compute_yield_stress_fully_intact(s, mp);
  const float yield_stress_fully_damaged =
      yield_compute_yield_stress_fully_damaged(s, mp);

  /* Combine yield stresses based on how damaged the material is. */
  float yield_stress =
      fmaxf(yield_compute_damaged_yield_stress(
                yield_stress_fully_intact, yield_stress_fully_damaged, s, mp),
            0.f);

  /* Apply weakening to yield stress. */
  yield_weakening_apply_to_yield_stress(&yield_stress, s, mp);

  return yield_stress;
}

/**
 * @brief Apply the yield stress to a symmetric matrix.
 *
 * Scales the matrix back by Y / sqrt(J_2) when sqrt(J_2) exceeds the yield
 * stress.
 *
 * @param M The symmetric matrix to be modified.
 * @param deviatoric_stress_tensor The deviatoric stress tensor to check
 * against the yield stress.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return 1 if the yield stress limited the matrix, i.e. the particle is
 * yielding, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, const struct sym_matrix deviatoric_stress_tensor,
    struct strength_state *s, const struct mat_params *mp) {

  const float yield_stress = strength_state_get_yield_stress(s, mp);

  /* Yield criterion sqrt(J_2) >= Y. */
  float f = 1.f;
  const float J_2 =
      strength_compute_deviatoric_sym_matrix_J_2(deviatoric_stress_tensor);
  if (J_2 > 0.f) {
    f = fminf(yield_stress / sqrtf(J_2), 1.f);
  }

  /* If the yield stress is not limiting, early return 0. */
  if (f >= 1.f) {
    return 0;
  }

  /* Reduce elements of symmetric matrix based on yield stress */
  M->xx *= f;
  M->yy *= f;
  M->zz *= f;
  M->xy *= f;
  M->xz *= f;
  M->yz *= f;

  return 1;
}

#endif /* SWIFT_YIELD_STRESS_COLLINS_H */
