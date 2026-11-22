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
 * @brief No yield stress method. There is no yield surface, so the yield
 * stress isn't computed or applied.
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
 * Damage scales the deviatoric stress tensor and negative pressures by (1 - D),
 * so that fully damaged material acts as a fluid. See Benz&Asphaug95 Eqn. 13.
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float yield_damage_factor(
    struct strength_state *s, const struct mat_params *mp) {

  return 1.f - s->damage.total;
}

/**
 * @brief Apply the yield stress to a symmetric matrix.
 *
 * Empty function when configuring without yield stress method.
 *
 * @param M The symmetric matrix to be modified.
 * @param deviatoric_stress_tensor The deviatoric stress tensor to check
 * against the yield stress.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return 0, since there is no yield surface.
 */
__attribute__((always_inline)) INLINE static int
yield_apply_yield_stress_to_sym_matrix(
    struct sym_matrix *M, const struct sym_matrix deviatoric_stress_tensor,
    struct strength_state *s, const struct mat_params *mp) {

  return 0;
}

#endif /* SWIFT_YIELD_STRESS_NONE_H */
