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
#ifndef SWIFT_ARTIFICIAL_STRESS_NONE_H
#define SWIFT_ARTIFICIAL_STRESS_NONE_H

/**
 * @file strength/artificial_stress/artificial_stress_none.h
 * @brief No artificial stress method.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Sets the artificial stress properties of the strength state from the
 * particle.
 *
 * Nothing extra is needed in the state for this method.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void artif_stress_set_state(
    struct strength_state *s, const struct part *restrict p) {}

/**
 * @brief Apply artificial stress to pairwise stress tensors.
 *
 * Empty function when configuring without artificial stress.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interaction with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interaction with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
 * @param mp_i The material parameters of the first particle.
 * @param mp_j The material parameters of the second particle.
 */
__attribute__((always_inline)) INLINE static void artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r,
    const struct mat_params *mp_i, const struct mat_params *mp_j) {}

#endif /* SWIFT_ARTIFICIAL_STRESS_NONE_H */
