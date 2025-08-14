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
#ifndef SWIFT_DAMAGE_TENSILE_NONE_H
#define SWIFT_DAMAGE_TENSILE_NONE_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_none.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Calculates the rate of cbrt(damage) accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,  const int number_of_flaws, const float activation_thresholds[40], // ### Change this length
    const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage) {}

/**
 * @brief Calculates the rate of damage accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_dD_dt(
    struct part *restrict p, float *tensile_dD_dt, const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage, const float tensile_damage) {}

/**
 * @brief Steps particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_tensile_apply_timestep_to_tensile_damage(
    float *tensile_damage, const float tensile_cbrtD_dt, 
    const int number_of_activated_flaws, const int number_of_flaws, const float dt_therm) {}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_tensile_evolve(
    struct part *restrict p, float *tensile_damage, const struct sym_matrix stress_tensor, const int mat_id, 
    const float mass, const float density, const float damage, const float dt_therm) {}

#endif /* SWIFT_DAMAGE_TENSILE_NONE_H */