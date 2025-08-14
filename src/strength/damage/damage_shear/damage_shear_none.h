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
#ifndef SWIFT_DAMAGE_SHEAR_NONE_H
#define SWIFT_DAMAGE_SHEAR_NONE_H

/**
 * @file strength/damage/damage_shear/damage_shear_none.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

__attribute__((always_inline)) INLINE static float damage_get_shear_damage(const struct part *restrict p) {

  return 0.f;
}

__attribute__((always_inline)) INLINE static float damage_get_shear_damage_full(const struct xpart *restrict xp) {

  return 0.f;
}

__attribute__((always_inline)) INLINE static void damage_set_shear_damage(struct part *restrict p, const float shear_damage) {}

__attribute__((always_inline)) INLINE static void damage_set_shear_damage_full(struct xpart *restrict xp, const float shear_damage_full) {}

/**
 * @brief Calculates the rate of damage accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const struct sym_matrix deviatoric_stress_tensor, const int mat_id, const float density, const float u, const float yield_stress, const float shear_damage) {}

/**
 * @brief Steps particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt, const float dt_therm) {}

/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_shear_evolve(
    struct part *restrict p, float *shear_damage, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float density, const float u, const float yield_stress, const float dt_therm) {}

#endif /* SWIFT_DAMAGE_SHEAR_NONE_H */