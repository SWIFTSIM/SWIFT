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
#ifndef SWIFT_DAMAGE_NONE_H
#define SWIFT_DAMAGE_NONE_H

/**
 * @file strength/damage/damage_none.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

__attribute__((always_inline)) INLINE static float strength_get_damage(struct part *restrict p) {

  return 0.f;
}

__attribute__((always_inline)) INLINE static float strength_get_damage_full(const struct xpart *restrict xp) {

  return 0.f;
}

__attribute__((always_inline)) INLINE static void strength_compute_timestep_damage(
    const struct part *restrict p, float *dt_cfl) {}

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure and applying damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_compute_stress_tensor(
    struct sym_matrix *stress_tensor, const struct sym_matrix damaged_deviatoric_stress_tensor, const float pressure, const float damage) {}

__attribute__((always_inline)) INLINE static void damage_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp) {}

/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_evolve(
    struct part *restrict p, float *tensile_damage, float *shear_damage, float *damage,
    const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_predict_evolve(
    struct part *restrict p, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_kick_evolve(
    struct part *restrict p, struct xpart *restrict xp, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {}

/**
 * @brief Calculate time derivative of damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_compute_dD_dt(
    struct part *restrict p, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress) {}

__attribute__((always_inline)) INLINE static void damage_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {}


#endif /* SWIFT_DAMAGE_NONE_H */
