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

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure and applying damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_stress_tensor(struct part *restrict p,
    struct sym_matrix stress_tensor, struct sym_matrix deviatoric_stress_tensor, const float pressure, const float damage) {}

__attribute__((always_inline)) INLINE static void damage_timestep(
    const struct part *restrict p, float dt_cfl) {}

__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values_extra_damage(
    struct part *restrict p, const struct xpart *restrict xp) {}

/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage(
    struct part *restrict p, float *tensile_damage, float *shear_damage,
    float *damage, struct sym_matrix deviatoric_stress_tensor, const float yield_stress, 
    const float density, const float u, const float dt_therm) {}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_evolve_damage(
    struct part *restrict p, const float yield_stress, const float density, const float u, const float dt_therm) {}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_kick_evolve_damage(
    struct part *restrict p, struct xpart *restrict xp, const float yield_stress, const float density, const float u, const float dt_therm) {}

/**
 * @brief Calculate time derivative of damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void update_dD_dt(
    struct part *restrict p, const int phase_state, const float density, const float u) {}

__attribute__((always_inline)) INLINE static void hydro_first_init_part_damage(
    struct part *restrict p, struct xpart *restrict xp) {}


#endif /* SWIFT_DAMAGE_NONE_H */
