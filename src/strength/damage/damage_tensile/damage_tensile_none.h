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
#ifndef SWIFT_DAMAGE_TENSILE_NONE_H
#define SWIFT_DAMAGE_TENSILE_NONE_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_none.h
 * @brief No tensile fracture models.
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Get tensile damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_tensile_damage(const struct part *restrict p) {

  return 0.f;
}

/**
 * @brief Get tensile damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_tensile_damage_full(const struct xpart *restrict xp) {

  return 0.f;
}

/**
 * @brief Set drift-time tensile damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_tensile_damage(struct part *restrict p, const float tensile_damage) {}

/**
 * @brief Set kick-time tensile damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_tensile_damage_full(struct xpart *restrict xp, const float tensile_damage_full) {}

/**
 * @brief Calculates the rate of cbrt(damage) accumulation due to tension.
 *
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param number_of_flaws The total number of flaws.
 * @param activation_thresholds The activation thresholds fo flaws.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,  const int number_of_flaws, const float activation_thresholds[40], // ### Change this length
    const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage) {}

/**
 * @brief Calculates the rate of damage accumulation due to tension
 *
 * @param tensile_dD_dt The rate of tensile damage accumulation.
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 * @param tensile_damage The tensile damage.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_dD_dt(
    float *tensile_dD_dt, struct part *restrict p,  const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage, const float tensile_damage) {}

/**
 * @brief Steps tensile damage by applying time-step to a tensile_cbrtD_dt.
 *
 * @param tensile_damage The tensile damage.
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param number_of_flaws The total number of flaws.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_apply_timestep_to_tensile_damage(
    float *tensile_damage, const float tensile_cbrtD_dt,
    const int number_of_activated_flaws, const int number_of_flaws, const float dt_therm) {}

/**
 * @brief Evolves tensile damage.
 *
 * @param tensile_damage The tensile damage.
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_evolve(
    float *tensile_damage, struct part *restrict p, const struct sym_matrix stress_tensor, const int mat_id,
    const float mass, const float density, const float damage, const float dt_therm) {}

#endif /* SWIFT_DAMAGE_TENSILE_NONE_H */
