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
#ifndef SWIFT_DAMAGE_SHEAR_NONE_H
#define SWIFT_DAMAGE_SHEAR_NONE_H

/**
 * @file strength/damage/damage_shear/damage_shear_none.h
 * @brief No shear fracture models.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Get shear damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage(const struct part *restrict p) {

  return 0.f;
}

/**
 * @brief Get shear damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage_full(const struct xpart *restrict xp) {

  return 0.f;
}

/**
 * @brief Set drift-time shear damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage(struct part *restrict p, const float shear_damage) {}

/**
 * @brief Set kick-time shear damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage_full(struct xpart *restrict xp, const float shear_damage_full) {}

/**
 * @brief Compute the rate of shear damage accumulation.
 *
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param is_above_yield_criterion Whether the yield criterion is reached.s
 * @param strain_rate_tensor The strain rate tensor.
 * @param mat_id The material ID.
 * @param pressure The pressure.
 * @param shear_damage The shear damage.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const int is_above_yield_criterion, const struct sym_matrix strain_rate_tensor,
    const int mat_id, const float pressure, const float shear_damage) {

  /* Set the rate of shear damage accumulation to zero. */
  *shear_dD_dt = 0.f;
}

/**
 * @brief Steps shear damage by applying time-step to a shear_dD_dt.
 *
 * @param shear_damage The shear damage.
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt, const float dt_therm) {}

/**
 * @brief Evolves shear damage.
 *
 * Carries out all calculations required to step shear damage in time.
 *
 * @param shear_damage The shear damage.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param strain_rate_tensor The strain rate tensor.
 * @param mat_id The material ID.
 * @param pressure The pressure.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_evolve(
    float *shear_damage, const int is_above_yield_criterion, const struct sym_matrix strain_rate_tensor,
    const int mat_id, const float pressure, const float dt_therm) {}

#endif /* SWIFT_DAMAGE_SHEAR_NONE_H */
