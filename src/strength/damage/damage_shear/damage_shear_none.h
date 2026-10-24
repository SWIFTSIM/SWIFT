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
 * @brief Sets the shear damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_shear_set_state(
    struct strength_state *s, const struct part *restrict p) {}

/**
 * @brief Compute the rate of shear damage accumulation.
 *
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param shear_damage The shear damage.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const float shear_damage,
    const struct strength_state *s, const struct mat_params *mp) {

  /* Set the rate of shear damage accumulation to zero. */
  *shear_dD_dt = 0.f;
}

/**
 * @brief Steps shear damage by applying time-step to a shear_dD_dt.
 *
 * @param shear_damage The shear damage.
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt,
    const struct strength_state *s, const struct mat_params *mp,
    const float dt_therm) {}

/**
 * @brief Evolves shear damage.
 *
 * Carries out all calculations required to step shear damage in time.
 *
 * @param d The damage of the particle.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 *
 * @return The shear damage, which is always zero.
 */
__attribute__((always_inline)) INLINE static float damage_shear_evolve(
    struct strength_damage *d, const struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  return 0.f;
}

/**
 * @brief Computes the rate of shear damage accumulation used to set the
 * damage time-step.
 *
 * @param d The damage of the particle.
 * @param damage_timestep_factor The growth in damage allowed per time-step.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return The rate of shear damage accumulation, which is always zero.
 */
__attribute__((always_inline)) INLINE static float damage_shear_compute_timestep_dD_dt(
    const struct strength_damage *d, const float damage_timestep_factor,
    const struct strength_state *s, const struct mat_params *mp) {

  return 0.f;
}

#endif /* SWIFT_DAMAGE_SHEAR_NONE_H */
