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
#include "math.h"

/**
 * @brief Sets the tensile damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

/**
 * @brief Checks the flaws of a particle read from the ICs.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_first_init_part(
    const struct part *restrict p) {}

/**
 * @brief Calculates the rate of cbrt(damage) accumulation due to tension.
 *
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,
    struct strength_state *s, const struct mat_params *mp) {

  /* Set the rate of tensile cbrt(damage) accumulation and the number of
   * activated flaws to zero. */
  *tensile_cbrtD_dt = 0.f;
  *number_of_activated_flaws = 0;
}

/**
 * @brief Steps tensile damage by applying time-step to a tensile_cbrtD_dt.
 *
 * @param tensile_damage The tensile damage.
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_apply_timestep_to_tensile_damage(
    float *tensile_damage, const float tensile_cbrtD_dt,
    const int number_of_activated_flaws, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {}

/**
 * @brief Evolves tensile damage.
 *
 * @param d The damage of the particle.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_evolve(
    struct strength_damage *d, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  return 0.f;
}

/**
 * @brief Computes the rate of tensile damage accumulation used to set the
 * damage time-step.
 *
 * @param d The damage of the particle.
 * @param damage_timestep_factor The growth in damage allowed per time-step.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return The rate of tensile damage accumulation, which is always zero.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_compute_timestep_dD_dt(
    const struct strength_damage *d, const float damage_timestep_factor,
    struct strength_state *s, const struct mat_params *mp) {

  return 0.f;
}

#endif /* SWIFT_DAMAGE_TENSILE_NONE_H */
