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
 * @param damage_rates The rates of damage accumulation.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    struct strength_damage_rates *damage_rates, struct strength_state *s,
    const struct mat_params *mp) {}

/**
 * @brief Evolves tensile damage.
 *
 * @param d The damage of the particle.
 * @param damage_rates The rates of damage accumulation.
 * @param dt_therm The time-step duration.
 *
 * @return The tensile damage, which is always zero.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_evolve(
    struct strength_damage *d, const struct strength_damage_rates *damage_rates,
    const float dt_therm) {

  return 0.f;
}

/**
 * @brief Limits the time-step by the growth of cbrt(tensile damage).
 *
 * @param dt_cfl The hydro + strength time-step.
 * @param damage_rates The rates of damage accumulation.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_timestep(
    float *dt_cfl, const struct strength_damage_rates *damage_rates,
    const struct strength_props *strength_props) {}

#endif /* SWIFT_DAMAGE_TENSILE_NONE_H */
