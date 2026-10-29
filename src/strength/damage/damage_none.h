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
#ifndef SWIFT_DAMAGE_NONE_H
#define SWIFT_DAMAGE_NONE_H

/**
 * @file strength/damage/damage_none.h
 * @brief No fracture models.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Get damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static struct strength_damage strength_get_damage(const struct part *restrict p) {

  return (struct strength_damage){0};
}

/**
 * @brief Get damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static struct strength_damage strength_get_damage_full(const struct xpart *restrict xp) {

  return (struct strength_damage){0};
}

/**
 * @brief Sets the damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

/**
 * @brief Computes the damage time-step of a given particle
 *
 * Calculates a time-step based on the particle's rate of damage accumulation.
 * If this time-step is smaller than dt_cfl, dt_cfl gets overwritten to this
 * damage time-step.
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void strength_compute_timestep_damage(
    float *dt_cfl, const struct part *restrict p,
    const struct strength_props *strength_props) {}

/**
 * @brief Sets the values of additional particle damage properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_damage(
    struct part *restrict p, const struct xpart *restrict xp) {}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle of interest.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_predict_evolve(
    struct part *restrict p, const float dt_therm) {}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_kick_evolve(
    const struct part *restrict p, struct xpart *restrict xp,
    const float dt_therm) {}

/**
 * @brief Computes the rates of damage accumulation.
 *
 * @param p The particle of interest.
 * @param s The strength state of the particle, from its drifted values.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_compute_rates(
    struct part *restrict p, struct strength_state *s,
    const struct mat_params *mp) {}

/**
 * @brief Initialises the damage properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_first_init_part_damage(
    struct part *restrict p, struct xpart *restrict xp) {}


#endif /* SWIFT_DAMAGE_NONE_H */
