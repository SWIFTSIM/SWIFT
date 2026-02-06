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
#ifndef SWIFT_REMIX_STRENGTH_NONE_H
#define SWIFT_REMIX_STRENGTH_NONE_H

/**
 * @file REMIX/strength/none/hydro_strength.h
 * @brief REMIX implementation of SPH with no material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength.h"

/**
 * @brief Updates the hydro (+ strength) time-step based on strength methods.
 *
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void hydro_compute_timestep_strength(
    float *dt_cfl, const struct part *restrict p, const struct hydro_props *restrict hydro_properties) {}

/**
 * @brief Updates the max wave speed based on strength methods.
 *
 * @param wave_speed The wave speed to be updated.
 * @param p The particle of interest.
 * @param soundspeed The sound speed.
 * @param density The sound density.
 */
__attribute__((always_inline)) INLINE static void
hydro_compute_max_wave_speed_strength(float *wave_speed, const struct part *restrict p, const float soundspeed, const float density) {}

/**
 * @brief Prepares extra strength parameters for a particle for the density
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_init_part_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the density calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_density_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the gradient
 * calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_gradient_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the gradient calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_gradient_strength(struct part *restrict p) {}

/**
 * @brief Prepares extra strength parameters for a particle for the force
 * calculation.
 *
 * @param p The particle to act upon
 * @param density The density
 * @param u The specific internal energy
 */
__attribute__((always_inline)) INLINE static void
hydro_prepare_force_strength(struct part *restrict p,
                                   const float density, const float u) {}

/**
 * @brief Resets strength time derivative fields in preparation
  * for the sums taking  place in the various force tasks.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_reset_acceleration_strength(struct part *restrict p) {}

/**
 * @brief Finishes extra strength parts of the force calculation.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_end_force_strength(struct part *restrict p) {}

/**
 * @brief Sets the values of additional particle strength properties at a
 * kick time
 *
 * @param p The particle.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values_strength(
    struct part *restrict p, const struct xpart *restrict xp) {}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At beginning of hydro function, before hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_strength_beginning(
    struct part *restrict p, const float dt_therm) {}

/**
 * @brief Predict additional particle strength properties forward in time when
 * drifting. At end of hydro function, after hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void hydro_predict_strength_end(
    struct part *restrict p, const float dt_therm) {}

/**
 * @brief Kick the additional particle strength properties.
 * At beginning of hydro function, before hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_strength_beginning(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm) {}

/**
 * @brief Kick the additional particle strength properties
 * At end of hydro function, after hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void hydro_kick_strength_end(
    struct part *restrict p, struct xpart *restrict xp, float dt_therm) {}

/**
 * @brief Initialises the strength properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle to act upon
 * @param xp The extended particle data to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_first_init_part_strength(
    struct part *restrict p, struct xpart *restrict xp) {}

#endif /* SWIFT_REMIX_STRENGTH_NONE_H */
