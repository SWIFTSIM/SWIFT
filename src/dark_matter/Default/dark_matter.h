/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Camila Correa (camila.correa@cea.fr)
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
#ifndef SWIFT_DEFAULT_DARK_MATTER_H
#define SWIFT_DEFAULT_DARK_MATTER_H

/* Config parameters. */
#include <config.h>
#include "approx_math.h"
#include "dark_matter_part.h"
#include "dark_matter_iact.h"

/**
 * @brief Prepares a dm-particle for its interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void dark_matter_init_dmpart(struct dmpart* dmp) {}


/**
 * @brief Prepares a dm-particle for its SIDM interactions
 *
 * @param gp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_velocities(struct dmpart* dmp, float dt_drift) {}


/**
 * @brief Initialises the dark matter particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param gp The particle to act upon.
 * @param sidm_properties Properties of the self-interacting dark matter model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_first_init_dmpart(
     struct dmpart* dmp, const struct sidm_props* sidm_props) {}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param gp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void dark_matter_predict_extra(
               struct dmpart* restrict dmp, float dt_drift) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param gp The particle.
 */
__attribute__((always_inline)) INLINE static void dark_matter_reset_predicted_values(
                struct dmpart* restrict dmp) {}


/**
 * @brief Finishes the calculation of density on dark matter particles
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_end_density(
    struct dmpart* dmp, const struct cosmology* cosmo, const struct sidm_props* sidm_props,
    const double dt) {}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param sp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void dark_matter_part_has_no_neighbours(
    struct dmpart* restrict dmp, const struct cosmology* cosmo) {}



/**
 * @brief Updates #dmparts velocities
 *
 * @param dmp #dmpart
 *
 */
__attribute__((always_inline)) INLINE static void sidm_kick_to_dmpart(struct dmpart *restrict dmp) {}


/**
 * @brief Computes the dark matter time-step of a given particle
 *
 * This function returns the time-step of a particle given its DM-dynamical
 * state?
 *
 * @param p Pointer to the particle data
 * @param sidm_properties The SIDM parameters
 */
__attribute__((always_inline)) INLINE static float dark_matter_compute_timestep(
    const struct dmpart *restrict dmp, const struct sidm_props* sidm_props,
    const struct cosmology* cosmo) {

  return FLT_MAX;
}

/**
 * @brief Kick the additional variables
 *
 * @param dmp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void dark_matter_kick_extra(struct dmpart* dmp, float dt) {}


#endif

