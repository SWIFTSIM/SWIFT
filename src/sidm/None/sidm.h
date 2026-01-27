/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Katy Proctor (katy.proctor@fysik.su.se)
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
#ifndef SWIFT_NONE_SIDM_H
#define SWIFT_NONE_SIDM_H

/* Local includes */
#include "sidm_part.h"
#include "sidm_properties.h"

/**
 * @brief Initialises the si-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param sip The particle to act upon
 * @param sidm_properties The properties of the SIDM model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param scale_factor The current scale-factor (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 */
__attribute__((always_inline)) INLINE static void sidm_first_init_sipart(
    struct sipart *sip, const struct sidm_props *sidm_properties) {
  sip->time_bin = 0;
}

/**
 * @brief Prepares an si-particle for its interactions
 *
 * @param sip The particle to act upon
 */
__attribute__((always_inline)) INLINE static void sidm_init_sipart(
    struct sipart *sip) {}

/**
 * @brief Finishes the calculation of density on SIDM
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_end_density(
    struct sipart *sip) {}

/**
 * @brief Sets all particle fields to sensible values when the #sipart has 0
 * ngbs.
 *
 * @param sip The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void sidm_sipart_has_no_neighbours(
    struct sipart *restrict sip, const struct cosmology *cosmo) {}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * @param sip The #sipart.
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void sidm_predict_extra(
    struct sipart *restrict sip, float dt_drift) {}

/**
 * @brief Returns the comoving density of a particle
 *
 * @param sip The si-particle of interest
 */
__attribute__((always_inline)) INLINE static float sidm_get_comoving_density(
    const struct sipart *restrict sip) {

  return sip->rho;
}

/**
 * @brief Kick the additional variables
 *
 * @param sip The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void sidm_kick_extra(
    struct sipart *sip, float dt) {}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param sip The particle.
 */
__attribute__((always_inline)) INLINE static void sidm_reset_predicted_values(
    struct sipart *restrict sip) {}

/**
 * @brief Computes the time-step of a given SIDM particle.
 *
 * @param sip Pointer to the si-particle data.
 * @param sidm_properties Properties of the SIDM model.
 * @param with_cosmology Are we running with cosmological time integration.
 * @param cosmo The current cosmological model (used if running with
 * cosmology).
 * @param time The current time (used if running without cosmology).
 * @param time_base The time base.
 */
__attribute__((always_inline)) INLINE static float sidm_compute_timestep(
    const struct sipart *const sip, const struct sidm_props *sidm_properties,
    const int with_cosmology, const struct cosmology *cosmo, const double time,
    const double time_base) {
  return FLT_MAX;
}

/**
 * @brief Operations performed when an SIDM particle gets removed from the
 * simulation volume.
 *
 * @param sip The si-particle.
 * @param time The simulation time.
 */
__attribute__((always_inline)) INLINE static void sidm_remove_sipart(
    const struct sipart *sip, const double time) {}

#endif /* SWIFT_NONE_SIDM_H */
