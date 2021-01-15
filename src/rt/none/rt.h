/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_NONE_H
#define SWIFT_RT_NONE_H

#include "rt_thermochemistry.h"

/**
 * @file src/rt/none/rt.h
 * @brief Main header file for no radiative transfer scheme.
 */

/**
 * @brief Initialisation of the RT density loop related particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {}

/**
 * @brief Reset of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {}

/**
 * @brief First initialisation of the RT extra hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p) {}

/**
 * @brief Initialisation of the RT density loop related particle data.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {}

/**
 * @brief First initialisation of the RT extra star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Update the photon number of a particle, i.e. compute
 *        E^{n+1} = E^n + dt * dE_* / dt
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p) {}

/**
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is initialized
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param time current system time
 * @param star_age age of the star *at the end of the step*
 * @param dt star time step
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp, double time,
                                 double star_age, double dt) {}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_finalise_gradient(
    struct part* restrict p) {}

/**
 * @brief finishes up the transport step by actually doing the time integration
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p) {}

/**
 * @brief Wraps around rt_do_thermochemistry function
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_tchem(
    struct part* restrict p) {

  rt_do_thermochemistry(p);
}

#endif /* SWIFT_RT_NONE_H */
