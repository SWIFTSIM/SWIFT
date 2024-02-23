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
#ifndef SWIFT_DEFAULT_DARK_MATTER_IACT_H
#define SWIFT_DEFAULT_DARK_MATTER_IACT_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "random.h"
#include "dark_matter.h"
#include "kernel_dark_matter.h"
#include "dark_matter_logger.h"
#include "timestep_sync_part.h"
#include "units.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_dark_matter_density(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H) {}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_density(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    const struct dmpart* pj, float a, float H) {}

/**
 * @brief Perform the 'kick' operation on both #gparts
 *
 * @param gpj #gpart
 * @param gpi #gpart
 * @param ti_current Current integer time (for random numbers).
 *
 */
__attribute__((always_inline)) INLINE static void sidm_do_kick(struct dmpart *restrict pi,
                                                               struct dmpart *restrict pj, const integertime_t ti_current,
                                                               /*struct sidm_history* sidm_history,*/
                                                               const struct cosmology* cosmo,
                                                               const struct sidm_props* sidm_props,
                                                               const struct unit_system* us) {}


/**
 * @brief Interaction between two dark matter particles during force loop
 * It calculates the probability of DM particles i & j of scattering within the next time step
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param ti_current Current integer time (for random numbers).
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, float time_base, integertime_t t_current,
    const struct cosmology* cosmo, const int with_cosmology, const struct sidm_props* sidm_props,
    const struct unit_system* us) {}

/**
 * @brief Interaction between two dark matter particles during force loop (non-symmetric).
 * It calculates the probability of DM particles i & j of scattering within the next time step
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not active).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_dark_matter_sidm(
    float r2, const float* dx, float hi, float hj, struct dmpart* pi,
    struct dmpart* pj, float a, float H, float time_base, integertime_t t_current,
    const struct cosmology* cosmo, const int with_cosmology, const struct sidm_props* sidm_props,
    const struct unit_system* us) {}


#endif

