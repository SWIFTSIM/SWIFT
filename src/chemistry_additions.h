/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Yolan Uyttenhove (yolan.uyttenhove@ugent.be)
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
#ifndef SWIFT_CHEMISTRY_ADDITIONS_H
#define SWIFT_CHEMISTRY_ADDITIONS_H

/**
 * @file src/chemistry_additions.h
 * @brief Branches between the different additional functions required outside
 * of the chemistry files (e.g. in hydro loops);
 * Specifically for functions used for advection of tracked elements for
 * hydro schemes with mass fluxes.
 **/

/* Config parameters. */
#include <config.h>

#ifdef HYDRO_DOES_MASS_FLUX
/* Import the right chemistry definition */
#if defined(CHEMISTRY_AGORA)
#include "./chemistry/AGORA/chemistry_additions.h"
#elif defined(CHEMISTRY_EAGLE)
#include "./chemistry/EAGLE/chemistry_additions.h"
#elif defined(CHEMISTRY_GEAR)
#include "./chemistry/GEAR/chemistry_additions.h"
#elif defined(CHEMISTRY_NONE)
#include "./chemistry/none/chemistry_additions.h"
#elif defined(CHEMISTRY_NONE)
#include "./chemistry/QLA/chemistry_additions.h"
#else
#error "Metal advection unimpmlemented for selected chemistry scheme!"
#endif
#else
/**
 * @brief Extra operations done during the kick. This needs to be
 * done before the particle mass is updated in the hydro_kick_extra.
 *
 * @param p Particle to act upon.
 * @param dt_therm Thermal energy time-step @f$\frac{dt}{a^2}@f$.
 * @param dt_grav Gravity time-step @f$\frac{dt}{a}@f$.
 * @param dt_hydro Hydro acceleration time-step
 * @f$\frac{dt}{a^{3(\gamma{}-1)}}@f$.
 * @param dt_kick_corr Gravity correction time-step @f$adt@f$.
 * @param cosmo Cosmology.
 * @param hydro_props Additional hydro properties.
 */
__attribute__((always_inline)) INLINE static void chemistry_kick_extra(
    struct part* p, float dt_therm, float dt_grav, float dt_hydro,
    float dt_kick_corr, const struct cosmology* cosmo,
    const struct hydro_props* hydro_props) {}

/**
 * @brief update metal mass fluxes between two interacting particles during
 * hydro_iact_(non)sym(...) calls.
 *
 * @param pi first interacting particle
 * @param pj second interacting particle
 * @param mass_flux the mass flux between these two particles.
 * @param flux_dt the time-step over which the fluxes are exchanged
 * @param mode 0: non-symmetric interaction, update i only. 1: symmetric
 * interaction.
 **/
__attribute__((always_inline)) INLINE static void runner_iact_chemistry_fluxes(
    struct part* restrict pi, struct part* restrict pj, float mass_flux,
    float flux_dt, int mode) {}
#endif

#endif  // SWIFT_CHEMISTRY_ADDITIONS_H
