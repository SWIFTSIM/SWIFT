/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#ifndef SWIFT_GEAR_RT_FLUX_H
#define SWIFT_GEAR_RT_FLUX_H

/**
 * @file src/rt/GEAR/rt_flux.h
 * @brief Functions related to compute the interparticle flux term of the
 * conservation law. This is its own file so we can switch between Riemann
 * solvers more easily.
 */

/**
 * @brief Compute the flux between a left state Qleft and a right
 * state Qright along the direction of the unit vector n_unit
 * through a surface of size Anorm.
 *
 * @param Qleft left state
 * @param Qright right state
 * @param n_unit unit vector of the direction of the surface
 * @param Anorm size of the surface through which the flux goes
 * @param totflux the resulting flux
 */
__attribute__((always_inline)) INLINE static void rt_compute_flux(
    const float Qleft[4], const float Qright[4], const float n_unit[3],
    const float Anorm, float totflux[4]) {}
#endif /* SWIFT_GEAR_RT_FLUX_H */
