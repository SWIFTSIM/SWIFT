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

#ifndef SWIFT_GEAR_RT_RIEMANN_HLL_H
#define SWIFT_GEAR_RT_RIEMANN_HLL_H

/**
 * @file src/rt/GEAR/rt_riemann_HLL.h
 * @brief  The Hartmann-Lax-van Leer Riemann solver for the moments of the
 * radiative transfer equation following Rosdahl et al 2013
 * */

/**
 * @brief Solve the Riemann problem for the RT equations and return the
 * flux at the interface.
 *
 * @param UL left state
 * @param UR right state
 * @param flux_half the resulting flux at the interface
 * @param n_unit the unit vector perpendicular to the "intercell" surface.
 */
__attribute__((always_inline)) INLINE static void rt_riemann_solve_for_flux(
    const float UL[4], const float UR[4], float flux_half[4],
    const float n_unit[3]) {

  error("RT HLL Riemann Solver is not implemented yet, sorry :(");
}

#endif /* SWIFT_GEAR_RT_RIEMANN_HLL_H */
