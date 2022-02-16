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

#if defined(RT_RIEMANN_SOLVER_GLF)
#include "rt_riemann_GLF.h"
#elif defined(RT_RIEMANN_SOLVER_HLL)
#include "rt_riemann_HLL.h"
#else
#error "No valid choice of RT Riemann solver has been selected"
#endif

/**
 * @file src/rt/GEAR/rt_flux.h
 * @brief Functions related to compute the interparticle flux term of the
 * conservation law. This is its own file so we can switch between Riemann
 * solvers more easily.
 */

/**
 * @brief Reset the fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void rt_part_reset_fluxes(
    struct part* restrict p) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.flux[g].energy = 0.f;
    p->rt_data.flux[g].flux[0] = 0.f;
    p->rt_data.flux[g].flux[1] = 0.f;
    p->rt_data.flux[g].flux[2] = 0.f;
  }
}

/**
 * @brief Compute the flux between a left state Qleft and a right
 * state Qright along the direction of the unit vector n_unit
 * through a surface of size Anorm.
 *
 * @param UL left state
 * @param UR right state
 * @param n_unit unit vector of the direction of the surface
 * @param Anorm size of the surface through which the flux goes
 * @param fluxes the resulting flux
 */
__attribute__((always_inline)) INLINE static void rt_compute_flux(
    float UL[4], float UR[4], const float n_unit[3], const float Anorm,
    float fluxes[4]) {

  rt_check_unphysical_density(&UL[0], &UL[1], 2);
  rt_check_unphysical_density(&UR[0], &UR[1], 2);

  rt_riemann_solve_for_flux(UL, UR, fluxes, n_unit);

  /* get the actual flux */
  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
}

#endif /* SWIFT_GEAR_RT_FLUX_H */
