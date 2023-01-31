/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#include "rt_unphysical.h"

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
 * @brief Compute the flux between a left state UL and a right
 * state UR along the direction of the unit vector n_unit
 * through a surface of size Anorm.
 *
 * @param UL left state (energy density, fluxes)
 * @param UR right state (energy density, fluxes)
 * @param n_unit unit vector of the direction of the surface
 * @param Anorm size of the surface through which the flux goes
 * @param fluxes (return) the resulting flux
 */
__attribute__((always_inline)) INLINE static void rt_compute_flux(
    float UL[4], float UR[4], const float n_unit[3], const float Anorm,
    float fluxes[4]) {

  /* Unphysical check not necessary here.
   * It's done in gradients_predict as well. */

  const float FLnorm = sqrtf(UL[1] * UL[1] + UL[2] * UL[2] + UL[3] * UL[3]);
  const float FRnorm = sqrtf(UR[1] * UR[1] + UR[2] * UR[2] + UR[3] * UR[3]);

  /* Get the fluxes in the hyperbolic conservation law sense. */
  float hyperFluxL[4][3];
  rt_get_hyperbolic_flux(UL, FLnorm, hyperFluxL);
  float hyperFluxR[4][3];
  rt_get_hyperbolic_flux(UR, FRnorm, hyperFluxR);

#ifdef SWIFT_RT_DEBUG_CHECKS
  rt_check_unphysical_hyperbolic_flux(hyperFluxL);
  rt_check_unphysical_hyperbolic_flux(hyperFluxR);
#endif

  rt_riemann_solve_for_flux(UL, UR, FLnorm, FRnorm, hyperFluxL, hyperFluxR,
                            n_unit, fluxes);

  /* get the actual flux */
  fluxes[0] *= Anorm;
  fluxes[1] *= Anorm;
  fluxes[2] *= Anorm;
  fluxes[3] *= Anorm;
}

/**
 * @brief reset the mass fluxes of constituent species for a particle.
 *
 * @param p particle to work on.
 **/
__attribute__((always_inline)) INLINE static void rt_part_reset_mass_fluxes(
    struct part* restrict p) {
  p->rt_data.mass_flux.HI = 0.f;
  p->rt_data.mass_flux.HII = 0.f;
  p->rt_data.mass_flux.HeI = 0.f;
  p->rt_data.mass_flux.HeII = 0.f;
  p->rt_data.mass_flux.HeIII = 0.f;
}

#endif /* SWIFT_GEAR_RT_FLUX_H */
