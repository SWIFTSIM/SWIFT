/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_GEAR_CHEMISTRY_FLUX_H
#define SWIFT_GEAR_CHEMISTRY_FLUX_H

#include "chemistry_riemann_HLL.h"
#include "chemistry_unphysical.h"

/**
 * @brief Reset the hydrodynamical fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_part_reset_hydro_fluxes(
    struct part* restrict p) {

  p->flux.momentum[0] = 0.0f;
  p->flux.momentum[1] = 0.0f;
  p->flux.momentum[2] = 0.0f;
  p->flux.energy = 0.0f;
}

/**
 * @brief Get the fluxes for the given particle.
 *
 * @param p Particle.
 * @param flux Fluxes for the particle (array of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void hydro_part_get_fluxes(
    const struct part* restrict p, float* flux) {

  flux[1] = p->flux.momentum[0];
  flux[2] = p->flux.momentum[1];
  flux[3] = p->flux.momentum[2];
  flux[4] = p->flux.energy;
}

/**
 * @brief Compute the flux for the Riemann problem with the given left and right
 * state, and interface normal, surface area and velocity.
 *
 * @param WL Left state variables.
 * @param WR Right state variables.
 * @param n_unit Unit vector of the interface.
 * @param vLR Velocity of the interface.
 * @param Anorm Surface area of the interface.
 * @param fluxes Array to store the result in (of size 5 or more).
 */
__attribute__((always_inline)) INLINE static void chemitry_compute_flux(
    const struct part* restrict p
    const float UL, const float UR, const float* n_unit,
    const float Anorm, float F_diff[3], float metal_flux) {

  /* TODO */
  /* While solving the Riemann problem, we shall get a scalar because of the
     scalar product betwee F_diff_ij^* and A_ij */
  /* Note: We should give the particle to the riemann solver, because we need
     to compute lambda_+/- from the particles' properties */  
  
  chemistry_riemann_solve_flux(p, UL, UR, n_unit, metal_flux);

  metal_flux *= Anorm;
}


#endif /* SWIFT_GEAR_CHEMISTRY_FLUX_H  */
