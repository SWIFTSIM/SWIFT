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
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
      const struct part* restrict pi, const struct part* restrict pj,
      const float UL, const float UR, const float n_unit[3],
      const float Anorm, const float F_diff_i[3], const float F_diff_j[3], float metal_flux) {

  /* TODO */
  /* While solving the Riemann problem, we shall get a scalar because of the
     scalar product betwee F_diff_ij^* and A_ij */
  /* Note: We should give the particle to the riemann solver, because we need
     to compute lambda_+/- from the particles' properties */  
  
  chemistry_riemann_solve_for_flux(pi, pj, UL, UR, n_unit, Anorm, F_diff_i, F_diff_j, metal_flux);

  metal_flux *= Anorm;
}


#endif /* SWIFT_GEAR_CHEMISTRY_FLUX_H  */
