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

#include "chemistry_unphysical.h"
#include "chemistry_getters.h"
#include "chemistry_riemann_HLL.h"
/* #include "chemistry_unphysical.h" */

/**
 * @brief Reset the hydrodynamical fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_part_reset_chemistry_fluxes(
    struct part* restrict p) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.diffusion_flux[i] = 0.0;
  }
}

/**
 * @brief Get the fluxes for the given particle.
 *
 * @param p Particle.
 * @param flux Fluxes for the particle (array of size 5 or more).
 */
__attribute__((always_inline)) INLINE void chemistry_part_get_fluxes(const struct part* restrict p, int metal, double* flux) {
  *flux = p->chemistry_data.diffusion_flux[metal];
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
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
      const struct part* restrict pi, const struct part* restrict pj,
      const double UL, const double UR, const float n_unit[3], int g,
      const float Anorm, const double F_diff_i[3], const double F_diff_j[3], double* metal_flux) {

  /* While solving the Riemann problem, we shall get a scalar because of the
     scalar product betwee F_diff_ij^* and A_ij */
  chemistry_riemann_solve_for_flux(pi, pj, UL, UR, n_unit, g, Anorm, F_diff_i, F_diff_j, metal_flux);

  *metal_flux *= Anorm;
}


#endif /* SWIFT_GEAR_CHEMISTRY_FLUX_H  */
