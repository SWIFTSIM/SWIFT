/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H
#define SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H

#include "../chemistry_riemann_HLL.h"
#include "../chemistry_unphysical.h"

/**
 * @file src/chemistry/GEAR_MF_DIFFUSION/hyperbolic/chemistry_flux.h
 * @brief Header dealing with hyperbolic diffusion fluxes.
 */

/**
 * @brief Reset the metal mass fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_reset_fluxes(struct part* restrict p) {
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.metal_mass_riemann[i] = 0.0;
    p->chemistry_data.flux_riemann[metal][0] = 0.0;
    p->chemistry_data.flux_riemann[metal][1] = 0.0;
    p->chemistry_data.flux_riemann[metal][2] = 0.0;
  }
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the left of the interparticle interface.
 *
 * The sign convention is that a positive total flux is subtracted from the
 * left state, and added to the right state, based on how we chose the unit
 * vector.
 *
 * @param p Particle.
 * @param metal Metal specie.
 * @param fluxes Fluxes accross the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_fluxes_left(struct part* restrict p,
				  const int metal,
				  const double fluxes[4],
				  const float dt) {
  p->chemistry_data.metal_mass_riemann[metal] -= fluxes[0] * dt;
  p->chemistry_data.flux_riemann[metal][0] -= fluxes[1] * dt;
  p->chemistry_data.flux_riemann[metal][1] -= fluxes[2] * dt;
  p->chemistry_data.flux_riemann[metal][2] -= fluxes[3] * dt;
}

/**
 * @brief Update the fluxes for the particle with the given contributions,
 * assuming the particle is to the right of the interparticle interface.
 *
 * The sign convention is that a positive total flux is subtracted from the
 * left state, and added to the right state, based on how we chose the unit
 * vector.
 *
 * @param p Particle.
 * @param metal Metal specie.
 * @param fluxes Fluxes accross the interface.
 * @param dt Time step for the flux exchange.
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_update_fluxes_right(struct part* restrict p,
				   const int metal,
				   const double fluxes[4],
				   const float dt) {
  p->chemistry_data.metal_mass_riemann[metal] += fluxes[0] * dt;
  p->chemistry_data.flux_riemann[metal][0] += fluxes[1] * dt;
  p->chemistry_data.flux_riemann[metal][1] += fluxes[2] * dt;
  p->chemistry_data.flux_riemann[metal][2] += fluxes[3] * dt;
}

// TODO: Rewrite this function (and the riemann solver part to use the
// hyperbolic version)
/**
 * @brief Compute the metal mass flux for the Riemann problem with the given
 * left and right state, and interface normal, surface area and velocity.
 *
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param pi The #part pi.
 * @param pj The #part pj.
 * @param metal Metal specie.
 * @param UL Left diffusion state variables (in physical units).
 * @param UR Right diffusion state variables (in physical units).
 * @param WL Left state variables (in physical units).
 * @param WR Right state variables (in physical units).
 * @param n_unit Unit vector of the interface.
 * @param Anorm Surface area of the interface (in physical units).
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 * @param fluxes (return) The resulting flux at the interface (of size 1).
 */
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
    const float dx[3], const struct part* restrict pi,
    const struct part* restrict pj, const int metal,
    const double UL[4], const double UR[4],
    const float WL[5], const float WR[5], const float n_unit[3],
    const float Anorm, const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo, double fluxes[4]) {

  /* Use the predicted fluxes. They improve metal mass conservation and reduce
     artificial diffusion. */
  double F_diff_i[3] = {
      pi->chemistry_data.hyperbolic_flux[metal].F_diff_pred[0],
      pi->chemistry_data.hyperbolic_flux[metal].F_diff_pred[1],
      pi->chemistry_data.hyperbolic_flux[metal].F_diff_pred[2]};
  double F_diff_j[3] = {
      pj->chemistry_data.hyperbolic_flux[metal].F_diff_pred[0],
      pj->chemistry_data.hyperbolic_flux[metal].F_diff_pred[1],
      pj->chemistry_data.hyperbolic_flux[metal].F_diff_pred[2]};

  /* Check that the fluxes are meaningful */
  chemistry_check_unphysical_diffusion_flux(F_diff_i);
  chemistry_check_unphysical_diffusion_flux(F_diff_j);

  /* While solving the Riemann problem, we shall get a scalar because of the
     scalar product betwee F_diff_ij^* and A_ij */
  chemistry_riemann_solve_for_flux(dx, pi, pj, UL, UR, WL, WR, F_diff_i,
				   F_diff_j, Anorm, n_unit, metal, chem_data,
				   cosmo, metal_flux);

  /* Anorm is already in physical units here. */
  *metal_flux *= Anorm;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION_FLUX_H */
