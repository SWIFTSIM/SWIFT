/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H

#include "chemistry_getters.h"
#include "chemistry_riemann_HLL.h"
#include "chemistry_unphysical.h"

/**
 * @brief Reset the diffusion fluxes for the given particle.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void
chemistry_reset_chemistry_fluxes(struct part* restrict p) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; ++i) {
    p->chemistry_data.diffusion_flux[i] = 0.0;
  }
}

/**
 * @brief Get the diffusion fluxes for the given particle.
 *
 * Notice that this function returns the solution to the Riemann problem. Hence
 * the flux is 1D.
 *
 * @param p Particle.
 * @param metal Metal index.
 * @param flux Fluxes for the particle (array of size 1).
 */
__attribute__((always_inline)) INLINE void chemistry_get_fluxes(
    const struct part* restrict p, int metal, double* flux) {
  *flux = p->chemistry_data.diffusion_flux[metal];
}

/**
 * @brief Compute the diffusion flux of given metal group.
 *
 * F_diss = - K * \nabla \otimes q
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param F_diff (return) Array to write diffusion flux component into
 */
__attribute__((always_inline)) INLINE static void
chemistry_compute_physical_diffusion_flux(
    const struct part* restrict p, int metal, double F_diff[3],
    const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo) {

  /* In physical units */
  const double kappa = p->chemistry_data.kappa;

  /* The gradient needs to be converted to physical units:
                                 grad_p = a^{-1} * grad_c.
     The metallicity is already physical (Z_p = Z_c = Z). */

  if (chem_data->diffusion_mode == isotropic_constant ||
      chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Isotropic diffusion: K = kappa * I_3. */
    F_diff[0] = -kappa * p->chemistry_data.gradients.Z[metal][0] * cosmo->a_inv;
    F_diff[1] = -kappa * p->chemistry_data.gradients.Z[metal][1] * cosmo->a_inv;
    F_diff[2] = -kappa * p->chemistry_data.gradients.Z[metal][2] * cosmo->a_inv;
  } else {
    /* Initialise to the flux to 0 */
    F_diff[0] = 0.0;
    F_diff[1] = 0.0;
    F_diff[2] = 0.0;

    /* Compute diffusion matrix K */
    double K[3][3];
    chemistry_get_physical_matrix_K(p, K, chem_data, cosmo);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        F_diff[i] +=
            K[i][j] * p->chemistry_data.gradients.Z[metal][j] * cosmo->a_inv;
      }
    } /* End of matrix multiplication */
  } /* end of if else diffusion_mode */
}

/**
 * @brief Compute the flux for the Riemann problem with the given left and right
 * state, and interface normal, surface area and velocity.
 *
 * @param pi The #part pi.
 * @param pj The #part pj.
 * @param UL Left diffusion state variables (in physical units).
 * @param UR Right diffusion state variables (in physical units).
 * @param WL Left state variables (in physical units).
 * @param WR Right state variables (in physical units).
 * @param n_unit Unit vector of the interface.
 * @param Anorm Surface area of the interface (in physical units).
 * @param metal Metal specie.
 * @param metal_flux Array to store the result in (of size 1).
 */
__attribute__((always_inline)) INLINE static void chemistry_compute_flux(
    const struct part* restrict pi, const struct part* restrict pj,
    const double UL, const double UR, const float WL[5], const float WR[5],
    const float n_unit[3], const float Anorm, int metal, double* metal_flux,
    const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo) {

  /* Note: F_diff_R and F_diff_L are computed with a first order
         reconstruction */
  /* Get the diffusion flux */
  double F_diff_i[3], F_diff_j[3];
  chemistry_compute_physical_diffusion_flux(pi, metal, F_diff_i, chem_data,
                                            cosmo);
  chemistry_compute_physical_diffusion_flux(pj, metal, F_diff_j, chem_data,
                                            cosmo);

#ifdef SWIFT_DEBUG_CHECKS
  chemistry_check_unphysical_diffusion_flux(F_diff_i);
  chemistry_check_unphysical_diffusion_flux(F_diff_j);
#endif

  /* While solving the Riemann problem, we shall get a scalar because of the
     scalar product betwee F_diff_ij^* and A_ij */
  chemistry_riemann_solve_for_flux(pi, pj, UL, UR, WL, WR, F_diff_i, F_diff_j,
                                   Anorm, n_unit, metal, metal_flux, chem_data,
                                   cosmo);

  /* Anorm_physical = a^2 A_comoving. Hence we are missing a a^2 factor. */
  *metal_flux *= Anorm;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_FLUX_H  */