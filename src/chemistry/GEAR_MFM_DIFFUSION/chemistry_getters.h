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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GETTERS_H

#include "const.h"
#include "hydro.h"
#include "part.h"

/**
 * @brief Get a 1-element state vector U containing the metal mass density
 * a specific metal group.
 *
 * @param p Particle.
 * @param group Index of photon group
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_get_diffusion_state_vector(const struct part *restrict p,
                                          int group, double *U) {

  /* The state vector is 1D and contains the metal density. */
  *U = p->rho * p->chemistry_data.metal_mass[group] / hydro_get_mass(p);

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
chemistry_part_compute_diffusion_flux(const struct part *restrict p, int metal,
                                      double F_diff[3]) {

  const double kappa = p->chemistry_data.kappa;

  /* For isotropic diffusion, \grad U = \nabla \otimes q = \grad n_Z.
     Note: K = kappa * I_3 for isotropic diffusion. */
  F_diff[0] = -kappa * p->chemistry_data.gradients[metal].nabla_otimes_q[0];
  F_diff[1] = -kappa * p->chemistry_data.gradients[metal].nabla_otimes_q[1];
  F_diff[2] = -kappa * p->chemistry_data.gradients[metal].nabla_otimes_q[2];
}

/**
 * @brief Get the gradients of metal mass density a given metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param dFx (return) Array to write flux component gradient into
 */
__attribute__((always_inline)) INLINE static void chemistry_part_get_gradients(
    const struct part *restrict p, int metal, double dF[3]) {

  /* For isotropic diffusion, \grad U = \nabla \otimes q = \grad n_Z */
  dF[0] = p->chemistry_data.gradients[metal].nabla_otimes_q[0];
  dF[1] = p->chemistry_data.gradients[metal].nabla_otimes_q[1];
  dF[2] = p->chemistry_data.gradients[metal].nabla_otimes_q[2];
}

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
chemistry_part_geometry_well_behaved(const struct part *restrict p) {

  return p->chemistry_data.geometry.wcorr > const_gizmo_min_wcorr;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GETTERS_H  */
