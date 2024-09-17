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
 * @brief Get a  metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double
chemistry_part_get_metal_density(const struct part *restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] /
         p->chemistry_data.geometry.volume;
}

/**
 * @brief Get a 1-element state vector U containing the metal mass density
 * a specific metal group.
 *
 * @TODO: Rewrite this to remove the pointer... We can use return instead.
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static void
chemistry_part_get_diffusion_state_vector(const struct part *restrict p,
                                          int metal, double *U) {

  /* The state vector is 1D and contains the metal density. */
  *U = chemistry_part_get_metal_density(p, metal);
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
 * @brief Get the particle volume.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float chemistry_part_get_volume(
    const struct part *restrict p) {
  return p->chemistry_data.geometry.volume;
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

/**
 * @brief Get matrix K Frobenius norm.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static double
chemistry_compute_matrix_K_norm(const struct part *restrict p) {
  /* For isotropic cases: K = \kappa * I_3 */
  return sqrtf(3.0) * fabsf(p->chemistry_data.kappa);
}

/**
 * @brief Compute the particle parabolic timestep proportional to h^2.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_parabolic_timestep(const struct part *restrict p) {

  const struct chemistry_part_data chd = p->chemistry_data;

  /* Note: The State vector is U = (metal_density_1, metal_density_2, etc),
     which is also = q. Hence, the last term in eq (15) is unity. */

  float norm_q = 0.0;
  float norm_nabla_otimes_q = 0.0;
  float expression = 0.0;
  const float norm_matrix_K = chemistry_compute_matrix_K_norm(p);
  const float delta_x = kernel_gamma * p->h;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_q += p->chemistry_data.metal_mass[i] * p->chemistry_data.metal_mass[i];

    for (int j = 0; j < 3; j++) {
      /* Compute the Froebnius norm of \nabla \otimes q */
      norm_nabla_otimes_q += chd.gradients[i].nabla_otimes_q[j] *
                             chd.gradients[i].nabla_otimes_q[j];
    }
  }
  /* Take the sqrt and divide by the volume to get a density */
  norm_q = sqrtf(norm_q) / chd.geometry.volume;

  /* If the norm of q (metal density = 0), then use the following formula */
  if (norm_q == 0) {
    warning("norm q = 0");
    return delta_x*delta_x / norm_matrix_K;
  }

  /* Take the sqrt to get the norm */
  norm_nabla_otimes_q = sqrtf(norm_nabla_otimes_q);

  /* Finish the computations */
  expression = norm_nabla_otimes_q / norm_q + 1.0 / delta_x;

  return 1.0 / (norm_matrix_K * expression * expression);
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_supertimestep(const struct part *restrict p,
                                const struct chemistry_global_data *cd) {

  const struct chemistry_part_data chd = p->chemistry_data;

  /* Note: The State vector is U = (metal_density_1, metal_density_2, etc),
     which is also = q. Hence, the last term in eq (15) is unity. */

  float norm_U = 0.0;
  float norm_nabla_otimes_q = 0.0;
  const float norm_matrix_K = chemistry_compute_parabolic_timestep(p);
  const float delta_x = kernel_gamma * p->h;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_U += p->chemistry_data.metal_mass[i] * p->chemistry_data.metal_mass[i];

    for (int j = 0; j < 3; j++) {
      /* Compute the Froebnius norm of \nabla \otimes q */
      norm_nabla_otimes_q += chd.gradients[i].nabla_otimes_q[j] *
                             chd.gradients[i].nabla_otimes_q[j];
    }
  }

  /* Take the sqrt and divide by the volume to get a density */
  norm_U = sqrtf(norm_U) / chd.geometry.volume;

  /* Take the sqrt to get the norm */
  norm_nabla_otimes_q = sqrtf(norm_nabla_otimes_q);

  return cd->C_CFL_chemistry * delta_x * norm_U /
         (norm_matrix_K * norm_nabla_otimes_q);
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_subtimestep(const struct part *restrict p,
                              const struct chemistry_global_data *cd) {
  const struct chemistry_part_data chd = p->chemistry_data;
  const float cos_argument = M_PI *
                             (2.0 * chd.timesteps.current_substep - 1.0) /
                             (2.0 * cd->N_substeps);
  const float expression = (1 + cd->nu) - (1 - cd->nu) * cos(cos_argument);
  return chd.timesteps.explicit_timestep / expression;
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_minimal_timestep_from_all_modules(const struct part *restrict p,
                              const struct chemistry_global_data *cd) {
  /* Don't do this. Create a fct that does get_part_timestep() first part
     (without integer time conversion) without chemistry. Then call
     chemistry_timestep(), do supertimestepping here. Finally, take the
     min. See how we can get around to not go below the min timestep and above
     the mac timestep.
  */

  return 0.0;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GETTERS_H  */
