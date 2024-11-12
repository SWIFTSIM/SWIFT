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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H

#include "chemistry_getters.h"
#include "chemistry_struct.h"

/**
 * @brief Compute the particle parabolic timestep proportional to h^2.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_parabolic_timestep(
    const struct part *restrict p,
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute diffusion matrix K */
  double K[3][3];
  chemistry_get_physical_matrix_K(p, K, chem_data, cosmo);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...), and q = (Z_1, Z_2,
     ...). Hence, the term norm(U)/norm(q) in eq (15) is abs(rho). */
  const float norm_U_over_norm_q =
      cosmo->a3_inv * fabs(chemistry_get_comoving_density(p));

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h * cosmo->a;
  float norm_q = 0.0;
  float norm_nabla_q = 0.0;
  float expression = 0.0;

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0) {
    return FLT_MAX;
  }

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_q += chemistry_get_metal_mass_fraction(p, i) *
              chemistry_get_metal_mass_fraction(p, i);

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q = Grad Z */
      norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt. Add the missing a^{-1} to have physical gradients */
  norm_q = sqrtf(norm_q);  // Physical
  norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv;

  /* If the norm of q (metal density = 0), then use the following
     expression. Notice that if norm q = 0, the true timestep muste be 0... */
  if (norm_q == 0) {
    return delta_x * delta_x / norm_matrix_K * norm_U_over_norm_q;
  }

  /* Compute the expression in the square bracket in eq (15). Notice that I
     rewrote it to avoid division by 0 when norm_nabla_1 = 0. */
  expression = norm_q * delta_x / (norm_nabla_q * delta_x + norm_q);

  /* This trick comes from Gizmo. This is a mix between eq (15) and eq (D1) */
  expression = max(expression, delta_x);

  return expression * expression / norm_matrix_K * norm_U_over_norm_q;
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * This is equation (10) in Alexiades, Amiez and Gremaud (1996).
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float chemistry_get_supertimestep(
    const struct part *restrict p, const struct chemistry_global_data *cd,
    float timestep_explicit) {
  const float N = cd->N_substeps;
  const float nu = cd->nu;
  const float nu_plus_term = pow(1 + sqrtf(nu), 2.0 * N);
  const float nu_minus_term = pow(1 - sqrtf(nu), 2.0 * N);
  const float left_term =
      (nu_plus_term - nu_minus_term) / (nu_plus_term + nu_minus_term);
  return timestep_explicit * N / (2.0 * sqrt(nu)) * left_term;
}

/**
 * @brief Compute the particle supertimestep with using a CFL-like condition,
 * proportioanl to h.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_CFL_supertimestep(const struct part *restrict p,
                                    const struct chemistry_global_data *cd,
                                    const struct cosmology *cosmo) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute diffusion matrix K */
  double K[3][3];
  chemistry_get_physical_matrix_K(p, K, cd, cosmo);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h * cosmo->a;

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...). */
  float norm_U = 0.0;
  float norm_nabla_q = 0.0;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_U += chemistry_get_physical_metal_density(p, i, cosmo) *
              chemistry_get_physical_metal_density(p, i, cosmo);

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q */
      norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt and convert to physical units */
  norm_U = sqrtf(norm_U);
  norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv;

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0 || norm_nabla_q == 0.0) {
    return FLT_MAX;
  }

  return cd->C_CFL_chemistry * delta_x * norm_U /
         (norm_matrix_K * norm_nabla_q);
}

/**
 * @brief Compute the particle supertimestep proportional to h.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_subtimestep(const struct part *restrict p,
                              const struct chemistry_global_data *cd,
                              int current_substep_number) {
  const struct chemistry_part_data *chd = &p->chemistry_data;
  const float cos_argument =
      M_PI * (2.0 * current_substep_number - 1.0) / (2.0 * cd->N_substeps);
  const float expression = (1 + cd->nu) - (1 - cd->nu) * cos(cos_argument);
  return chd->timesteps.explicit_timestep / expression;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_TIMESTEPS_H  */
