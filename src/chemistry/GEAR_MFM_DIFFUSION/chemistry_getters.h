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

#include "chemistry_struct.h"
#include "const.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "part.h"

/**
 * @brief Check if the gradient matrix for this particle is well behaved.
 *
 * @param p Particle.
 * @return 1 if the gradient matrix is well behaved, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
chemistry_geometry_well_behaved(const struct part *restrict p) {
  return p->chemistry_data.geometry.wcorr > const_gizmo_min_wcorr;
}

/**
 * @brief Get the particle volume.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float chemistry_get_volume(
    const struct part *restrict p) {
  return p->chemistry_data.geometry.volume;
}

/**
 * @brief Get a  metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double chemistry_get_metal_density(
    const struct part *restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / chemistry_get_volume(p);
}

/**
 * @brief Get a  metal mass fraction from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_metal_mass_fraction(const struct part *restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / hydro_get_mass(p);
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
chemistry_get_diffusion_state_vector(const struct part *restrict p, int metal,
                                     double *U) {
  /* The state vector is 1D and contains the metal density. */
  *U = chemistry_get_metal_density(p, metal);
}

/**
 * @brief Get particle density.
 *
 * This function must be used for sensitive operations like computing
 * timesteps. At the beggining of a simulation, it can happen that the
 * particle's density is 0 (e.g. not read from ICs) and not yet updated. Since
 * timesteps computations and the diffusion coefficient require the density, we
 * need to estimate it. Otherwise we have null timesteps. This is particularly
 * true with MFM SPH.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float chemistry_get_density(
    const struct part *restrict p) {
  float rho = p->rho;

  if (rho == 0.0) {
    const float r_cubed =
        kernel_gamma * kernel_gamma * kernel_gamma * p->h * p->h * p->h;
    const float volume = 4.0 / 3.0 * M_PI * r_cubed;
    rho = hydro_get_mass(p) / volume;

    if (rho == 0.0) {
      rho = FLT_MIN;
      error("Density cannot be null!");
    }
  }
  return rho;
}


/**
 * @brief Get the shear tensor.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_get_shear_tensor(
    const struct part *restrict p, double S[3][3]) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        S[i][j] = 0.5*(p->chemistry_data.filtered.grad_v_tilde[i][j] +
		       p->chemistry_data.filtered.grad_v_tilde[j][i]);
      }
    }
}


/**
 * @brief Get the diffusion matrix K.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_get_matrix_K(
    const struct part *restrict p, double kappa, double K[3][3],
    const struct chemistry_global_data *chem_data) {
  if (chem_data->diffusion_mode == isotropic_constant ||
      chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* K = kappa * I */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K[i][j] = 0.0;
      }
    }
    K[0][0] = kappa;
    K[1][1] = kappa;
    K[2][2] = kappa;
  } else {
    /* K = kappa * S */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K[i][j] = kappa * (p->chemistry_data.filtered.grad_v_tilde[i][j] +
                           p->chemistry_data.filtered.grad_v_tilde[j][i]);
      }
    }
  }
}

/**
 * @brief Get matrix K Frobenius norm.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static double chemistry_get_matrix_norm(
    const double K[3][3]) {
  float norm = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      norm += K[i][j] * K[i][j];
    }
  }
  return sqrtf(norm);
}

/**
 * @brief Compute the diffusion coefficient of the particle.
 *
 * This must be called in chemistry_prepare_force() to have the values of the
 * density and the matrix S (which depends on grad_v_tilde).
 *
 * Note: The diffusion coefficient depends on the particle's density. If the
 * density is 0, then the coefficient is 0 as well and the timestep for
 * chemistry will be 0.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static double
chemistry_compute_diffusion_coefficient(
    struct part *restrict p, const struct chemistry_global_data *chem_data) {

  float rho = p->chemistry_data.filtered.rho;

  if (rho == 0.0) {
    rho = chemistry_get_density(p);
  }

  if (chem_data->diffusion_mode == isotropic_constant) {
    return chem_data->diffusion_coefficient;
  } else if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Compute K with kappa = 1 ==> get matrix S */
    double S[3][3];
    chemistry_get_matrix_K(p, 1.0, S, chem_data);

    /* In the smagorinsky model, we remove the trace from S */
    const double trace = S[0][0] + S[1][1] + S[2][2];

    S[0][0] -= trace;
    S[1][1] -= trace;
    S[2][2] -= trace;

    return chem_data->diffusion_coefficient * kernel_gamma2 * p->h * p->h *
           rho * chemistry_get_matrix_norm(S);
  } else {
    /* Note that this is multiplied by the matrix S to get the full matrix K */
    return chem_data->diffusion_coefficient * kernel_gamma2 * p->h * p->h * rho;
  }
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
chemistry_compute_diffusion_flux(const struct part *restrict p, int metal,
                                 const struct chemistry_global_data *chem_data,
                                 double F_diff[3]) {

  const double kappa = p->chemistry_data.kappa;

  if (chem_data->diffusion_mode == isotropic_constant ||
      chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Isotropic diffusion: K = kappa * I_3 for isotropic diffusion. */
    F_diff[0] = -kappa * p->chemistry_data.gradients.Z[metal][0];
    F_diff[1] = -kappa * p->chemistry_data.gradients.Z[metal][1];
    F_diff[2] = -kappa * p->chemistry_data.gradients.Z[metal][2];
  } else {
    F_diff[0] = 0.0;
    F_diff[1] = 0.0;
    F_diff[2] = 0.0;

    /* Compute diffusion matrix K_star = 0.5*(KR + KL) */
    double K[3][3];
    chemistry_get_matrix_K(p, p->chemistry_data.kappa, K, chem_data);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        F_diff[i] += K[i][j] * p->chemistry_data.gradients.Z[metal][j];
      }
    } /* End of matrix multiplication */
  } /* end of if else diffusion_mode */
}

/**
 * @brief Get the gradients of metal mass density a given metal group.
 *
 * Get grad U = grad rho_Z.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param dvx x velocity gradient (of size 3 or more).
 * @param dvy y velocity gradient (of size 3 or more).
 * @param dvz z velocity gradient (of size 3 or more).
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_diffusion_gradients(const struct part *restrict p, int metal,
                                  const float grad_rho[3], double dF[3]) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* We have U = rho_Z and q = Z.
     But we computed Grad Z and not Grad (rho*Z).
     However, Grad (rho*Z) = Z*Grad_rho + rho*Grad_Z
     We can estimate grad_rho = (rho_max_ij - rho_min_ij) * dx[3] / (r*r). */

  const double Z = chemistry_get_metal_mass_fraction(p, metal);

  /* For isotropic diffusion, \grad U = \nabla \otimes q = \grad n_Z */
  dF[0] = chd->gradients.Z[metal][0] * p->rho + grad_rho[0] * Z;
  dF[1] = chd->gradients.Z[metal][1] * p->rho + grad_rho[1] * Z;
  dF[2] = chd->gradients.Z[metal][2] * p->rho + grad_rho[2] * Z;
}

/**
 * @brief Get the velocity gradients.
 *
 * @param p Particle.
 * @param dvx x velocity gradient (of size 3 or more).
 * @param dvy y velocity gradient (of size 3 or more).
 * @param dvz z velocity gradient (of size 3 or more).
 */
__attribute__((always_inline)) INLINE static void chemistry_get_hydro_gradients(
    const struct part *restrict p, float dvx[3], float dvy[3], float dvz[3]) {

  dvx[0] = p->chemistry_data.gradients.v[0][0];
  dvx[1] = p->chemistry_data.gradients.v[0][1];
  dvx[2] = p->chemistry_data.gradients.v[0][2];
  dvy[0] = p->chemistry_data.gradients.v[1][0];
  dvy[1] = p->chemistry_data.gradients.v[1][1];
  dvy[2] = p->chemistry_data.gradients.v[1][2];
  dvz[0] = p->chemistry_data.gradients.v[2][0];
  dvz[1] = p->chemistry_data.gradients.v[2][1];
  dvz[2] = p->chemistry_data.gradients.v[2][2];
}

/**
 * @brief Compute the particle parabolic timestep proportional to h^2.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static float
chemistry_compute_parabolic_timestep(
    const struct part *restrict p,
    const struct chemistry_global_data *chem_data) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute diffusion matrix K */
  double K[3][3];
  chemistry_get_matrix_K(p, chd->kappa, K, chem_data);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...), and q = (Z_1, Z_2,
     ...). Hence, the term norm(U)/norm(q) in eq (15) is abs(rho). */
  const float norm_U_over_norm_q = fabs(chemistry_get_density(p));

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h;
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

  /* Take the sqrt */
  norm_q = sqrtf(norm_q);
  norm_nabla_q = sqrtf(norm_nabla_q);

  /* If the norm of q (metal density = 0), then use the following
     expression. Notice that if norm q = 0, the true timestep muste be 0... */
  if (norm_q == 0) {
    return delta_x * delta_x / norm_matrix_K * norm_U_over_norm_q;
  }

  /* Finish the computations */
  expression = norm_q * delta_x / (norm_nabla_q * delta_x + norm_q);

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
                                    const struct chemistry_global_data *cd) {

  const struct chemistry_part_data *chd = &p->chemistry_data;

  /* Compute diffusion matrix K */
  double K[3][3];
  chemistry_get_matrix_K(p, chd->kappa, K, cd);
  const float norm_matrix_K = chemistry_get_matrix_norm(K);

  /* Some helpful variables */
  const float delta_x = kernel_gamma * p->h;

  /* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...). */
  float norm_U = 0.0;
  float norm_nabla_otimes_q = 0.0;

  /* Compute the norms */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    norm_U += p->chemistry_data.metal_mass[i] * p->chemistry_data.metal_mass[i];

    for (int j = 0; j < 3; j++) {
      /* Compute the Frobenius norm of \nabla \otimes q */
      norm_nabla_otimes_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j];
    }
  }

  /* Take the sqrt and divide by the volume to get a density */
  norm_U = sqrtf(norm_U) / chemistry_get_volume(p);

  /* Take the sqrt to get the norm */
  norm_nabla_otimes_q = sqrtf(norm_nabla_otimes_q);

  /* Prevent pathological cases */
  if (norm_matrix_K == 0.0 || norm_nabla_otimes_q == 0.0) {
    return FLT_MAX;
  }

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
                              const struct chemistry_global_data *cd,
                              int current_substep_number) {
  const struct chemistry_part_data *chd = &p->chemistry_data;
  const float cos_argument =
      M_PI * (2.0 * current_substep_number - 1.0) / (2.0 * cd->N_substeps);
  const float expression = (1 + cd->nu) - (1 - cd->nu) * cos(cos_argument);
  return chd->timesteps.explicit_timestep / expression;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_GETTERS_H  */
