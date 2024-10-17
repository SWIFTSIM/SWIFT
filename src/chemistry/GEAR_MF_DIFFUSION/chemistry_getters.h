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
#ifndef SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GETTERS_H
#define SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GETTERS_H

#include "chemistry_struct.h"
#include "const.h"
#include "cosmology.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "part.h"

#ifdef HAVE_LIBGSL
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#endif

/**
 * @brief Get metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double chemistry_get_comoving_metal_density(
    const struct part *restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / p->geometry.volume;
}

/**
 * @brief Get metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double chemistry_get_physical_metal_density(
  const struct part *restrict p, int metal, const struct cosmology *cosmo) {
  return cosmo->a3_inv*chemistry_get_comoving_metal_density(p, metal);
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
 * @brief Get a 1-element state vector U containing the metal mass density (in
 * comoving units) of a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param U Pointer to the array in which the result needs to be stored
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_comoving_diffusion_state_vector(const struct part *restrict p, int metal) {
  /* The state vector is 1D and contains the metal density. */
  return chemistry_get_comoving_metal_density(p, metal);
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
__attribute__((always_inline)) INLINE static float chemistry_get_comoving_density(
    const struct part *restrict p) {
  float rho = hydro_get_comoving_density(p);

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
 * @brief Compute the physical diffusion coefficient of the particle.
 *
 * This must be called in chemistry_prepare_force() to have the values of the
 * density and the matrix S (which depends on grad_v_tilde).
 *
 * Note: The diffusion coefficient depends on the particle's density. If the
 * density is 0, then the coefficient is 0 as well and the timestep for
 * chemistry will be 0.
 *
 * @param p Particle.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_compute_diffusion_coefficient(
  struct part *restrict p, const struct chemistry_global_data *chem_data,
  const struct cosmology *cosmo) {

  float rho = p->chemistry_data.filtered.rho;

  /* In case the filtered density is 0, e.g. during the fake-timestep,
     approximate the density */
  if (rho == 0.0) {
    rho = chemistry_get_comoving_density(p);
  }

  /* Convert density to physical units. */
  rho *= cosmo->a3_inv;

  /* Convert smoothing length to physical units */
  const double h2_p = cosmo->a*cosmo->a * p->h * p->h;

  if (chem_data->diffusion_mode == isotropic_constant) {
    return chem_data->diffusion_coefficient;
  } else if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Get the physical shear tensor */
    double S[3][3];
    chemistry_get_physical_shear_tensor(p, S, cosmo);

    /* In the smagorinsky model, we remove the trace from S */
    const double trace = S[0][0] + S[1][1] + S[2][2];

    S[0][0] -= trace;
    S[1][1] -= trace;
    S[2][2] -= trace;

    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho * chemistry_get_matrix_norm(S);
  } else {
    /* Note that this is multiplied by the matrix S to get the full matrix K */
    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho;
  }
}

/**
 * @brief Get the physical shear tensor.
 *
 * @param p Particle.
 * @param S (return) Pointer to a 3x3 matrix.
 */
__attribute__((always_inline)) INLINE static void chemistry_get_physical_shear_tensor(
  const struct part *restrict p, double S[3][3], const struct cosmology* cosmo) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      S[i][j] = 0.5 * (p->chemistry_data.filtered.grad_v_tilde[i][j] +
                       p->chemistry_data.filtered.grad_v_tilde[j][i]);

      /* Convert to physical units.
	 In Swift, v_peculiar = v_c / a ; grad_p = a^{-1} grad_c. */
      S[i][j] *= cosmo->a2_inv;
    }
  }
}

/**
 * @brief Regularize the shear tensor as described in Balarac et al. (2013).
 *
 * This needs to be checked.
 *
 * @param S (return) Pointer to a 3x3 matrix shear tensor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_regularize_shear_tensor(double S[3][3]) {
#ifdef HAVE_LIBGSL
  /* Create the necessary GSL objects to hold the eigenvalues and eigenvectors
   */
  gsl_matrix *S_matrix = gsl_matrix_alloc(3, 3);  // Matrix for input/output S
  gsl_vector *eigenvalues =
      gsl_vector_alloc(3);  // Vector to hold the eigenvalues
  gsl_matrix *eigenvectors =
      gsl_matrix_alloc(3, 3);  // Matrix to hold the eigenvectors
  gsl_eigen_symmv_workspace *workspace =
      gsl_eigen_symmv_alloc(3);  // Workspace for eigen computation

  /* Fill S_matrix with the values from S */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      gsl_matrix_set(S_matrix, i, j, S[i][j]);
    }
  }

  /* Compute the eigenvalues and eigenvectors. S is symmetric by construction.
   */
  gsl_eigen_symmv(S_matrix, eigenvalues, eigenvectors, workspace);

  /* Zero-initialize the matrix S to store S_minus */
  double S_minus[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  /* Compute S_ij^minus using the sum of min(0, lambda^(k)) * e_i^(k) * e_j^(k)
   */
  for (int k = 0; k < 3; k++) {
    double lambda_k =
        gsl_vector_get(eigenvalues, k);        // Get the k-th eigenvalue
    double lambda_k_minus = min(0, lambda_k);  // Take min(0, lambda^(k))

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        double e_ik = gsl_matrix_get(eigenvectors, i, k);
        double e_jk = gsl_matrix_get(eigenvectors, j, k);
        S_minus[i][j] += lambda_k_minus * e_ik * e_jk;
      }
    }
  }

  // Copy S_minus back into S (overwriting it)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      S[i][j] = S_minus[i][j];
    }
  }

  // Free the allocated memory
  gsl_matrix_free(S_matrix);
  gsl_vector_free(eigenvalues);
  gsl_matrix_free(eigenvectors);
  gsl_eigen_symmv_free(workspace);
#else
  error(
      "Code not compiled with GSL. Can't compute eigenvalues of the filtered "
      "shear tensor.");
#endif
}

/**
 * @brief Get the physical diffusion matrix K.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void chemistry_get_physical_matrix_K(
    const struct part *restrict p, double K[3][3],
    const struct chemistry_global_data *chem_data,
    const struct cosmology* cosmo) {
  if (chem_data->diffusion_mode == isotropic_constant ||
      chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* K = kappa * I */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        K[i][j] = 0.0;
      }
    }
    /* kappa is already in physical units */
    K[0][0] = p->chemistry_data.kappa;
    K[1][1] = p->chemistry_data.kappa;
    K[2][2] = p->chemistry_data.kappa;

  } else {
    /* Get the full shear tensor */
    chemistry_get_physical_shear_tensor(p, K, cosmo);

    /* This takes way too much time, probably because we allocate and
       deallocate the workspace too often. Comment it for now */
    /* Now regularize the shear tensor by considering only the negative
       eigenvalues (Balarac et al. (2013)). This is now called the S_minus
       matrix. */
    /* chemistry_regularize_shear_tensor(K); */

    /* K = kappa * S_minus */
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
	/* kappa is already in physical units */
        K[i][j] *= p->chemistry_data.kappa;
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

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GETTERS_H  */
