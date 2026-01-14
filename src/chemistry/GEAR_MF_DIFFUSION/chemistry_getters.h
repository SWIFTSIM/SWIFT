/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@ealumni.pfl.ch)
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
#include "chemistry_utils.h"
#include "const.h"
#include "cosmology.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "part.h"

/**
 * @brief Get metal mass fraction from a specific metal specie.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_metal_mass_fraction(const struct part *restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / hydro_get_mass(p);
}

/**
 * @brief Get comoving metal density from a specific metal specie.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_comoving_metal_density(const struct part *restrict p, int metal) {
  return chemistry_get_metal_mass_fraction(p, metal) *
         hydro_get_comoving_density(p);
}

/**
 * @brief Get the physical metal density from a specific metal specie.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_physical_metal_density(const struct part *restrict p, int metal,
                                     const struct cosmology *cosmo) {
  return cosmo->a3_inv * chemistry_get_comoving_metal_density(p, metal);
}

/**
 * @brief Get metal mass from a specific metal specie.
 *
 * This function sets the metal mass to 0 if metal_mass negative to prevent
 * contaminating other modules.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_part_corrected_metal_mass(const struct part *restrict p,
                                        int metal) {
  /* Clamp negative values to 0.0. We tolerate negative values in the diffusion
     solver, not in the physics module */
  double mZ = max(p->chemistry_data.metal_mass[metal], 0.0);
  const double m = hydro_get_mass(p);

  if (mZ >= m) {
    /* Reduce it */
    mZ = 0.01*m;
  }
  return mZ;
}

/**
 * @brief Get the traceless physical shear tensor.
 *
 * @param p Particle.
 * @param cosmo The current cosmological model.
 * @param S (return) Pointer to a 3x3 matrix.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_physical_shear_tensor(const struct part *restrict p,
                                    const struct cosmology *cosmo,
                                    double S[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      S[i][j] = 0.5 * (p->chemistry_data.filtered.grad_v_tilde[i][j] +
                       p->chemistry_data.filtered.grad_v_tilde[j][i]);

      /* Convert to physical units.
         In Swift, v_peculiar = v_c / a ; grad_p = a^{-1} grad_c. */
      S[i][j] *= cosmo->a2_inv;
    }
  }

  /* Now add the Hubble flow term to the diagonal elements. Note that because
     the shear tensore is symmetric, the Hubble flow adds only once. */
  S[0][0] += cosmo->H;
  S[1][1] += cosmo->H;
  S[2][2] += cosmo->H;

  /* The trace encode volume changes (compression and expansion). We do not
     consider these are turbulences so we remove them. */
  const double trace = S[0][0] + S[1][1] + S[2][2];
  S[0][0] -= trace / 3.0;
  S[1][1] -= trace / 3.0;
  S[2][2] -= trace / 3.0;
}

/**
 * @brief Regularize the shear tensor as described in Balarac et al. (2013) to
 * be positive definite
 *
 * We want the diffusion to *smooth* out metals:
 *             d (rho Z)/dt = - div F, F = - K grad Z.
 * So K must be positive definite. If K is negative definite, it will flip
 * the sign of the flux and thus the diffusion *sharpens* metals instead of
 * smoothing them. The flux must always points down the gradient to smooth.
 *
 * TODO: Check this function. Add a unit test for the diagonalisation.
 *
 * @param S (return) Pointer to a 3x3 matrix shear tensor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_regularize_shear_tensor(double S[3][3]) {
  double eigenvalues[3] = {0.0};
  double eigenvector0[3] = {0.0};
  double eigenvector1[3] = {0.0};
  double eigenvector2[3] = {0.0};

  /* Compute the eigenvalues and eigenvectors. S is symmetric by construction.
   */
  chemistry_utils_diagonalize_3x3(S, eigenvalues, eigenvector0, eigenvector1,
                                  eigenvector2);

  const double eigenvectors[3][3] = {
      {eigenvector0[0], eigenvector1[0], eigenvector2[0]},
      {eigenvector0[1], eigenvector1[1], eigenvector2[1]},
      {eigenvector0[2], eigenvector1[2], eigenvector2[2]}};
  double S_plus[3][3] = {{0.0}};

  /* Compute S_plus as the sum of max(0, lambda^(k)) * e_i^(k) * e_j^(k) */
  for (int k = 0; k < 3; k++) {
    const double lambda_k = eigenvalues[k];  // Get the k-th eigenvalue
    const double lambda_k_plus =
        fmax(0.0, lambda_k);  // Take max(0, lambda^(k))

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        const double e_ik = eigenvectors[i][k];
        const double e_jk = eigenvectors[j][k];
        S_plus[i][j] += lambda_k_plus * e_ik * e_jk;
      }
    }
  }

  /* Copy S_plus back into S (overwriting it) */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      S[i][j] = S_plus[i][j];
    }
  }
}

/**
 * @brief Get the physical diffusion matrix K.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 * @param K (return) Pointer to a 3x3 diffusion tensor.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_physical_matrix_K(const struct part *restrict p,
                                const struct chemistry_global_data *chem_data,
                                const struct cosmology *cosmo, double K[3][3]) {
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
    chemistry_get_physical_shear_tensor(p, cosmo, K);

    /* Now regularize the shear tensor by considering only the positive
       eigenvalues (Balarac et al. (2013)). This is now called the S_plus
       matrix. */
    chemistry_regularize_shear_tensor(K);

    /* K = kappa * S_plus */
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
 * @param K Pointer to a 3x3 diffusion tensor.
 */
__attribute__((always_inline)) INLINE static double chemistry_get_matrix_norm(
    const double K[3][3]) {
  double norm = 0.0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      norm += K[i][j] * K[i][j];
    }
  }
  return sqrt(norm);
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
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_diffusion_coefficient(
    struct part *restrict p, const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {
  /* Convert density to physical units. */
  const float rho = p->chemistry_data.filtered.rho * cosmo->a3_inv;

  /* Convert smoothing length to physical units */
  const double h2_p = cosmo->a * cosmo->a * p->h * p->h;

  if (chem_data->diffusion_mode == isotropic_constant) {
    return chem_data->diffusion_coefficient;
  } else if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Get the physical shear tensor */
    double S[3][3];
    chemistry_get_physical_shear_tensor(p, cosmo, S);

    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho *
           chemistry_get_matrix_norm(S);
  } else {
    /* Note that this is multiplied by the matrix S to get the full matrix K */
    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho;
  }
}

/**
 * @brief Get the gradients of metal mass fraction a given metal specie.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dF Metal mass fraction gradient (of size 3).
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_metal_mass_fraction_gradients(const struct part *restrict p,
                                            int metal, double dF[3]) {
  dF[0] = p->chemistry_data.gradients.Z[metal][0];
  dF[1] = p->chemistry_data.gradients.Z[metal][1];
  dF[2] = p->chemistry_data.gradients.Z[metal][2];
}

/**
 * @brief Get the gradients of metal density a given metal specie.
 *
 * Get grad U = grad rho_Z.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dF Metal mass density gradient (of size 3).
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_metal_density_gradients(const struct part *restrict p, int metal,
                                      double dF[3]) {
  dF[0] = p->chemistry_data.gradients.rhoZ[metal][0];
  dF[1] = p->chemistry_data.gradients.rhoZ[metal][1];
  dF[2] = p->chemistry_data.gradients.rhoZ[metal][2];
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
 * @brief Compute the parabolic diffusion flux of given metal specie.
 *
 * F_diss = - K * \nabla \otimes q
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param F_diff (return) Array to write diffusion flux component into.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_physical_parabolic_flux(
    const struct part *restrict p, int metal, double F_diff[3],
    const struct chemistry_global_data *chem_data,
    const struct cosmology *cosmo) {

  /* In physical units */
  const double kappa = p->chemistry_data.kappa;

  /* The gradient needs to be converted to physical units:
                                 grad_p = a^{-1} * grad_c.
     The metallicity is already physical (Z_p = Z_c = Z). */

  if (chem_data->diffusion_mode == isotropic_constant) {
    /* Here, we use grad (rho Z) */
    const double a_inv_4 =
        cosmo->a_inv * cosmo->a_inv * cosmo->a_inv * cosmo->a_inv;

    /* a^-3 for density and a^-1 for gradient */
    F_diff[0] = -kappa * p->chemistry_data.gradients.rhoZ[metal][0] * a_inv_4;
    F_diff[1] = -kappa * p->chemistry_data.gradients.rhoZ[metal][1] * a_inv_4;
    F_diff[2] = -kappa * p->chemistry_data.gradients.rhoZ[metal][2] * a_inv_4;
  } else if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Here, we use grad Z */
    F_diff[0] = -kappa * p->chemistry_data.gradients.Z[metal][0] * cosmo->a_inv;
    F_diff[1] = -kappa * p->chemistry_data.gradients.Z[metal][1] * cosmo->a_inv;
    F_diff[2] = -kappa * p->chemistry_data.gradients.Z[metal][2] * cosmo->a_inv;
  } else {
    /* Here, we use grad Z */
    /* Initialise to the flux to 0 */
    F_diff[0] = 0.0;
    F_diff[1] = 0.0;
    F_diff[2] = 0.0;

    /* Compute diffusion matrix K */
    double K[3][3];
    chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        F_diff[i] -=
            K[i][j] * p->chemistry_data.gradients.Z[metal][j] * cosmo->a_inv;
      }
    } /* End of matrix multiplication */
  } /* end of if else diffusion_mode */
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction(const struct part *restrict p) {
  const float m = hydro_get_mass(p);
  float m_Z_tot = 0.0;
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    m_Z_tot += chemistry_get_part_corrected_metal_mass(p, i);
  }

  /* If this really happens, reduce the mass. */
  if (m_Z_tot > m) {
    m_Z_tot = 0.9*m;
  }

  return m_Z_tot / m;
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * This is unused in GEAR MF diffusion. --> return NULL
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const *
chemistry_get_metal_mass_fraction_for_feedback(const struct part *restrict p) {
  error("Not implemented");
  return NULL;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction_for_feedback(
    const struct part *restrict p) {
  return chemistry_get_total_metal_mass_fraction(p);
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the star
 * particle to be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_star_total_metal_mass_fraction_for_feedback(
    const struct spart *restrict sp) {
  float Z_tot = 0.0;
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    Z_tot += sp->chemistry_data.metal_mass_fraction[i];
  }
  return Z_tot;
}

/**
 * @brief Returns the total iron mass fraction of the star particle to be used
 * in feedback/enrichment related routines.
 *
 * We assume iron to be stored at index 0.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_star_total_iron_mass_fraction_for_feedback(
    const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction[0];
}

/**
 * @brief Returns the total iron mass fraction of the sink particle to be used
 * in feedback/enrichment related routines.
 *
 * We assume iron to be stored at index 0.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_sink_total_iron_mass_fraction_for_feedback(
    const struct sink *restrict sink) {

  return sink->chemistry_data.metal_mass_fraction[0];
}

/**
 * @brief Returns the abundances (metal mass fraction) of the star particle to
 * be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const *
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart *restrict sp) {

  return sp->chemistry_data.metal_mass_fraction;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the gas
 * particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_total_metal_mass_fraction_for_cooling(
    const struct part *restrict p) {
  return chemistry_get_total_metal_mass_fraction(p);
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the gas
 * particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const *
chemistry_get_metal_mass_fraction_for_cooling(const struct part *restrict p) {
  error("This function is not used in GEAR");
  return p->chemistry_data.metal_mass;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the gas
 * particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_total_metal_mass_fraction_for_star_formation(
    const struct part *restrict p) {
  return chemistry_get_total_metal_mass_fraction(p);
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the gas
 * particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const *
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part *restrict p) {
  error("This function is not used in GEAR");
  return p->chemistry_data.metal_mass;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the gas
 * particle to be used in the stats related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_for_stats(const struct part *restrict p) {
  return chemistry_get_total_metal_mass_fraction(p) * hydro_get_mass(p);
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the star
 * particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart *restrict sp) {
  float Z_tot = 0.0;
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    Z_tot += sp->chemistry_data.metal_mass_fraction[i];
  }
  return Z_tot * sp->mass;
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the black hole
 * particle to be used in the stats related routines.
 *
 * @param bp Pointer to the BH particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart *restrict bp) {
  error("No BH yet in GEAR");
  return 0.f;
}

/* Import the right chemistry_part_data definition */
#if defined(CHEMISTRY_GEAR_MF_HYPERBOLIC_DIFFUSION)
#include "hyperbolic/chemistry_getters.h"
#endif

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GETTERS_H  */
