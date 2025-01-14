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
 * @brief Get comoving metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_comoving_metal_density(const struct part* restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / p->geometry.volume;
}

/**
 * @brief Get the physical metal density from a specific metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_physical_metal_density(const struct part* restrict p, int metal,
                                     const struct cosmology* cosmo) {
  return cosmo->a3_inv * chemistry_get_comoving_metal_density(p, metal);
}

/**
 * @brief Get metal mass fraction from a specific metal specie.
 *
 * @param p Particle.
 * @param metal Index of metal specie
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_metal_mass_fraction(const struct part* restrict p, int metal) {
  return p->chemistry_data.metal_mass[metal] / hydro_get_mass(p);
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
__attribute__((always_inline)) INLINE static float
chemistry_get_comoving_density(const struct part* restrict p) {
  float rho = hydro_get_comoving_density(p);

  if (rho == 0.0) {
    const float r_cubed =
        kernel_gamma * kernel_gamma * kernel_gamma * p->h * p->h * p->h;
    const float volume = 4.0 / 3.0 * M_PI * r_cubed;
    rho = hydro_get_mass(p) / volume;
  }
  return rho;
}

/**
 * @brief Get the physical shear tensor.
 *
 * @param p Particle.
 * @param cosmo The current cosmological model.
 * @param S (return) Pointer to a 3x3 matrix.
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_physical_shear_tensor(const struct part* restrict p,
                                    const struct cosmology* cosmo,
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
  double S_minus[3][3] = {{0.0}};

  /* Compute S_minus as the sum of min(0, lambda^(k)) * e_i^(k) * e_j^(k) */
  for (int k = 0; k < 3; k++) {
    const double lambda_k = eigenvalues[k];  // Get the k-th eigenvalue
    const double lambda_k_minus =
        fmin(0.0, lambda_k);  // Take min(0, lambda^(k))

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        const double e_ik = eigenvectors[i][k];
        const double e_jk = eigenvectors[j][k];
        S_minus[i][j] += lambda_k_minus * e_ik * e_jk;
      }
    }
  }

  /* Copy S_minus back into S (overwriting it) */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      S[i][j] = S_minus[i][j];
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
chemistry_get_physical_matrix_K(const struct part* restrict p,
                                const struct chemistry_global_data* chem_data,
                                const struct cosmology* cosmo, double K[3][3]) {
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

    /* Now regularize the shear tensor by considering only the negative
       eigenvalues (Balarac et al. (2013)). This is now called the S_minus
       matrix. */
    chemistry_regularize_shear_tensor(K);

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
 * @param K Pointer to a 3x3 diffusion tensor.
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
chemistry_compute_diffusion_coefficient(
    struct part* restrict p, const struct chemistry_global_data* chem_data,
    const struct cosmology* cosmo) {

  float rho = p->chemistry_data.filtered.rho;

  /* In case the filtered density is 0, e.g. during the fake-timestep,
     approximate the density */
  if (rho == 0.0) {
    rho = chemistry_get_comoving_density(p);
  }

  /* Convert density to physical units. */
  rho *= cosmo->a3_inv;

  /* Convert smoothing length to physical units */
  const double h2_p = cosmo->a * cosmo->a * p->h * p->h;

  if (chem_data->diffusion_mode == isotropic_constant) {
    return chem_data->diffusion_coefficient;
  } else if (chem_data->diffusion_mode == isotropic_smagorinsky) {
    /* Get the physical shear tensor */
    double S[3][3];
    chemistry_get_physical_shear_tensor(p, cosmo, S);

    /* In the smagorinsky model, we remove the trace from S */
    const double trace = S[0][0] + S[1][1] + S[2][2];

    S[0][0] -= trace;
    S[1][1] -= trace;
    S[2][2] -= trace;

    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho *
           chemistry_get_matrix_norm(S);
  } else {
    /* Note that this is multiplied by the matrix S to get the full matrix K */
    return chem_data->diffusion_coefficient * kernel_gamma2 * h2_p * rho;
  }
}

/**
 * @brief Get the gradients of metal mass density a given metal group.
 *
 * Get grad U = grad rho_Z.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dF Metal mass fraction gradient (of size 3).
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_metal_mass_fraction_gradients(const struct part* restrict p,
                                            int metal, double dF[3]) {
  dF[0] = p->chemistry_data.gradients.Z[metal][0];
  dF[1] = p->chemistry_data.gradients.Z[metal][1];
  dF[2] = p->chemistry_data.gradients.Z[metal][2];
}

/**
 * @brief Get the gradients of metal mass fraction a given metal group.
 *
 * @param p Particle.
 * @param metal Index of metal specie.
 * @param dF Metal mass density gradient (of size 3).
 */
__attribute__((always_inline)) INLINE static void
chemistry_get_metal_density_gradients(const struct part* restrict p, int metal,
				      double dF[3]) {

  const struct chemistry_part_data* chd = &p->chemistry_data;

  /* We have U = rho_Z and q = Z.  But we computed Grad Z and Grad rho, not
     Grad (rho*Z). However, Grad (rho*Z) = Z*Grad_rho + rho*Grad_Z */
  const double rho = hydro_get_comoving_density(p);
  const double Z = chemistry_get_metal_mass_fraction(p, metal);

  /* For isotropic diffusion, \grad U = \nabla \otimes q = \grad n_Z */
  dF[0] = chd->gradients.Z[metal][0] * rho + chd->gradients.rho[0] * Z;
  dF[1] = chd->gradients.Z[metal][1] * rho + chd->gradients.rho[1] * Z;
  dF[2] = chd->gradients.Z[metal][2] * rho + chd->gradients.rho[2] * Z;
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
    const struct part* restrict p, float dvx[3], float dvy[3], float dvz[3]) {

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
 * @brief Get the physical hyperbolic diffusion soundspeed.
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static double
chemistry_get_physical_hyperbolic_soundspeed(const struct part* restrict p,
					     const struct chemistry_global_data* chem_data,
					     const struct cosmology* cosmo) {
#if defined(GEAR_MF_HYPERBOLIC_DIFFUSION)
  if (chem_data->diffusion_mode == isotropic_constant) {
    return chem_data->diffusion_coefficient/chem_data->tau;
  } else {
    /* const struct chemistry_part_data *chd = &p->chemistry_data; */

    /* Compute diffusion matrix K */
    /* double K[3][3]; */
    /* chemistry_get_physical_matrix_K(p, chem_data, cosmo, K); */
    /* const float norm_matrix_K = chemistry_get_matrix_norm(K); */

    /* /\* Note: The State vector is U = (rho*Z_1,rho*Z_2, ...). *\/ */
    /* float norm_U = 0.0; */
    /* float norm_nabla_q = 0.0; */

    /* /\* Compute the norms *\/ */
    /* for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) { */
    /*   norm_U += chemistry_get_physical_metal_density(p, i, cosmo) * */
    /* 	chemistry_get_physical_metal_density(p, i, cosmo); */

    /*   for (int j = 0; j < 3; j++) { */
    /* 	/\* Compute the Frobenius norm of \nabla \otimes q *\/ */
    /* 	norm_nabla_q += chd->gradients.Z[i][j] * chd->gradients.Z[i][j]; */
    /*   } */
    /* } */

    /* /\* Take the sqrt and convert to physical units *\/ */
    /* norm_U = sqrtf(norm_U); */
    /* norm_nabla_q = sqrtf(norm_nabla_q) * cosmo->a_inv; */

    /* /\* Prevent pathological cases *\/ */
    /* if (norm_U == 0.0) { */
    /*   /\* Limit by the soundspeed *\/ */
    /*   return hydro_get_physical_soundspeed(p, cosmo); */
    /* } */

    /* return norm_matrix_K * norm_nabla_q/norm_U ; */
    return hydro_get_physical_soundspeed(p, cosmo);
  }
#else
  return hydro_get_physical_soundspeed(p, cosmo);
#endif
}

/**
 * @brief Get the physical hyperbolic diffusion relaxation time.
 *
 * @param p Particle.
 * @param chem_data The global properties of the chemistry scheme.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
chemistry_compute_physical_tau(const struct part* restrict p,
			       const struct chemistry_global_data* chem_data,
			       const struct cosmology* cosmo) {
#if defined(GEAR_MF_HYPERBOLIC_DIFFUSION)
  if (chem_data->diffusion_mode != isotropic_constant) {
    /* Compute the diffusion matrix K */
    double K[3][3];
    chemistry_get_physical_matrix_K(p, chem_data, cosmo, K);
    const float norm_matrix_K = chemistry_get_matrix_norm(K);

    /* Get soundspeed */
    const double c_hyp = chemistry_get_physical_hyperbolic_soundspeed(p, chem_data, cosmo);
    return  norm_matrix_K / (c_hyp*c_hyp);
  } else {
    return chem_data->tau;
  }
#else
  return 0.0;
#endif
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_total_metal_mass_fraction(const struct part* restrict p) {
  float m_Z_tot = 0.0;
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    m_Z_tot += p->chemistry_data.metal_mass[i];
  }
  return m_Z_tot / hydro_get_mass(p);
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the
 * gas particle to be used in feedback/enrichment related routines.
 *
 * This is unused in GEAR MF diffusion. --> return NULL
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float const*
chemistry_get_metal_mass_fraction_for_feedback(const struct part* restrict p) {
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
    const struct part* restrict p) {
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
    const struct spart* restrict sp) {
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
    const struct spart* restrict sp) {

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
    const struct sink* restrict sink) {

  return sink->chemistry_data.metal_mass_fraction[0];
}

/**
 * @brief Returns the abundances (metal mass fraction) of the star particle to
 * be used in feedback/enrichment related routines.
 *
 * @param sp Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_star_metal_mass_fraction_for_feedback(
    const struct spart* restrict sp) {

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
    const struct part* restrict p) {
  return chemistry_get_total_metal_mass_fraction(p);
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the gas
 * particle to be used in cooling related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_metal_mass_fraction_for_cooling(const struct part* restrict p) {
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
    const struct part* restrict p) {
  return chemistry_get_total_metal_mass_fraction(p);
}

/**
 * @brief Returns the abundance array (metal mass fractions) of the gas
 * particle to be used in star formation related routines.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static double const*
chemistry_get_metal_mass_fraction_for_star_formation(
    const struct part* restrict p) {
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
chemistry_get_total_metal_mass_for_stats(const struct part* restrict p) {
  return chemistry_get_total_metal_mass_fraction(p) * hydro_get_mass(p);
}

/**
 * @brief Returns the total metallicity (metal mass fraction) of the star
 * particle to be used in the stats related routines.
 *
 * @param sp Pointer to the star particle data.
 */
__attribute__((always_inline)) INLINE static float
chemistry_get_star_total_metal_mass_for_stats(const struct spart* restrict sp) {
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
chemistry_get_bh_total_metal_mass_for_stats(const struct bpart* restrict bp) {
  error("No BH yet in GEAR");
  return 0.f;
}

#endif /* SWIFT_CHEMISTRY_GEAR_MF_DIFFUSION_GETTERS_H  */
