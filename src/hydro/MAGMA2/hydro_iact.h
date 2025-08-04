/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk) &
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Doug Rennehan (douglas.rennehan@gmail.com)
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
#ifndef SWIFT_MAGMA2_HYDRO_IACT_H
#define SWIFT_MAGMA2_HYDRO_IACT_H

/**
 * @file MAGMA2/hydro_iact.h
 * @brief Density-Energy non-conservative implementation of SPH,
 *        with added MAGMA2 physics (Rosswog 2020) (interaction routines)
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Kernel weights to be filled */
  float wi, wj, wi_dx, wj_dx;

  const hydro_real_t r = sqrt(r2);

  /* Get the masses. */
  const hydro_real_t mi = pi->mass;
  const hydro_real_t mj = pj->mass;

  /* Compute density of pi. */
  const hydro_real_t hi_inv = 1. / hi;
  const float xi = r * hi_inv;

  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + xi * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  adaptive_softening_add_correction_term(pi, xi, hi_inv, mj);

  /* Compute density of pj. */
  const hydro_real_t hj_inv = 1. / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + xj * wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + xj * wj_dx);

  adaptive_softening_add_correction_term(pj, xj, hj_inv, mi);

  /* Now we need to compute the derivative terms */
  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;
  const hydro_real_t faci = mj * wi_dx * r_inv;
  const hydro_real_t facj = mi * wj_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */
  const hydro_real_t du = pi->u - pj->u;
  pi->gradients.u_aux[0] += du * dx[0] * faci;
  pi->gradients.u_aux[1] += du * dx[1] * faci;
  pi->gradients.u_aux[2] += du * dx[2] * faci;

  pj->gradients.u_aux[0] += du * dx[0] * facj;
  pj->gradients.u_aux[1] += du * dx[1] * facj;
  pj->gradients.u_aux[2] += du * dx[2] * facj;

  pi->gradients.u_aux_norm[0] += dx[0] * dx[0] * faci;
  pi->gradients.u_aux_norm[1] += dx[1] * dx[1] * faci;
  pi->gradients.u_aux_norm[2] += dx[2] * dx[2] * faci;

  pj->gradients.u_aux_norm[0] += dx[0] * dx[0] * facj;
  pj->gradients.u_aux_norm[1] += dx[1] * dx[1] * facj;
  pj->gradients.u_aux_norm[2] += dx[2] * dx[2] * facj;

  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * faci;
      pj->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * facj;

      pi->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * faci;
      pj->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * facj;
    }
  }

#ifdef MAGMA2_DEBUG_CHECKS
  /* Number of neighbors */
  pi->num_ngb++;
  pj->num_ngb++;
#endif

#ifdef hydro_props_use_adiabatic_correction
  /* Needed for the adiabatic kernel correction factor */
  pi->gradients.adiabatic_f_numerator += mj * r2 * wi;
  pj->gradients.adiabatic_f_numerator += mi * r2 * wj;
#endif

}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {

  /* Kernel weights to be filled */
  float wi, wi_dx;

  /* Get the masses. */
  const hydro_real_t mj = pj->mass;

  /* Get r and r inverse. */
  const hydro_real_t r = sqrt(r2);

  const hydro_real_t h_inv = 1. / hi;
  const float xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + xi * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + xi * wi_dx);

  adaptive_softening_add_correction_term(pi, xi, h_inv, mj);

  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;
  const hydro_real_t faci = mj * wi_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */

  const hydro_real_t du = pi->u - pj->u;
  pi->gradients.u_aux[0] += du * dx[0] * faci;
  pi->gradients.u_aux[1] += du * dx[1] * faci;
  pi->gradients.u_aux[2] += du * dx[2] * faci;

  pi->gradients.u_aux_norm[0] += dx[0] * dx[0] * faci;
  pi->gradients.u_aux_norm[1] += dx[1] * dx[1] * faci;
  pi->gradients.u_aux_norm[2] += dx[2] * dx[2] * faci;

  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * faci;
      pi->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * faci;
    }
  }

#ifdef MAGMA2_DEBUG_CHECKS
  /* Neighbour number */
  pi->num_ngb++;
#endif

#ifdef hydro_props_use_adiabatic_correction
  /* Needed for the adiabatic kernel correction factor */
  pi->gradients.adiabatic_f_numerator += mj * r2 * wi;
#endif

}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Get particle properties */
  const hydro_real_t mi = hydro_get_mass(pi);
  const hydro_real_t mj = hydro_get_mass(pj);

  const hydro_real_t rhoi = hydro_get_comoving_density(pi);
  const hydro_real_t rhoj = hydro_get_comoving_density(pj);

  const hydro_real_t rhoi_inv = 1. / rhoi;
  const hydro_real_t rhoj_inv = 1. / rhoj;

  const hydro_real_t r = sqrt(r2);

  float wi, wi_dx, wj, wj_dx;

  const float xi = r / hi;
  const float xj = r / hj;

  kernel_deval(xi, &wi, &wi_dx);
  kernel_deval(xj, &wj, &wj_dx);

  const hydro_real_t faci = mj * rhoj_inv * wi;
  const hydro_real_t facj = mi * rhoi_inv * wj;
  
  /* Compute all of the first-order gradients, and second-order gradients */

  /* Rosswog 2020 Equation 18 gradients. In the paper he uses (vj - vi) and
   * (rj - ri), however this is symmetric so no sign problems. */

  /* Internal energy gradient */
  const hydro_real_t du = pi->u - pj->u;

  pi->gradients.u[0] += du * dx[0] * faci;
  pi->gradients.u[1] += du * dx[1] * faci;
  pi->gradients.u[2] += du * dx[2] * faci;

  pj->gradients.u[0] += du * dx[0] * facj;
  pj->gradients.u[1] += du * dx[1] * facj;
  pj->gradients.u[2] += du * dx[2] * facj;

  /* Velocity gradients */
  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};


  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      const hydro_real_t du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[k];
      pi->gradients.u_hessian[k][i] += du_k * dx[i] * faci;
      pj->gradients.u_hessian[k][i] += du_k * dx[i] * facj;

      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.correction_matrix[k][i] += dx[k] * dx[i] * faci;
      pj->gradients.correction_matrix[k][i] += dx[k] * dx[i] * facj;

      /* Indices in Rosswog 2020 are i for dv and k for dx. In this loop,
       * they are swapped just because correction_matrix is computed with
       * the paper indices. */
      pi->gradients.velocity_tensor[k][i] += dv[k] * dx[i] * faci;
      pj->gradients.velocity_tensor[k][i] += dv[k] * dx[i] * facj;

      const hydro_real_t dv_grad_ki = pi->gradients.velocity_tensor_aux[k][i] -
                                      pj->gradients.velocity_tensor_aux[k][i];

      /* Equation 19 indices:
       * Index i: velocity direction
       * Index k: gradient direction
       *
       * Our indices:
       * Index k: velocity direction
       * Index i: gradient direction
       * Index j: second derivative gradient direction
       */
      for (int j = 0; j < 3; j++) {
        pi->gradients.velocity_hessian[k][i][j] += dv_grad_ki * dx[j] * faci;
        pj->gradients.velocity_hessian[k][i][j] += dv_grad_ki * dx[j] * facj;
      }
    }
  }

#ifdef hydro_props_use_adiabatic_correction
  /* Correction terms for div v */
  pi->gradients.adiabatic_f_denominator += mj * rhoj_inv * r2 * wi;
  pj->gradients.adiabatic_f_denominator += mi * rhoi_inv * r2 * wj;
#endif

}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Get particle properties */
  const hydro_real_t mj = hydro_get_mass(pj);
  const hydro_real_t rhoj = hydro_get_comoving_density(pj);
  const hydro_real_t rhoj_inv = 1. / rhoj;
  const hydro_real_t r = sqrt(r2);
  float wi, wi_dx;
  const float xi = r / hi;
  kernel_deval(xi, &wi, &wi_dx);
  const hydro_real_t faci = mj * rhoj_inv * wi;

  /* Compute all of the first-order gradients, and second-order gradients */

  /* Rosswog 2020 Equation 18 gradients. In the paper he uses (vj - vi) and
   * (rj - ri), however this is symmetric so no sign problems. */

  /* Internal energy gradient */
  const hydro_real_t du = pi->u - pj->u;

  pi->gradients.u[0] += du * dx[0] * faci;
  pi->gradients.u[1] += du * dx[1] * faci;
  pi->gradients.u[2] += du * dx[2] * faci;

  /* Velocity gradients */
  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};


  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      const hydro_real_t du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[i];
      pi->gradients.u_hessian[k][i] += du_k * dx[i] * faci;

      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.correction_matrix[k][i] += dx[k] * dx[i] * faci;

      /* Indices in Rosswog 2020 are i for dv and k for dx. In this loop,
       * they are swapped just because correction_matrix is computed with
       * the paper indices. */
      pi->gradients.velocity_tensor[k][i] += dv[k] * dx[i] * faci;

      const hydro_real_t dv_grad_ki = pi->gradients.velocity_tensor_aux[k][i] - 
                                      pj->gradients.velocity_tensor_aux[k][i];

      /* Equation 19 indices:
       * Index i: velocity direction
       * Index k: gradient direction
       *
       * Our indices:
       * Index k: velocity direction
       * Index i: gradient direction
       * Index j: second derivative gradient direction
       */
      for (int j = 0; j < 3; j++) {
        pi->gradients.velocity_hessian[k][i][j] += dv_grad_ki * dx[j] * faci;
      }
    }
  }

#ifdef hydro_props_use_adiabatic_correction
  /* Correction terms for div v */
  pi->gradients.adiabatic_f_denominator += mj * rhoj_inv * r2 * wi;
#endif

}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, struct part* restrict pj, const float a,
    const float H) {

  /* Cosmological factors entering the EoMs */
  const hydro_real_t fac_mu = pow_three_gamma_minus_five_over_two(a);
  const hydro_real_t a_Hubble = a * H;
  const hydro_real_t a2_Hubble = a * a_Hubble;

  const hydro_real_t r = sqrt(r2);
  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const hydro_real_t mj = pj->mass;
  const hydro_real_t mi = pi->mass;

  const hydro_real_t rhoi = pi->rho;
  const hydro_real_t rhoj = pj->rho;
  const hydro_real_t rhoi_inv = 1. / rhoi;
  const hydro_real_t rhoj_inv = 1. / rhoj;
  const hydro_real_t rhoij_inv = rhoi_inv * rhoj_inv;

  const hydro_real_t pressurei = pi->force.pressure;
  const hydro_real_t pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const hydro_real_t hi_inv = 1. / hi;
  const hydro_real_t hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const hydro_real_t hj_inv = 1. / hj;
  const hydro_real_t hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* For dh/dt and fall-back SPH */
  const hydro_real_t hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
  const hydro_real_t hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
  const hydro_real_t wi_dr = hi_inv_dim_plus_one * wi_dx;
  const hydro_real_t wj_dr = hj_inv_dim_plus_one * wj_dx;

  /* Peculiar velocity difference vector */
  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  hydro_real_t G_ij[3] = {0., 0., 0.};

  /* These are set whether or not we fall back onto SPH gradients */
  hydro_real_t visc_acc_term = 0.;
  hydro_real_t sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  hydro_real_t sph_du_term_i = pressurei * rhoij_inv;
  hydro_real_t sph_du_term_j = pressurej * rhoij_inv;
  hydro_real_t visc_du_term = 0.;
  hydro_real_t cond_du_term = 0.;

  /* Correction factor must be well-conditioned. */
  const char C_well_conditioned = 
      (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned);
  /* First order density-independent velocity gradients */
  const char D_well_conditioned = 
      (pi->gradients.D_well_conditioned && pj->gradients.D_well_conditioned);

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  if (C_well_conditioned && D_well_conditioned) {

    /* Separation vectors and swapped */
    const hydro_real_t dx_ij[3] = {dx[0], dx[1], dx[2]};
    const hydro_real_t dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* Corrected gradients */
    hydro_real_t G_i[3] = {0., 0., 0.};
    hydro_real_t G_j[3] = {0., 0., 0.};
    /* Rosswog 2020 Eqs 4 & 5 use dx_ji for both */
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx_ji, G_i);
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx_ji, G_j);

    G_i[0] *= wi * hi_inv_dim;
    G_i[1] *= wi * hi_inv_dim;
    G_i[2] *= wi * hi_inv_dim;

    G_j[0] *= wj * hj_inv_dim;
    G_j[1] *= wj * hj_inv_dim;
    G_j[2] *= wj * hj_inv_dim;

    /* Compute second order reconstruction of velocity between pi & pj */
    hydro_real_t vi_reconstructed[3] = {0., 0., 0.};
    hydro_real_t vj_reconstructed[3] = {0., 0., 0.};

    /* Important: Rosswog 2020 h_i and h_j are without kernel_gamma. Therefore,
     * use xi and xj in the slope limiting procedure. */

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const hydro_real_t phi_ij_vec = 
        hydro_vector_van_leer_phi(pi->gradients.velocity_tensor, 
                                  pj->gradients.velocity_tensor,
                                  dx_ij, xi, xj);

    /* Need these recast in case of switching precision */
    const hydro_real_t v_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
    const hydro_real_t v_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ji, v_i, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ij, v_j, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);


    /* Artificial viscosity */


    const hydro_real_t dv_ij[3] = {vi_reconstructed[0] - vj_reconstructed[0],
                                   vi_reconstructed[1] - vj_reconstructed[1],
                                   vi_reconstructed[2] - vj_reconstructed[2]};

    hydro_real_t mu_i = 0.;
    hydro_real_t mu_j = 0.;

    /* Get the acceleration term, depends on the weighting scheme */
    visc_acc_term = 
        hydro_get_visc_acc_term_and_mu(pi, pj, dv_ij, dx_ij, r2, r_inv, 
                                       fac_mu, a2_Hubble, &mu_i, &mu_j);

    /* Split heating between the two particles */
    visc_du_term = 0.5 * visc_acc_term;

    /* Viscous signal velocity */
    const hydro_real_t v_sig_visc = 
        hydro_get_visc_signal_velocity(dx, pi, pj, mu_i, mu_j, 
                                       const_viscosity_beta);


    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;

    /* First order internal energy gradients are well-conditioned */
    if (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned) {
      /* Compute global Van Leer limiter (scalar, not component-wise) */
      const hydro_real_t phi_ij_scalar = 
          hydro_scalar_van_leer_phi(pi->gradients.u, 
                                    pj->gradients.u,
                                    dx_ij, xi, xj);

      /* dx_ji for particle i and dx_ij for particle j */
      hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ji, 
                                               (hydro_real_t)pi->u, 
                                               pi->gradients.u,
                                               pi->gradients.u_hessian,
                                               &ui_reconstructed);

      hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ij, 
                                               (hydro_real_t)pj->u, 
                                               pj->gradients.u,
                                               pj->gradients.u_hessian,
                                               &uj_reconstructed);
    }

    const hydro_real_t rho_ij = 0.5 * (rhoi + rhoj);
    const hydro_real_t rho_ij_inv = 1. / rho_ij;
    const hydro_real_t dv_Hubble[3] = {dv[0] + a_Hubble * dx[0],
                                       dv[1] + a_Hubble * dx[1],
                                       dv[2] + a_Hubble * dx[2]};
    const hydro_real_t dv_Hubble_norm = hydro_vec3_norm(dv_Hubble);

    /* Signal velocity is the sum of pressure and speed contributions */
    const hydro_real_t v_sig_speed = fac_mu * dv_Hubble_norm;
    const hydro_real_t v_sig_pressure = 
        sqrt(fabs(pressurei - pressurej) * rho_ij_inv);
    const hydro_real_t v_sig_cond = 0.5 * (v_sig_speed + v_sig_pressure);

    const hydro_real_t du_ij = ui_reconstructed - uj_reconstructed;

    /* Add conductivity to the specific energy */
    cond_du_term = const_conductivity_alpha * v_sig_cond * du_ij * rho_ij_inv;


    /* Finalize everything with the correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    hydro_get_average_kernel_gradient(pi, pj, G_i, G_j, G_ij);

    /* For conduction term */
    const hydro_real_t G_ij_norm = hydro_vec3_norm(G_ij);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const hydro_real_t dv_dot_G_ij = hydro_vec3_vec3_dot(dv, G_ij);

    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij;
    cond_du_term *= -G_ij_norm; /* Eq. 24 Rosswog 2020 */


    /* Get the time derivative for h. */
    

    pi->force.h_dt -= dv_dot_G_ij;
    pj->force.h_dt -= dv_dot_G_ij;


    /* Timestepping */


    /* New signal velocity */
    const hydro_real_t v_sig_max = max(v_sig_visc, v_sig_cond);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, v_sig_max);
    pj->v_sig_max = max(pj->v_sig_max, v_sig_max);

    /* Average softening in kernel */
    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = min(pi->h_min, h_ij);
    pj->h_min = min(pj->h_min, h_ij);

    /* New timestep estimate */
    const hydro_real_t dt_min_i = pi->h_min / pi->v_sig_max;
    const hydro_real_t dt_min_j = pj->h_min / pj->v_sig_max;
    pi->dt_min = min(pi->dt_min, dt_min_i);
    pj->dt_min = min(pj->dt_min, dt_min_j);
  }
  else {

    /* Dot product with the distance vector to align
     * with dW/dr (dW/dr = r_ij / |r_ij| dot grad W */
    G_ij[0] = dx[0];
    G_ij[1] = dx[1];
    G_ij[2] = dx[2];

    /* Compute dv dot dr. */
    const hydro_real_t dvdr = hydro_vec3_vec3_dot(dv, G_ij);

    /* Includes the hubble flow term; not used for du/dt */
    const hydro_real_t dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const hydro_real_t omega_ij = min(dvdr_Hubble, 0.);
    /* This is 0 or negative */
    const hydro_real_t mu_ij = fac_mu * r_inv * omega_ij;

    /* Construct the full viscosity term */
    const hydro_real_t rho_ij = rhoi + rhoj;
    const hydro_real_t alpha = const_viscosity_alpha;
    const hydro_real_t visc =
        omega_ij < 0.
            ? (-0.25 * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                   mu_ij +
               const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5 * rho_ij)
            : 0.;

    /* New signal velocity */
    const hydro_real_t new_v_sig =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, new_v_sig);
    pj->v_sig_max = max(pj->v_sig_max, new_v_sig);

    /* Minimum softening in kernel */
    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = min(pi->h_min, h_ij);
    pj->h_min = min(pj->h_min, h_ij);

    /* New timestep estimate */
    const hydro_real_t dt_min_i = pi->h_min / pi->v_sig_max;
    const hydro_real_t dt_min_j = pj->h_min / pj->v_sig_max;
    pi->dt_min = min(pi->dt_min, dt_min_i);
    pj->dt_min = min(pj->dt_min, dt_min_j);

    visc_acc_term = visc;
    visc_du_term = 0.5 * visc_acc_term;

    /* Variable smoothing length term */
    const hydro_real_t f_i = pi->force.f;
    const hydro_real_t f_j = pj->force.f;
    const hydro_real_t kernel_gradient = 
        0.5 * (f_i * wi_dr + f_j * wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;

    sph_du_term_i *= dvdr * kernel_gradient;
    sph_du_term_j *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= dvdr * r_inv * wi_dr;
    pj->force.h_dt -= dvdr * r_inv * wj_dr;
  }

  /* Assemble the acceleration */
  const hydro_real_t acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  pj->a_hydro[0] += mi * acc * G_ij[0];
  pj->a_hydro[1] += mi * acc * G_ij[1];
  pj->a_hydro[2] += mi * acc * G_ij[2];

  /* Get the time derivative for u. */

  /* Assemble the energy equation term */
  const hydro_real_t du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;
  const hydro_real_t du_dt_j = sph_du_term_j + visc_du_term - cond_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of part*icle i.
 * @param hj Comoving smoothing-length of part*icle j.
 * @param pi First part*icle.
 * @param pj Second part*icle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part* restrict pi, const struct part* restrict pj, const float a,
    const float H) {

  /* Cosmological factors entering the EoMs */
  const hydro_real_t fac_mu = pow_three_gamma_minus_five_over_two(a);
  const hydro_real_t a_Hubble = a * H;
  const hydro_real_t a2_Hubble = a * a_Hubble;

  const hydro_real_t r = sqrtf(r2);
  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const hydro_real_t mj = pj->mass;

  const hydro_real_t rhoi = pi->rho;
  const hydro_real_t rhoj = pj->rho;
  const hydro_real_t rhoi_inv = 1. / rhoi;
  const hydro_real_t rhoj_inv = 1. / rhoj;
  const hydro_real_t rhoij_inv = rhoi_inv * rhoj_inv;

  const hydro_real_t pressurei = pi->force.pressure;
  const hydro_real_t pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const hydro_real_t hi_inv = 1.0 / hi;
  const hydro_real_t hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const hydro_real_t hj_inv = 1.0 / hj;
  const hydro_real_t hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* For dh/dt and fall-back SPH */
  const hydro_real_t hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
  const hydro_real_t wi_dr = hi_inv_dim_plus_one * wi_dx;
  
  /* Peculiar velocity difference vector */
  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  hydro_real_t G_ij[3] = {0., 0., 0.};

  /* These are set whether or not we fall back onto SPH gradients */
  hydro_real_t visc_acc_term = 0.;
  hydro_real_t sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  hydro_real_t sph_du_term_i = pressurei * rhoij_inv;
  hydro_real_t visc_du_term = 0.;
  hydro_real_t cond_du_term = 0.;

  /* Correction factor must be well-conditioned. */
  const char C_well_conditioned = 
      (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned);
  /* First order density-independent velocity gradients */
  const char D_well_conditioned = 
      (pi->gradients.D_well_conditioned && pj->gradients.D_well_conditioned);

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  if (C_well_conditioned && D_well_conditioned) {

    /* Separation vectors and swapped */
    const hydro_real_t dx_ij[3] = {dx[0], dx[1], dx[2]};
    const hydro_real_t dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* Corrected gradients */
    hydro_real_t G_i[3] = {0., 0., 0.};
    hydro_real_t G_j[3] = {0., 0., 0.};
    /* Rosswog 2020 Eqs 4 & 5 use dx_ji for both */
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx_ji, G_i);
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx_ji, G_j);

    G_i[0] *= wi * hi_inv_dim;
    G_i[1] *= wi * hi_inv_dim;
    G_i[2] *= wi * hi_inv_dim;

    G_j[0] *= wj * hj_inv_dim;
    G_j[1] *= wj * hj_inv_dim;
    G_j[2] *= wj * hj_inv_dim;

    /* Compute second order reconstruction of velocity between pi & pj */
    hydro_real_t vi_reconstructed[3] = {0., 0., 0.};
    hydro_real_t vj_reconstructed[3] = {0., 0., 0.};

    /* Important: Rosswog 2020 h_i and h_j are without kernel_gamma. Therefore,
     * use xi and xj in the slope limiting procedure. */

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const hydro_real_t phi_ij_vec = 
        hydro_vector_van_leer_phi(pi->gradients.velocity_tensor, 
                                  pj->gradients.velocity_tensor,
                                  dx_ij, xi, xj);

    /* Need these recast in case of switching precision */
    const hydro_real_t v_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
    const hydro_real_t v_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ji, v_i, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ij, v_j, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);


    /* Artificial viscosity */


    const hydro_real_t dv_ij[3] = {vi_reconstructed[0] - vj_reconstructed[0],
                                   vi_reconstructed[1] - vj_reconstructed[1],
                                   vi_reconstructed[2] - vj_reconstructed[2]};

    hydro_real_t mu_i = 0.;
    hydro_real_t mu_j = 0.;

    /* Get the acceleration term, depends on the weighting scheme */
    visc_acc_term = 
        hydro_get_visc_acc_term_and_mu(pi, pj, dv_ij, dx_ij, r2, r_inv, 
                                       fac_mu, a2_Hubble, &mu_i, &mu_j);
    visc_du_term = 0.5 * visc_acc_term;

        /* Viscous signal velocity */
    const hydro_real_t v_sig_visc = 
        hydro_get_visc_signal_velocity(dx, pi, pj, mu_i, mu_j, 
                                       const_viscosity_beta);


    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;

    /* First order internal energy gradients are well-conditioned */
    if (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned) {
      /* Compute global Van Leer limiter (scalar, not component-wise) */
      const hydro_real_t phi_ij_scalar = 
          hydro_scalar_van_leer_phi(pi->gradients.u, 
                                    pj->gradients.u,
                                    dx_ij, xi, xj);

      /* dx_ji for particle i and dx_ij for particle j */
      hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ji, 
                                               (hydro_real_t)pi->u, 
                                               pi->gradients.u,
                                               pi->gradients.u_hessian,
                                               &ui_reconstructed);

      hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ij, 
                                               (hydro_real_t)pj->u, 
                                               pj->gradients.u,
                                               pj->gradients.u_hessian,
                                               &uj_reconstructed);
    }

    const hydro_real_t rho_ij = 0.5 * (rhoi + rhoj);
    const hydro_real_t rho_ij_inv = 1. / rho_ij;
    const hydro_real_t dv_Hubble[3] = {dv[0] + a_Hubble * dx[0],
                                       dv[1] + a_Hubble * dx[1],
                                       dv[2] + a_Hubble * dx[2]};
    const hydro_real_t dv_Hubble_norm = hydro_vec3_norm(dv_Hubble);

    /* Signal velocity is the sum of pressure and speed contributions */
    const hydro_real_t v_sig_speed = fac_mu * dv_Hubble_norm;
    const hydro_real_t v_sig_pressure = 
        sqrt(fabs(pressurei - pressurej) * rho_ij_inv);
    const hydro_real_t v_sig_cond = 0.5 * (v_sig_speed + v_sig_pressure);

    const hydro_real_t du_ij = ui_reconstructed - uj_reconstructed;

    /* Add conductivity to the specific energy */
    cond_du_term = const_conductivity_alpha * v_sig_cond * du_ij * rho_ij_inv;


    /* Finalize the viscosity and conductivity with correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    hydro_get_average_kernel_gradient(pi, pj, G_i, G_j, G_ij);

    const hydro_real_t G_ij_norm = hydro_vec3_norm(G_ij);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const hydro_real_t dv_dot_G_ij = hydro_vec3_vec3_dot(dv, G_ij);                    

    sph_du_term_i *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij;
    cond_du_term *= -G_ij_norm; /* Eq. 24 Rosswog 2020 */


    /* Get the time derivative for h. */


    pi->force.h_dt -= dv_dot_G_ij;


    /* Timestepping */


    /* New signal velocity */
    const hydro_real_t v_sig_max = max(v_sig_visc, v_sig_cond);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, v_sig_max);

    /* Compute new timestep */
    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = min(pi->h_min, h_ij);

    const hydro_real_t dt_min_i = pi->h_min / pi->v_sig_max;
    pi->dt_min = min(pi->dt_min, dt_min_i);
  }
  else {
    
    /* Fallback SPH points in the direction between the particles. */
    G_ij[0] = dx[0];
    G_ij[1] = dx[1];
    G_ij[2] = dx[2];

    /* Compute dv dot dr. */
    const hydro_real_t dvdr = hydro_vec3_vec3_dot(dv, G_ij);

    /* Includes the hubble flow term; not used for du/dt */
    const hydro_real_t dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const hydro_real_t omega_ij = min(dvdr_Hubble, 0.);
    /* This is 0 or negative */
    const hydro_real_t mu_ij = fac_mu * r_inv * omega_ij;

    /* Construct the full viscosity term */
    const hydro_real_t rho_ij = rhoi + rhoj;
    const hydro_real_t alpha = const_viscosity_alpha;
    const hydro_real_t visc =
        omega_ij < 0.
            ? (-0.25 * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                  mu_ij +
              const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5 * rho_ij)
            : 0.;

    /* New signal velocity */
    const hydro_real_t new_v_sig =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, new_v_sig);

    /* Minimum softening in kernel */
    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = min(pi->h_min, h_ij);

    /* New time-step estimate */
    const hydro_real_t dt_min_i = pi->h_min / pi->v_sig_max;
    pi->dt_min = min(pi->dt_min, dt_min_i);


    visc_acc_term = visc;
    visc_du_term = 0.5 * visc_acc_term;

    const hydro_real_t hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const hydro_real_t wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const hydro_real_t f_i = pi->force.f;
    const hydro_real_t f_j = pj->force.f;
    const hydro_real_t kernel_gradient = 
        0.5 * (f_i * wi_dr + f_j * wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= dvdr * r_inv * wi_dr;
  }

  /* Assemble the acceleration */
  const hydro_real_t acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  /* Assemble the energy equation term */
  const hydro_real_t du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

}

#endif /* SWIFT_MAGMA2_HYDRO_IACT_H */
