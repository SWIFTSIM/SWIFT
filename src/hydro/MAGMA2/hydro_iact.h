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
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, hi_inv, mj);

  /* Compute density of pj. */
  const hydro_real_t hj_inv = 1. / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  adaptive_softening_add_correction_term(pj, uj, hj_inv, mi);

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

  /* Number of neighbors */
  pi->num_ngb++;
  pj->num_ngb++;

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
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

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

  /* Neighbour number */
  pi->num_ngb++;
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

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

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
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);
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
  /* First order internal energy gradients */
  const char u_well_conditioned = 
      (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned);

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

    /* Need this for viscosity and conductivity */
    hydro_real_t dv_dot_dr_Hubble = 0.;
    hydro_vec3_vec3_dot(dv_ij, dx_ij, &dv_dot_dr_Hubble);
    dv_dot_dr_Hubble += a2_Hubble * r2;

#ifdef hydro_props_use_asymmetric_viscosity_mu
    const hydro_real_t dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const hydro_real_t eta_i[3] = {dx_ij[0] * hi_inv,
                                   dx_ij[1] * hi_inv, 
                                   dx_ij[2] * hi_inv};
    const hydro_real_t eta_j[3] = {dx_ji[0] * hj_inv,
                                   dx_ji[1] * hj_inv,
                                   dx_ji[2] * hj_inv};
    hydro_real_t eta_i2 = 0.;
    hydro_vec3_vec3_dot(eta_i, eta_i, &eta_i2);

    hydro_real_t eta_j2 = 0.;
    hydro_vec3_vec3_dot(eta_j, eta_j, &eta_j2);

    hydro_real_t dv_dot_eta_i = 0.;
    hydro_vec3_vec3_dot(dv_ij, eta_i, &dv_dot_eta_i);
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const hydro_real_t dv_dot_eta_i_Hubble = 
        dv_dot_eta_i + a2_Hubble * r2 * hi_inv;

    hydro_real_t dv_dot_eta_j = 0.;
    hydro_vec3_vec3_dot(dv_ji, eta_j, &dv_dot_eta_j);
    /* Scale Hubble flow by hj_inv so it is overall scaled */
    const hydro_real_t dv_dot_eta_j_Hubble = 
        dv_dot_eta_j + a2_Hubble * r2 * hj_inv;

    /* Is the flow converging? */
    const char conv_i = (dv_dot_eta_i_Hubble < 0.) ? 1 : 0;
    const char conv_j = (dv_dot_eta_j_Hubble < 0.) ? 1 : 0;

    /* Is the flow converging? If so, multiply mu by fac_mu. If not, zero. */
    const hydro_real_t conv = (conv_i && conv_j) ? fac_mu : 0.;

#ifdef MAGMA2_DEBUG_CHECKS
    if (conv_i != conv_j) {
      warning("Flow mismatch for particles pi=%lld and pj=%lld.\n"
              "conv_i = %d, conv_j = %d\n"
              "dv_dot_eta_i_Hubble = %g, dv_dot_eta_j_Hubble = %g\n"
              "vi_reconstructed = (%g, %g, %g), "
              "vj_reconstructed = (%g, %g, %g)\n"
              "dv_ij = (%g, %g, %g), dv_ji = (%g, %g, %g)\n"
              "eta_i = (%g, %g, %g), eta_j = (%g, %g, %g)\n"
              "eta_i2 = %g, eta_j2 = %g\n"
              "hi_inv = %g, hj_inv = %g\n"
              "a2_Hubble = %g, r2 = %g\n",
              pi->id, pj->id,
              conv_i, conv_j,
              dv_dot_eta_i_Hubble, dv_dot_eta_j_Hubble,
              vi_reconstructed[0], vi_reconstructed[1], vi_reconstructed[2],
              vj_reconstructed[0], vj_reconstructed[1], vj_reconstructed[2],
              dv_ij[0], dv_ij[1], dv_ij[2],
              dv_ji[0], dv_ji[1], dv_ji[2],
              eta_i[0], eta_i[1], eta_i[2],
              eta_j[0], eta_j[1], eta_j[2],
              eta_i2, eta_j2,
              hi_inv, hj_inv,
              a2_Hubble, r2);
    }
#endif

    /* mu_i and mu_j include the Hubble flow */
    const hydro_real_t mu_i = 
        conv * dv_dot_eta_i_Hubble  / (eta_i2 + const_viscosity_epsilon2);
    const hydro_real_t mu_j = 
        conv * dv_dot_eta_j_Hubble  / (eta_j2 + const_viscosity_epsilon2);
    
    const hydro_real_t Q_i_alpha = 
        -const_viscosity_alpha * pi->force.soundspeed * mu_i;
    const hydro_real_t Q_i_beta = const_viscosity_beta * mu_i * mu_i;
    const hydro_real_t Q_i = rhoi * (Q_i_alpha + Q_i_beta);

    const hydro_real_t Q_j_alpha =
        -const_viscosity_alpha * pj->force.soundspeed * mu_j;
    const hydro_real_t Q_j_beta = const_viscosity_beta * mu_j * mu_j;
    const hydro_real_t Q_j = rhoj * (Q_j_alpha + Q_j_beta);

    /* Add viscosity to the pressure */
    visc_acc_term = (Q_i + Q_j) * rhoij_inv;
#else
    const hydro_real_t conv = (dv_dot_dr_Hubble < 0.) ? fac_mu : 0.;
    const hydro_real_t mu_ij = conv * dv_dot_dr_Hubble * r_inv;
    const hydro_real_t c_ij = 
        0.5 * (pi->force.soundspeed + pj->force.soundspeed);

    const hydro_real_t Q_ij_alpha = -const_viscosity_alpha * c_ij * mu_ij;
    const hydro_real_t Q_ij_beta = const_viscosity_beta * mu_ij * mu_ij;
    const hydro_real_t Q_ij = (rhoi + rhoj) * (Q_ij_alpha + Q_ij_beta);
    
    /* Add viscosity to the pressure */
    visc_acc_term = Q_ij * rhoij_inv;
#endif

    /* Split heating between the two particles */
    visc_du_term = 0.5 * visc_acc_term;


    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;

    if (u_well_conditioned) {
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

    /* Signal velocity is the velocity difference projected with Hubble flow */
    const hydro_real_t v_sig_speed = fac_mu * fabs(dv_dot_dr_Hubble) * r_inv;
    const hydro_real_t v_sig_pressure = 
        sqrt(fabs(pressurei - pressurej) * rho_ij_inv);
    const hydro_real_t v_sig_cond = v_sig_speed + v_sig_pressure;

    const hydro_real_t du_ij = ui_reconstructed - uj_reconstructed;

    /* Add conductivity to the specific energy */
    cond_du_term = 
        const_conductivity_alpha * fabs(v_sig_cond) * du_ij * rho_ij_inv;


    /* Finalize everything with the correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5 * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5 * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5 * (G_i[2] + G_j[2]);

    hydro_real_t G_ij_norm = 0.;
    hydro_vec3_vec3_dot(G_ij, G_ij, &G_ij_norm);
    G_ij_norm = sqrt(G_ij_norm);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    hydro_real_t dv_dot_G_ij = 0.;
#ifdef hydro_props_use_second_order_velocities_in_divergence
    hydro_vec3_vec3_dot(dv_ij, G_ij, &dv_dot_G_ij);
#else
    hydro_vec3_vec3_dot(dv, G_ij, &dv_dot_G_ij);                      
#endif

    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij;
    cond_du_term *= -G_ij_norm; /* Eq. 24 Rosswog 2020 */

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
    pj->force.h_dt -= mi * dv_dot_G_ij;


    /* Timestepping */


    /* New signal velocity */
#ifdef hydro_props_use_asymmetric_viscosity_mu
    const hydro_real_t v_sig_visc_i =
        signal_velocity(dx, pi, pj, mu_i, const_viscosity_beta);
    const hydro_real_t v_sig_visc_j =
        signal_velocity(dx, pj, pi, mu_j, const_viscosity_beta);
    const hydro_real_t v_sig_visc = max(v_sig_visc_i, v_sig_visc_j);
#else
    const hydro_real_t v_sig_visc =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);
#endif
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
    hydro_real_t dvdr = 0.;
    hydro_vec3_vec3_dot(dv, G_ij, &dvdr);

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

    const hydro_real_t hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const hydro_real_t hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const hydro_real_t wi_dr = hi_inv_dim_plus_one * wi_dx;
    const hydro_real_t wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const hydro_real_t kernel_gradient = 0.5 * (wi_dr + wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;

    sph_du_term_i *= dvdr * kernel_gradient;
    sph_du_term_j *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dvdr * r_inv * wi_dr;
    pj->force.h_dt -= mi * dvdr * r_inv * wj_dr;
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
  /* First order internal energy gradients */
  const char u_well_conditioned = 
      (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned);

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

    /* Need this for viscosity and conductivity */
    hydro_real_t dv_dot_dr_Hubble = 0.;
    hydro_vec3_vec3_dot(dv_ij, dx_ij, &dv_dot_dr_Hubble);
    dv_dot_dr_Hubble += a2_Hubble * r2;

#ifdef hydro_props_use_asymmetric_viscosity_mu
    const hydro_real_t dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const hydro_real_t eta_i[3] = {dx_ij[0] * hi_inv,
                                   dx_ij[1] * hi_inv, 
                                   dx_ij[2] * hi_inv};
    const hydro_real_t eta_j[3] = {dx_ji[0] * hj_inv,
                                   dx_ji[1] * hj_inv,
                                   dx_ji[2] * hj_inv};

    hydro_real_t eta_i2 = 0.;
    hydro_vec3_vec3_dot(eta_i, eta_i, &eta_i2);

    hydro_real_t eta_j2 = 0.;
    hydro_vec3_vec3_dot(eta_j, eta_j, &eta_j2);

    hydro_real_t dv_dot_eta_i = 0.;
    hydro_vec3_vec3_dot(dv_ij, eta_i, &dv_dot_eta_i);
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const hydro_real_t dv_dot_eta_i_Hubble = 
        dv_dot_eta_i + a2_Hubble * r2 * hi_inv;

    hydro_real_t dv_dot_eta_j = 0.;
    hydro_vec3_vec3_dot(dv_ji, eta_j, &dv_dot_eta_j);
    /* Scale Hubble flow by hj_inv so it is overall scaled */
    const hydro_real_t dv_dot_eta_j_Hubble = 
        dv_dot_eta_j + a2_Hubble * r2 * hj_inv;

    /* Is the flow converging? */
    const char conv_i = (dv_dot_eta_i_Hubble < 0.) ? 1 : 0;
    const char conv_j = (dv_dot_eta_j_Hubble < 0.) ? 1 : 0;

    /* Is the flow converging? If yes, multiply by fac_mu. If no, zero. */
    const hydro_real_t conv = (conv_i && conv_j) ? fac_mu : 0.;

#ifdef MAGMA2_DEBUG_CHECKS
    if (conv_i != conv_j) {
      warning("Flow mismatch for particles pi=%lld and pj=%lld.\n"
              "conv_i = %d, conv_j = %d\n"
              "dv_dot_eta_i_Hubble = %g, dv_dot_eta_j_Hubble = %g\n"
              "vi_reconstructed = (%g, %g, %g), "
              "vj_reconstructed = (%g, %g, %g)\n"
              "dv_ij = (%g, %g, %g), dv_ji = (%g, %g, %g)\n"
              "eta_i = (%g, %g, %g), eta_j = (%g, %g, %g)\n"
              "eta_i2 = %g, eta_j2 = %g\n"
              "hi_inv = %g, hj_inv = %g\n"
              "a2_Hubble = %g, r2 = %g\n",
              pi->id, pj->id,
              conv_i, conv_j,
              dv_dot_eta_i_Hubble, dv_dot_eta_j_Hubble,
              vi_reconstructed[0], vi_reconstructed[1], vi_reconstructed[2],
              vj_reconstructed[0], vj_reconstructed[1], vj_reconstructed[2],
              dv_ij[0], dv_ij[1], dv_ij[2],
              dv_ji[0], dv_ji[1], dv_ji[2],
              eta_i[0], eta_i[1], eta_i[2],
              eta_j[0], eta_j[1], eta_j[2],
              eta_i2, eta_j2,
              hi_inv, hj_inv,
              a2_Hubble, r2);
    }
#endif

    /* mu_i and mu_j include the Hubble flow */
    const hydro_real_t mu_i = 
        conv * dv_dot_eta_i_Hubble  / (eta_i2 + const_viscosity_epsilon2);
    const hydro_real_t mu_j = 
        conv * dv_dot_eta_j_Hubble  / (eta_j2 + const_viscosity_epsilon2);
    
    const hydro_real_t Q_i_alpha = 
        -const_viscosity_alpha * pi->force.soundspeed * mu_i;
    const hydro_real_t Q_i_beta = const_viscosity_beta * mu_i * mu_i;
    const hydro_real_t Q_i = rhoi * (Q_i_alpha + Q_i_beta);

    const hydro_real_t Q_j_alpha =
        -const_viscosity_alpha * pj->force.soundspeed * mu_j;
    const hydro_real_t Q_j_beta = const_viscosity_beta * mu_j * mu_j;
    const hydro_real_t Q_j = rhoj * (Q_j_alpha + Q_j_beta);

    /* Add viscosity to the pressure */
    visc_acc_term = (Q_i + Q_j) * rhoij_inv;
#else
    const hydro_real_t conv = (dv_dot_dr_Hubble < 0.) ? fac_mu : 0.;
    const hydro_real_t mu_ij = conv * dv_dot_dr_Hubble * r_inv;
    const hydro_real_t c_ij = 
        0.5 * (pi->force.soundspeed + pj->force.soundspeed);

    const hydro_real_t Q_ij_alpha = -const_viscosity_alpha * c_ij * mu_ij;
    const hydro_real_t Q_ij_beta = const_viscosity_beta * mu_ij * mu_ij;
    const hydro_real_t Q_ij = (rhoi + rhoj) * (Q_ij_alpha + Q_ij_beta);
    
    /* Add viscosity to the pressure */
    visc_acc_term = Q_ij * rhoij_inv;
#endif

    visc_du_term = 0.5 * visc_acc_term;

  
    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;

    if (u_well_conditioned) {
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

    /* Signal velocity is the velocity difference projected with Hubble flow */
    const hydro_real_t v_sig_speed = fac_mu * fabs(dv_dot_dr_Hubble) * r_inv;
    const hydro_real_t v_sig_pressure = 
        sqrt(fabs(pressurei - pressurej) * rho_ij_inv);
    const hydro_real_t v_sig_cond = v_sig_speed + v_sig_pressure;

    const hydro_real_t du_ij = ui_reconstructed - uj_reconstructed;

    /* Add conductivity to the specific energy */
    cond_du_term = 
        const_conductivity_alpha * fabs(v_sig_cond) * du_ij * rho_ij_inv;


    /* Finalize the viscosity and conductivity with correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5 * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5 * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5 * (G_i[2] + G_j[2]);

    hydro_real_t G_ij_norm = 0.;
    hydro_vec3_vec3_dot(G_ij, G_ij, &G_ij_norm);
    G_ij_norm = sqrt(G_ij_norm);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    hydro_real_t dv_dot_G_ij = 0.;
#ifdef hydro_props_use_second_order_velocities_in_divergence
    hydro_vec3_vec3_dot(dv_ij, G_ij, &dv_dot_G_ij);
#else
    hydro_vec3_vec3_dot(dv, G_ij, &dv_dot_G_ij);                      
#endif

    sph_du_term_i *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij;
    cond_du_term *= -G_ij_norm; /* Eq. 24 Rosswog 2020 */

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;


    /* Timestepping */


    /* New signal velocity */
#ifdef hydro_props_use_asymmetric_viscosity_mu
    const hydro_real_t v_sig_visc_i =
        signal_velocity(dx, pi, pj, mu_i, const_viscosity_beta);
    const hydro_real_t v_sig_visc_j =
        signal_velocity(dx, pj, pi, mu_j, const_viscosity_beta);
    const hydro_real_t v_sig_visc = max(v_sig_visc_i, v_sig_visc_j);
#else
    const hydro_real_t v_sig_visc = 
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);
#endif
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
    hydro_real_t dvdr = 0.;
    hydro_vec3_vec3_dot(dv, G_ij, &dvdr);

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

    const hydro_real_t hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const hydro_real_t hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const hydro_real_t wi_dr = hi_inv_dim_plus_one * wi_dx;
    const hydro_real_t wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const hydro_real_t kernel_gradient = 0.5 * (wi_dr + wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dvdr * r_inv * wi_dr;
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
