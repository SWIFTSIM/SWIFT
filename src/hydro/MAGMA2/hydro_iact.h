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
 * @brief Density-Energy conservative implementation of SPH,
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
  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, hi_inv, mj);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  adaptive_softening_add_correction_term(pj, uj, hj_inv, mi);

  /* Now we need to compute the derivative terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */
  const float du = pi->u - pj->u;
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

  const float dv[3] = {pi->v[0] - pj->v[0],
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
  float wi, wi_dx;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */

  const float du = pi->u - pj->u;
  pi->gradients.u_aux[0] += du * dx[0] * faci;
  pi->gradients.u_aux[1] += du * dx[1] * faci;
  pi->gradients.u_aux[2] += du * dx[2] * faci;

  pi->gradients.u_aux_norm[0] += dx[0] * dx[0] * faci;
  pi->gradients.u_aux_norm[1] += dx[1] * dx[1] * faci;
  pi->gradients.u_aux_norm[2] += dx[2] * dx[2] * faci;

  const float dv[3] = {pi->v[0] - pj->v[0],
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
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  const float rhoi = hydro_get_comoving_density(pi);
  const float rhoj = hydro_get_comoving_density(pj);

  const float rhoi_inv = 1.f / rhoi;
  const float rhoj_inv = 1.f / rhoj;

  const float r = sqrtf(r2);

  float wi, wi_dx, wj, wj_dx;

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  const float faci = mj * rhoj_inv * wi;
  const float facj = mi * rhoi_inv * wj;

  /* Minimum smoothing length in the kernel */
  const float h_min = min(hi, hj);
  pi->h_min = min(pi->h_min, h_min);
  pj->h_min = min(pj->h_min, h_min);
  
  /* Compute all of the first-order gradients, and second-order gradients */

  /* Rosswog 2020 Equation 18 gradients. In the paper he uses (vj - vi) and
   * (rj - ri), however this is symmetric so no sign problems. */

  /* Internal energy gradient */
  const float du = pi->u - pj->u;

  pi->gradients.u[0] += du * dx[0] * faci;
  pi->gradients.u[1] += du * dx[1] * faci;
  pi->gradients.u[2] += du * dx[2] * faci;

  pj->gradients.u[0] += du * dx[0] * facj;
  pj->gradients.u[1] += du * dx[1] * facj;
  pj->gradients.u[2] += du * dx[2] * facj;

  /* Velocity gradients */
  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      const float du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[k];
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

      const float dv_grad_ki = pi->gradients.velocity_tensor_aux[k][i] -
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
  const float mj = hydro_get_mass(pj);
  const float rhoj = hydro_get_comoving_density(pj);
  const float rhoj_inv = 1.f / rhoj;
  const float r = sqrtf(r2);
  float wi, wi_dx;
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);
  const float faci = mj * rhoj_inv * wi;

  /* Minimum smoothing length in the kernel */
  const float h_min = min(hi, hj);
  pi->h_min = min(pi->h_min, h_min);

  /* Compute all of the first-order gradients, and second-order gradients */

  /* Rosswog 2020 Equation 18 gradients. In the paper he uses (vj - vi) and
   * (rj - ri), however this is symmetric so no sign problems. */

  /* Internal energy gradient */
  const float du = pi->u - pj->u;

  pi->gradients.u[0] += du * dx[0] * faci;
  pi->gradients.u[1] += du * dx[1] * faci;
  pi->gradients.u[2] += du * dx[2] * faci;

  /* Velocity gradients */
  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      const float du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[i];
      pi->gradients.u_hessian[k][i] += du_k * dx[i] * faci;

      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.correction_matrix[k][i] += dx[k] * dx[i] * faci;

      /* Indices in Rosswog 2020 are i for dv and k for dx. In this loop,
       * they are swapped just because correction_matrix is computed with
       * the paper indices. */
      pi->gradients.velocity_tensor[k][i] += dv[k] * dx[i] * faci;

      const float dv_grad_ki = pi->gradients.velocity_tensor_aux[k][i] - 
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
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float mi = pi->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float rhoi_inv = 1.f / rhoi;
  const float rhoj_inv = 1.f / rhoj;
  const float rhoij_inv = rhoi_inv * rhoj_inv;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  float G_ij[3] = {0.f, 0.f, 0.f};

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = 0.f;
  float sph_acc_term = (pressurei + pressurej) * rhoij_inv;

  float sph_du_term_i = pressurei * rhoij_inv;
  float sph_du_term_j = pressurej * rhoij_inv;
  float visc_du_term = 0.f;
  float cond_du_term = 0.f;

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  if (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned) {

    /* Corrected gradients */
    float G_i[3] = {0.f, 0.f, 0.f};
    float G_j[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx, G_i);
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx, G_j);

    /* Note: negative because dx is pj->x - pi->x in Rosswog 2020.
     * It is pj->x - pi->x for BOTH particles, and then sign flip
     * later */
    G_i[0] *= -wi * hi_inv_dim;
    G_i[1] *= -wi * hi_inv_dim;
    G_i[2] *= -wi * hi_inv_dim;

    G_j[0] *= -wj * hj_inv_dim;
    G_j[1] *= -wj * hj_inv_dim;
    G_j[2] *= -wj * hj_inv_dim;

    /* Velocity tensor dot product with separation vector for pi */
    float dvi_dx_dot_r[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pi->gradients.velocity_tensor, dx, dvi_dx_dot_r);

    /* Equation 17 has r_b - r_a, here we have pi - pj so the
     * sign is flipped (because of dx) */
    dvi_dx_dot_r[0] *= -0.5f;
    dvi_dx_dot_r[1] *= -0.5f;
    dvi_dx_dot_r[2] *= -0.5f;
    
    /* Velocity tensor dot product with separation vector for pi */
    float dvj_dx_dot_r[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pj->gradients.velocity_tensor, dx, dvj_dx_dot_r);

    /* Equation 17 has r_b - r_a, here we have pi - pj so the
     * sign is NOT flipped on dvj_dx_dot_r */
    dvj_dx_dot_r[0] *= 0.5f;
    dvj_dx_dot_r[1] *= 0.5f;
    dvj_dx_dot_r[2] *= 0.5f;

    /* Compute second order reconstruction of velocity between pi & pj */
    float vi_reconstructed[3] = {0.f, 0.f, 0.f};
    float vj_reconstructed[3] = {0.f, 0.f, 0.f};

    const float dx_ij[3] = {dx[0], dx[1], dx[2]};
    const float dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* TODO: make this constant */
    const float num_ngb_ij = 
        0.5f * ((float)pi->num_ngb + (float)pj->num_ngb);

    /* Rosswog 2020 h_i and h_j are actually kernel_gamma * h */
    const float hi_inv_kernel = hi_inv * kernel_gamma_inv;
    const float hj_inv_kernel = hj_inv * kernel_gamma_inv;
    const float eta_i_ij = r * hi_inv_kernel;
    const float eta_j_ij = r * hj_inv_kernel;

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const float phi_ij_vec = 
        hydro_vector_van_leer_phi(pi->gradients.velocity_tensor, 
                                  pj->gradients.velocity_tensor,
                                  dx_ij, eta_i_ij, eta_j_ij, num_ngb_ij);
    const float phi_ji_vec = 
        hydro_vector_van_leer_phi(pj->gradients.velocity_tensor,
                                  pi->gradients.velocity_tensor,
                                  dx_ji, eta_j_ij, eta_i_ij, num_ngb_ij);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ji, pi->v, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_ji_vec, dx_ij, pj->v, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);

    /* Artificial viscosity */
    const float dv_ij[3] = {vi_reconstructed[0] - vj_reconstructed[0],
                            vi_reconstructed[1] - vj_reconstructed[1],
                            vi_reconstructed[2] - vj_reconstructed[2]};
    const float dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const float eta_i[3] = {dx_ij[0] * hi_inv_kernel,
                            dx_ij[1] * hi_inv_kernel, 
                            dx_ij[2] * hi_inv_kernel};
    const float eta_j[3] = {dx_ji[0] * hj_inv_kernel,
                            dx_ji[1] * hj_inv_kernel,
                            dx_ji[2] * hj_inv_kernel};
    float eta_i2 = 0.f;
    hydro_vec3_vec3_dot(eta_i, eta_i, &eta_i2);

    float eta_j2 = 0.f;
    hydro_vec3_vec3_dot(eta_j, eta_j, &eta_j2);

    float dv_dot_eta_i = 0.f;
    hydro_vec3_vec3_dot(dv_ij, eta_i, &dv_dot_eta_i);
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const float dv_dot_eta_i_phys = 
        dv_dot_eta_i + a2_Hubble * r2 * hi_inv_kernel;

    float dv_dot_eta_j = 0.f;
    hydro_vec3_vec3_dot(dv_ji, eta_j, &dv_dot_eta_j);
    /* Scale Hubble flow by hj_inv so it is overall scaled */
    const float dv_dot_eta_j_phys = 
        dv_dot_eta_j + a2_Hubble * r2 * hj_inv_kernel;

    /* Is the flow converging? If so, multiply mu by fac_mu. If not, zero. */
    const float conv_i = (dv_dot_eta_i_phys < 0.f) ? fac_mu : 0.f;
    const float conv_j = (dv_dot_eta_j_phys < 0.f) ? fac_mu : 0.f;

    /* Both must be approaching */
    const float conv = (conv_i != 0.f && conv_j != 0.f) ? fac_mu : 0.f;

#ifdef MAGMA2_DEBUG_CHECKS
    if (conv_i != conv_j) {
      warning("Convergence factor mismatch for particles pi=%lld and pj=%lld.\n"
              "conv_i = %g, conv_j = %g\n"
              "dv_dot_eta_i_phys = %g, dv_dot_eta_j_phys = %g\n"
              "vi_reconstructed = (%g, %g, %g), vj_reconstructed = (%g, %g, %g)\n"
              "dv_ij = (%g, %g, %g), dv_ji = (%g, %g, %g)\n"
              "eta_i = (%g, %g, %g), eta_j = (%g, %g, %g)\n"
              "eta_i2 = %g, eta_j2 = %g\n"
              "hi_inv_kernel = %g, hj_inv_kernel = %g\n"
              "a2_Hubble = %g, r2 = %g\n",
              pi->id, pj->id,
              conv_i, conv_j,
              dv_dot_eta_i_phys, dv_dot_eta_j_phys,
              vi_reconstructed[0], vi_reconstructed[1], vi_reconstructed[2],
              vj_reconstructed[0], vj_reconstructed[1], vj_reconstructed[2],
              dv_ij[0], dv_ij[1], dv_ij[2],
              dv_ji[0], dv_ji[1], dv_ji[2],
              eta_i[0], eta_i[1], eta_i[2],
              eta_j[0], eta_j[1], eta_j[2],
              eta_i2, eta_j2,
              hi_inv_kernel, hj_inv_kernel,
              a2_Hubble, r2);
    }
#endif

    /* mu_i and mu_j include the Hubble flow */
    const float mu_i = 
        conv * dv_dot_eta_i_phys  / (eta_i2 + const_viscosity_epsilon2);
    const float mu_j = 
        conv * dv_dot_eta_j_phys  / (eta_j2 + const_viscosity_epsilon2);
    
    const float Q_i_alpha = 
        -const_viscosity_alpha * pi->force.soundspeed * mu_i;
    const float Q_i_beta = const_viscosity_beta * mu_i * mu_i;
    const float Q_i = rhoi * (Q_i_alpha + Q_i_beta);

    const float Q_j_alpha =
        -const_viscosity_alpha * pj->force.soundspeed * mu_j;
    const float Q_j_beta = const_viscosity_beta * mu_j * mu_j;
    const float Q_j = rhoj * (Q_j_alpha + Q_j_beta);

    /* Add viscosity to the pressure */
    visc_acc_term = (Q_i + Q_j) * rhoij_inv;
    visc_du_term = 0.5f * visc_acc_term;


    /* Artificial conductivity */


    float ui_reconstructed = 0.f;
    float uj_reconstructed = 0.f;

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const float phi_ij_scalar = 
        hydro_scalar_van_leer_phi(pi->gradients.u, 
                                  pj->gradients.u,
                                  dx_ij, eta_i_ij, eta_j_ij, num_ngb_ij);
    const float phi_ji_scalar = 
        hydro_scalar_van_leer_phi(pj->gradients.u,
                                  pi->gradients.u,
                                  dx_ji, eta_j_ij, eta_i_ij, num_ngb_ij);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ji, pi->u, 
                                             pi->gradients.u,
                                             pi->gradients.u_hessian,
                                             &ui_reconstructed);

    hydro_scalar_second_order_reconstruction(phi_ji_scalar, dx_ij, pj->u, 
                                             pj->gradients.u,
                                             pj->gradients.u_hessian,
                                             &uj_reconstructed);
  
    float dv_ij_dot_dv_ij = 0.f;
    hydro_vec3_vec3_dot(dv_ij, dv_ij, &dv_ij_dot_dv_ij);
    /* Signal velocity is the norm of the velocity difference vector */
    const float v_sig_cond = sqrtf(dv_ij_dot_dv_ij);

    const float du_ij = ui_reconstructed - uj_reconstructed;
    const float rho_ij = 0.5f * (rhoi + rhoj);

    /* Add conductivity to the specific energy */
    cond_du_term = -const_conductivity_alpha * v_sig_cond * du_ij / rho_ij;


    /* Finalize everything with the correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5f * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5f * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5f * (G_i[2] + G_j[2]);

    float G_ij_dot_G_ij = 0.f;
    hydro_vec3_vec3_dot(G_ij, G_ij, &G_ij_dot_G_ij);
    const float G_ij_norm = sqrtf(G_ij_dot_G_ij);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    /* TODO: Use second-order reconstructions here? */
    const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                               (pi->v[1] - pj->v[1]) * G_ij[1] +
                               (pi->v[2] - pj->v[2]) * G_ij[2];

    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij + a2_Hubble * r2;
    cond_du_term *= G_ij_norm;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
    pj->force.h_dt -= mi * dv_dot_G_ij;

    /* New signal velocity */
    const float v_sig_visc_i =
        signal_velocity(dx, pi, pj, mu_i, const_viscosity_beta);
    const float v_sig_visc_j =
        signal_velocity(dx, pj, pi, mu_j, const_viscosity_beta);
    const float v_sig_visc = max(v_sig_visc_i, v_sig_visc_j);
    const float v_sig = max(v_sig_visc, v_sig_cond);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, v_sig);
    pj->v_sig_max = max(pj->v_sig_max, v_sig);
  }
  else {

    /* Dot product with the distance vector to align
     * with dW/dr (dW/dr = r_ij / |r_ij| dot grad W */
    G_ij[0] = dx[0];
    G_ij[1] = dx[1];
    G_ij[2] = dx[2];

    /* Compute dv dot dr. */
    const float dv[3] = {pi->v[0] - pj->v[0],
                         pi->v[1] - pj->v[1],
                         pi->v[2] - pj->v[2]};
    float dvdr = 0.f;
    hydro_vec3_vec3_dot(dv, G_ij, &dvdr);

    /* Includes the hubble flow term; not used for du/dt */
    const float dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const float omega_ij = min(dvdr_Hubble, 0.f);
    const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

    /* Construct the full viscosity term */
    const float rho_ij = rhoi + rhoj;
    const float alpha = const_viscosity_alpha;
    const float visc =
        omega_ij < 0.f
            ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                   mu_ij +
               const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5f * rho_ij)
            : 0.f;

    /* New signal velocity */
    const float new_v_sig =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, new_v_sig);
    pj->v_sig_max = max(pj->v_sig_max, new_v_sig);

    visc_acc_term = visc;
    visc_du_term = 0.5f * visc_acc_term;

    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient = 0.5f * (wi_dr + wj_dr) * r_inv;

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
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  pj->a_hydro[0] += mi * acc * G_ij[0];
  pj->a_hydro[1] += mi * acc * G_ij[1];
  pj->a_hydro[2] += mi * acc * G_ij[2];

  /* Get the time derivative for u. */

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term - cond_du_term;

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
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;

  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float rhoi_inv = 1.f / rhoi;
  const float rhoj_inv = 1.f / rhoj;
  const float rhoij_inv = rhoi_inv * rhoj_inv;

  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj_inv_dim = pow_dimension(hj_inv);
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  float G_ij[3] = {0.f, 0.f, 0.f};

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = 0.f;
  float sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  float sph_du_term_i = pressurei * rhoij_inv;
  float visc_du_term = 0.f;
  float cond_du_term = 0.f;

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  if (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned) {

    /* Corrected gradients */
    float G_i[3] = {0.f, 0.f, 0.f};
    float G_j[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx, G_i);
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx, G_j);

    /* Note: negative because dx is pj->x - pi->x in Rosswog 2020.
     * It is pj->x - pi->x for BOTH particles, and then sign flip
     * later */
    G_i[0] *= -wi * hi_inv_dim;
    G_i[1] *= -wi * hi_inv_dim;
    G_i[2] *= -wi * hi_inv_dim;

    G_j[0] *= -wj * hj_inv_dim;
    G_j[1] *= -wj * hj_inv_dim;
    G_j[2] *= -wj * hj_inv_dim;

    /* Velocity tensor dot product with separation vector for pi */
    float dvi_dx_dot_r[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pi->gradients.velocity_tensor, dx, dvi_dx_dot_r);

    /* Equation 17 has r_b - r_a, here we have pi - pj so the
     * sign is flipped (because of dx) */
    dvi_dx_dot_r[0] *= -0.5f;
    dvi_dx_dot_r[1] *= -0.5f;
    dvi_dx_dot_r[2] *= -0.5f;
    
    /* Velocity tensor dot product with separation vector for pi */
    float dvj_dx_dot_r[3] = {0.f, 0.f, 0.f};
    hydro_mat3x3_vec3_dot(pj->gradients.velocity_tensor, dx, dvj_dx_dot_r);

    /* Equation 17 has r_b - r_a, here we have pi - pj so the
     * sign is NOT flipped on dvj_dx_dot_r */
    dvj_dx_dot_r[0] *= 0.5f;
    dvj_dx_dot_r[1] *= 0.5f;
    dvj_dx_dot_r[2] *= 0.5f;

    /* Compute second order reconstruction of velocity between pi & pj */
    float vi_reconstructed[3] = {0.f, 0.f, 0.f};
    float vj_reconstructed[3] = {0.f, 0.f, 0.f};

    const float dx_ij[3] = {dx[0], dx[1], dx[2]};
    const float dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* TODO: make this constant */
    const float num_ngb_ij = 
        0.5f * ((float)pi->num_ngb + (float)pj->num_ngb);

    /* Rosswog 2020 h_i and h_j are actually kernel_gamma * h */
    const float hi_inv_kernel = hi_inv * kernel_gamma_inv;
    const float hj_inv_kernel = hj_inv * kernel_gamma_inv;
    const float eta_i_ij = r * hi_inv_kernel;
    const float eta_j_ij = r * hj_inv_kernel;

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const float phi_ij_vec = 
        hydro_vector_van_leer_phi(pi->gradients.velocity_tensor, 
                                  pj->gradients.velocity_tensor,
                                  dx_ij, eta_i_ij, eta_j_ij, num_ngb_ij);
    const float phi_ji_vec = 
        hydro_vector_van_leer_phi(pj->gradients.velocity_tensor,
                                  pi->gradients.velocity_tensor,
                                  dx_ji, eta_j_ij, eta_i_ij, num_ngb_ij);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_ij_vec, dx_ji, pi->v, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_ji_vec, dx_ij, pj->v, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);

    /* Artificial viscosity */
    const float dv_ij[3] = {vi_reconstructed[0] - vj_reconstructed[0],
                            vi_reconstructed[1] - vj_reconstructed[1],
                            vi_reconstructed[2] - vj_reconstructed[2]};
    const float dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const float eta_i[3] = {dx_ij[0] * hi_inv_kernel,
                            dx_ij[1] * hi_inv_kernel, 
                            dx_ij[2] * hi_inv_kernel};
    const float eta_j[3] = {dx_ji[0] * hj_inv_kernel,
                            dx_ji[1] * hj_inv_kernel,
                            dx_ji[2] * hj_inv_kernel};

    float eta_i2 = 0.f;
    hydro_vec3_vec3_dot(eta_i, eta_i, &eta_i2);

    float eta_j2 = 0.f;
    hydro_vec3_vec3_dot(eta_j, eta_j, &eta_j2);

    float dv_dot_eta_i = 0.f;
    hydro_vec3_vec3_dot(dv_ij, eta_i, &dv_dot_eta_i);
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const float dv_dot_eta_i_phys = 
        dv_dot_eta_i + a2_Hubble * r2 * hi_inv_kernel;

    float dv_dot_eta_j = 0.f;
    hydro_vec3_vec3_dot(dv_ji, eta_j, &dv_dot_eta_j);
    /* Scale Hubble flow by hj_inv so it is overall scaled */
    const float dv_dot_eta_j_phys = 
        dv_dot_eta_j + a2_Hubble * r2 * hj_inv_kernel;

    /* Is the flow converging? */
    const float conv_i = (dv_dot_eta_i_phys < 0.f) ? fac_mu : 0.f;
    const float conv_j = (dv_dot_eta_j_phys < 0.f) ? fac_mu : 0.f;

    /* They both must be approaching. */
    const float conv = (conv_i != 0.f && conv_j != 0.f) ? fac_mu : 0.f;

#ifdef MAGMA2_DEBUG_CHECKS
    if (conv_i != conv_j) {
      warning("Convergence factor mismatch for particles pi=%lld and pj=%lld.\n"
              "conv_i = %g, conv_j = %g\n"
              "dv_dot_eta_i_phys = %g, dv_dot_eta_j_phys = %g\n"
              "vi_reconstructed = (%g, %g, %g), vj_reconstructed = (%g, %g, %g)\n"
              "dv_ij = (%g, %g, %g), dv_ji = (%g, %g, %g)\n"
              "eta_i = (%g, %g, %g), eta_j = (%g, %g, %g)\n"
              "eta_i2 = %g, eta_j2 = %g\n"
              "hi_inv_kernel = %g, hj_inv_kernel = %g\n"
              "a2_Hubble = %g, r2 = %g\n",
              pi->id, pj->id,
              conv_i, conv_j,
              dv_dot_eta_i_phys, dv_dot_eta_j_phys,
              vi_reconstructed[0], vi_reconstructed[1], vi_reconstructed[2],
              vj_reconstructed[0], vj_reconstructed[1], vj_reconstructed[2],
              dv_ij[0], dv_ij[1], dv_ij[2],
              dv_ji[0], dv_ji[1], dv_ji[2],
              eta_i[0], eta_i[1], eta_i[2],
              eta_j[0], eta_j[1], eta_j[2],
              eta_i2, eta_j2,
              hi_inv_kernel, hj_inv_kernel,
              a2_Hubble, r2);
    }
#endif

    /* mu_i and mu_j include the Hubble flow */
    const float mu_i = 
        conv * dv_dot_eta_i_phys  / (eta_i2 + const_viscosity_epsilon2);
    const float mu_j = 
        conv * dv_dot_eta_j_phys  / (eta_j2 + const_viscosity_epsilon2);
    
    const float Q_i_alpha = 
        -const_viscosity_alpha * pi->force.soundspeed * mu_i;
    const float Q_i_beta = const_viscosity_beta * mu_i * mu_i;
    const float Q_i = rhoi * (Q_i_alpha + Q_i_beta);

    const float Q_j_alpha =
        -const_viscosity_alpha * pj->force.soundspeed * mu_j;
    const float Q_j_beta = const_viscosity_beta * mu_j * mu_j;
    const float Q_j = rhoj * (Q_j_alpha + Q_j_beta);

    /* Add viscosity to the pressure */
    visc_acc_term = (Q_i + Q_j) * rhoij_inv;
    visc_du_term = 0.5f * visc_acc_term;

  
    /* Artificial conductivity */


    float ui_reconstructed = 0.f;
    float uj_reconstructed = 0.f;

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const float phi_ij_scalar = 
        hydro_scalar_van_leer_phi(pi->gradients.u, 
                                  pj->gradients.u,
                                  dx_ij, eta_i_ij, eta_j_ij, num_ngb_ij);
    const float phi_ji_scalar = 
        hydro_scalar_van_leer_phi(pj->gradients.u,
                                  pi->gradients.u,
                                  dx_ji, eta_j_ij, eta_i_ij, num_ngb_ij);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_scalar_second_order_reconstruction(phi_ij_scalar, dx_ji, pi->u, 
                                             pi->gradients.u,
                                             pi->gradients.u_hessian,
                                             &ui_reconstructed);

    hydro_scalar_second_order_reconstruction(phi_ji_scalar, dx_ij, pj->u, 
                                             pj->gradients.u,
                                             pj->gradients.u_hessian,
                                             &uj_reconstructed);
  
    float dv_ij_dot_dv_ij = 0.f;
    hydro_vec3_vec3_dot(dv_ij, dv_ij, &dv_ij_dot_dv_ij);
    const float v_sig_cond = sqrtf(dv_ij_dot_dv_ij);

    const float du_ij = ui_reconstructed - uj_reconstructed;
    const float rho_ij = 0.5f * (rhoi + rhoj);

    /* Add conductivity to the specific energy */
    cond_du_term = -const_conductivity_alpha * v_sig_cond * du_ij / rho_ij;


    /* Finalize the viscosity and conductivity with correct normalizations. */


    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5f * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5f * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5f * (G_i[2] + G_j[2]);

    float G_ij_dot_G_ij = 0.f;
    hydro_vec3_vec3_dot(G_ij, G_ij, &G_ij_dot_G_ij);
    const float G_ij_norm = sqrtf(G_ij_dot_G_ij);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    /* TODO: Use reconstructed velocities here? */
    const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                               (pi->v[1] - pj->v[1]) * G_ij[1] +
                               (pi->v[2] - pj->v[2]) * G_ij[2];

    sph_du_term_i *= dv_dot_G_ij;
    visc_du_term *= dv_dot_G_ij + a2_Hubble * r2;
    cond_du_term *= G_ij_norm; /* Eq. 24 Rosswog 2020 */

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;

    /* New signal velocity */
    const float v_sig_visc_i =
        signal_velocity(dx, pi, pj, mu_i, const_viscosity_beta);
    const float v_sig_visc_j =
        signal_velocity(dx, pj, pi, mu_j, const_viscosity_beta);
    const float v_sig_visc = max(v_sig_visc_i, v_sig_visc_j);
    const float v_sig = max(v_sig_visc, v_sig_cond);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, v_sig);
  }
  else {
    
    /* Fallback SPH points in the direction between the particles. */
    G_ij[0] = dx[0];
    G_ij[1] = dx[1];
    G_ij[2] = dx[2];

    /* Compute dv dot dr. */
    const float dv[3] = {pi->v[0] - pj->v[0],
                         pi->v[1] - pj->v[1],
                         pi->v[2] - pj->v[2]};
    float dvdr = 0.f;
    hydro_vec3_vec3_dot(dv, G_ij, &dvdr);

    /* Includes the hubble flow term; not used for du/dt */
    const float dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const float omega_ij = min(dvdr_Hubble, 0.f);
    const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

    /* Construct the full viscosity term */
    const float rho_ij = rhoi + rhoj;
    const float alpha = const_viscosity_alpha;
    const float visc =
        omega_ij < 0.f
            ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                  mu_ij +
              const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5f * rho_ij)
            : 0.f;

    /* New signal velocity */
    const float new_v_sig =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    /* Update if we need to */
    pi->v_sig_max = max(pi->v_sig_max, new_v_sig);

    visc_acc_term = visc;
    visc_du_term = 0.5f * visc_acc_term;

    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient = 0.5f * (wi_dr + wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dvdr * r_inv * wi_dr;
  }

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * G_ij[0];
  pi->a_hydro[1] -= mj * acc * G_ij[1];
  pi->a_hydro[2] -= mj * acc * G_ij[2];

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

}

#endif /* SWIFT_MAGMA2_HYDRO_IACT_H */
