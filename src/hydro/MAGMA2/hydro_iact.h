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

  /* For slope limiter */
  pi->gradients.kernel_size = fmax(r, pi->gradients.kernel_size);
  pj->gradients.kernel_size = fmax(r, pj->gradients.kernel_size);

  /* Now we need to compute the derivative terms */
  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;
  const hydro_real_t faci = mj * wi_dx * r_inv;
  const hydro_real_t facj = mi * wj_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */
  const hydro_real_t du = pi->u - pj->u;

  /* For slope limiter */
  pi->gradients.du_min = fmin(-du, pi->gradients.du_min);
  pi->gradients.du_max = fmax(-du, pi->gradients.du_max);

  pj->gradients.du_min = fmin(du, pj->gradients.du_min);
  pj->gradients.du_max = fmax(du, pj->gradients.du_max);

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

    /* For slope limiter */
    pi->gradients.dv_min[i] = fmin(-dv[i], pi->gradients.dv_min[i]);
    pi->gradients.dv_max[i] = fmax(-dv[i], pi->gradients.dv_max[i]);

    pj->gradients.dv_min[i] = fmin(dv[i], pj->gradients.dv_min[i]);
    pj->gradients.dv_max[i] = fmax(dv[i], pj->gradients.dv_max[i]);

    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * faci;
      pj->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * facj;

      pi->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * faci;
      pj->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * facj;
    }
  }

#ifdef MAGMA2_DEBUG_CHECKS
  /* Number of neighbors */
  pi->debug.num_ngb++;
  pj->debug.num_ngb++;
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

  /* For slope limiter */
  pi->gradients.kernel_size = fmax(r, pi->gradients.kernel_size);

  const hydro_real_t r_inv = r ? 1.0 / r : 0.0;
  const hydro_real_t faci = mj * wi_dx * r_inv;

  /* Equations 19 & 20 in Rosswog 2020. Compute the internal energy auxiliary 
   * vector and norm for the gradient */
  const hydro_real_t du = pi->u - pj->u;

  /* For slope limiter */
  pi->gradients.du_min = fmin(-du, pi->gradients.du_min);
  pi->gradients.du_max = fmax(-du, pi->gradients.du_max);

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

    /* For slope limiter */
    pi->gradients.dv_min[i] = fmin(-dv[i], pi->gradients.dv_min[i]);
    pi->gradients.dv_max[i] = fmax(-dv[i], pi->gradients.dv_max[i]);

    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] += dv[i] * dx[k] * faci;
      pi->gradients.velocity_tensor_aux_norm[i][k] += dx[i] * dx[k] * faci;
    }
  }

#ifdef MAGMA2_DEBUG_CHECKS
  /* Neighbour number */
  pi->debug.num_ngb++;
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
    const hydro_real_t du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[k];

    for (int i = 0; i < 3; i++) {
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
    const hydro_real_t du_k = pi->gradients.u_aux[k] - pj->gradients.u_aux[k];

    for (int i = 0; i < 3; i++) {
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
  const hydro_real_t a2_Hubble = a * a * H;

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
  const hydro_real_t dv_raw[3] = {dv[0], dv[1], dv[2]};

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  hydro_real_t G_ij[3] = {0., 0., 0.};
  hydro_real_t G_rad_ij[3] = {0., 0., 0.};
  hydro_real_t G_ij_norm = 0.;
  hydro_real_t G_rad_ij_norm = 0.;

  /* Separation vectors and swapped */
  const hydro_real_t dx_ij[3] = {dx[0], dx[1], dx[2]};
  const hydro_real_t dx_ji[3] = {-dx[0], -dx[1], -dx[2]};
  hydro_real_t dx_ij_hat[3] = {0., 0., 0.};
  hydro_vec3_unit(dx_ij, dx_ij_hat);

  /* These are set whether or not we fall back onto SPH gradients */
  hydro_real_t visc_acc_term = 0.;
  hydro_real_t sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  hydro_real_t sph_du_term_i = pressurei * rhoij_inv;
  hydro_real_t sph_du_term_j = pressurej * rhoij_inv;
  hydro_real_t visc_du_term = 0.;
  hydro_real_t cond_du_term = 0.;

  /* Correction factor must be well-conditioned. */
  const unsigned char C_well_conditioned = 
      (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned);
  /* First order density-independent velocity gradients */
  const unsigned char D_well_conditioned = 
      (pi->gradients.D_well_conditioned && pj->gradients.D_well_conditioned);
  const unsigned char u_well_conditioned = 
      (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned);

   /* Flag to revert to use high order gradients */
  unsigned char high_order_gradients_flag = 0;

  /* Always use high order gradients when both pi and pj are well-conditioned */
  if (C_well_conditioned && D_well_conditioned && u_well_conditioned) {

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

    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    hydro_get_average_kernel_gradient(pi, pj, G_i, G_j, G_ij);

    /* Compute from j perspective to ensure perfectly anti-symmetric */
    G_j[0] = 0.;
    G_j[1] = 0.;
    G_j[2] = 0.;
  
    G_i[0] = 0.;
    G_i[1] = 0.;
    G_i[2] = 0.;

    /* Swap G_i and G_j */
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx_ij, G_j);
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx_ij, G_i);

    G_j[0] *= wj * hj_inv_dim;
    G_j[1] *= wj * hj_inv_dim;
    G_j[2] *= wj * hj_inv_dim;

    G_i[0] *= wi * hi_inv_dim;
    G_i[1] *= wi * hi_inv_dim;
    G_i[2] *= wi * hi_inv_dim;

    hydro_real_t G_ji[3] = {0., 0., 0.};
    hydro_get_average_kernel_gradient(pj, pi, G_j, G_i, G_ji);

    /* Average the two estimators */
    G_ij[0] = 0.5 * (G_ij[0] - G_ji[0]);
    G_ij[1] = 0.5 * (G_ij[1] - G_ji[1]);
    G_ij[2] = 0.5 * (G_ij[2] - G_ji[2]);

    /* Check if G_ij is extremely misaligned with the radial direction */
    G_ij_norm = hydro_vec3_norm(G_ij);

    /* Get G_ij along the separation vector */
    hydro_real_t G_ij_dot_dx_ij_hat = hydro_vec3_vec3_dot(G_ij, dx_ij_hat);
    const hydro_real_t G_ij_dot_dx_ij_hat_abs = fabs(G_ij_dot_dx_ij_hat);

    /* Find the cos(theta) term between G and dx */
    hydro_real_t cosine_G_ij_dx_ij_hat = 
        G_ij_dot_dx_ij_hat_abs / (G_ij_norm + 1.e-10);

    /* Handle floating point errors */
    if (cosine_G_ij_dx_ij_hat > 1.) cosine_G_ij_dx_ij_hat = 1.;

    const unsigned char G_has_large_angle = 
        (cosine_G_ij_dx_ij_hat < const_viscosity_cosine_limit);
    const unsigned char G_in_wrong_direction = (G_ij_dot_dx_ij_hat > 0.);

    /* Good angle between separation and correct direction, good to go! */
    if (!G_has_large_angle && !G_in_wrong_direction) {

      /* Make sure we use the correct interaction */
      high_order_gradients_flag = 1;

#ifdef hydro_props_use_radial_artificial_terms
      /* G along the separation vector */
      G_rad_ij[0] = G_ij_dot_dx_ij_hat * dx_ij_hat[0];
      G_rad_ij[1] = G_ij_dot_dx_ij_hat * dx_ij_hat[1];
      G_rad_ij[2] = G_ij_dot_dx_ij_hat * dx_ij_hat[2];
#else
      G_rad_ij[0] = G_ij[0];
      G_rad_ij[1] = G_ij[1];
      G_rad_ij[2] = G_ij[2];
#endif
    }
    else {
      /* Revert back to standard separation vector */
      G_ij[0] = dx[0];
      G_ij[1] = dx[1];
      G_ij[2] = dx[2];
      G_ij_norm = hydro_vec3_norm(G_ij);

      /* Use the separation vector for SPH */
      G_rad_ij[0] = dx[0];
      G_rad_ij[1] = dx[1];
      G_rad_ij[2] = dx[2];
    }
  }

  /* MAGMA2-style gradients (Matrix-Inversion-2 SPH) */
  if (high_order_gradients_flag) {

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

    const hydro_real_t phi_ji_vec = 
        hydro_vector_van_leer_phi(pj->gradients.velocity_tensor, 
                                  pi->gradients.velocity_tensor,
                                  dx_ji, xj, xi);

    /* Make sure no floating point problems */
    hydro_real_t phi_vec_sym = fmin(phi_ij_vec, phi_ji_vec);
    phi_vec_sym = fmin(1., phi_vec_sym);
    phi_vec_sym = fmax(0., phi_vec_sym);

    /* Need these recast in case of switching precision */
    const hydro_real_t v_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
    const hydro_real_t v_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_vec_sym, dx_ji, v_i, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_vec_sym, dx_ij, v_j, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);

    const hydro_real_t dv_reconstructed[3] = {
        vi_reconstructed[0] - vj_reconstructed[0],
        vi_reconstructed[1] - vj_reconstructed[1],
        vi_reconstructed[2] - vj_reconstructed[2]
    };

    /* Get velocity difference, but limit reconstructed values */
    hydro_real_t dv_ij[3] = {0., 0., 0.};
    hydro_vec_minmod_limiter(dv_reconstructed, dv_raw,
                             pi->gradients.dv_min, pi->gradients.dv_max,
                             pj->gradients.dv_min, pj->gradients.dv_max,
                             dv_ij);

    /* Artificial viscosity */


    /* Get the acceleration term, depends on the weighting scheme */
    visc_acc_term = 
        hydro_get_visc_acc_term(pi, pj, dv_ij, dx_ij, fac_mu, a2_Hubble);

    /* Split heating between the two particles */
    visc_du_term = 0.5 * visc_acc_term;


    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;


    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const hydro_real_t phi_ij_scalar = 
        hydro_scalar_van_leer_phi(pi->gradients.u, 
                                  pj->gradients.u,
                                  dx_ij, xi, xj);
    const hydro_real_t phi_ji_scalar = 
        hydro_scalar_van_leer_phi(pj->gradients.u, 
                                  pi->gradients.u,
                                  dx_ji, xj, xi);

    /* Make sure no floating point problems */
    hydro_real_t phi_scalar_sym = fmin(phi_ij_scalar, phi_ji_scalar);
    phi_scalar_sym = fmin(1., phi_scalar_sym);
    phi_scalar_sym = fmax(0., phi_scalar_sym);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_scalar_second_order_reconstruction(phi_scalar_sym, dx_ji, 
                                             (hydro_real_t)pi->u, 
                                             pi->gradients.u,
                                             pi->gradients.u_hessian,
                                             &ui_reconstructed);

    hydro_scalar_second_order_reconstruction(phi_scalar_sym, dx_ij, 
                                             (hydro_real_t)pj->u, 
                                             pj->gradients.u,
                                             pj->gradients.u_hessian,
                                             &uj_reconstructed);

    const hydro_real_t rho_ij = 0.5 * (rhoi + rhoj);
    const hydro_real_t rho_ij_inv = 1. / rho_ij;
    const hydro_real_t dv_Hubble[3] = {
        dv[0] + a2_Hubble * dx_ij[0],
        dv[1] + a2_Hubble * dx_ij[1],
        dv[2] + a2_Hubble * dx_ij[2]
    };

    const hydro_real_t dv_Hubble_dot_dx_ij_hat =
        hydro_vec3_vec3_dot(dv_Hubble, dx_ij_hat);

    /* Limit art. cond. to only when information is communicable */
    const hydro_real_t c_ij = 
        0.5 * (pi->force.soundspeed + pj->force.soundspeed);
    const hydro_real_t v_sig_alpha = const_viscosity_alpha_prefactor * c_ij;

    /* Must connect the particles along the LOS */
    hydro_real_t mu_ij = fac_mu * dv_Hubble_dot_dx_ij_hat;
    const hydro_real_t v_sig_beta = const_viscosity_beta_prefactor * mu_ij;

    /* Skip conduction if expansion beats sound speed along LOS */
    if (v_sig_alpha > v_sig_beta) {

      const hydro_real_t dv_ij_Hubble[3] = {
          dv_ij[0] + a2_Hubble * dx_ij[0],
          dv_ij[1] + a2_Hubble * dx_ij[1],
          dv_ij[2] + a2_Hubble * dx_ij[2]
      };

      /* Signal velocity from speed contributions */
#ifdef hydro_props_use_radial_artificial_terms
      const hydro_real_t dv_ij_Hubble_dot_dx_ij_hat =
          hydro_vec3_vec3_dot(dv_ij_Hubble, dx_ij_hat);
      const hydro_real_t v_sig_speed = 
          fac_mu * fabs(dv_ij_Hubble_dot_dx_ij_hat);
#else
      const hydro_real_t v_sig_speed = fac_mu * hydro_vec3_norm(dv_ij_Hubble);
#endif

      /* Get spec. energy difference, but limit reconstructed values */
      const hydro_real_t du_raw = pi->u - pj->u;
      const hydro_real_t du_reconstructed = ui_reconstructed - uj_reconstructed;
      const hydro_real_t du_ij = 
          hydro_scalar_minmod_limiter(du_reconstructed, du_raw,
                                      pi->gradients.du_min, pi->gradients.du_max,
                                      pj->gradients.du_min, pj->gradients.du_max);

      const hydro_real_t alpha_cond = const_conductivity_alpha;
      const hydro_real_t delta_P = fabs(pressurei - pressurej);
      const hydro_real_t P_lim = delta_P / (pressurei + pressurej);

      /* Add conductivity to the specific energy */
      cond_du_term = alpha_cond * P_lim * v_sig_speed * du_ij * rho_ij_inv;
    }
    else {
      mu_ij = 0.;
    }


    /* Finalize everything with the correct normalizations. */


    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const hydro_real_t dv_dot_G_ij = hydro_vec3_vec3_dot(dv_raw, G_ij);

#ifdef hydro_props_use_radial_artificial_terms
    /* Compute Hubble flow along LOS */
    const hydro_real_t dv_Hubble_along_dx_ij[3] = {
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[0],
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[1],
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[2]
    };

    const hydro_real_t dv_Hubble_dot_G_rad_ij = 
        hydro_vec3_vec3_dot(dv_Hubble_along_dx_ij, G_rad_ij);
#else
    const hydro_real_t dv_Hubble_dot_G_rad_ij =
        hydro_vec3_vec3_dot(dv_Hubble, G_rad_ij);
#endif

    /* Evolve the heating terms using the velocity divergence */
    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    cond_du_term *= -G_rad_ij_norm; /* Eq. 24 Rosswog 2020 */
    visc_du_term *= dv_Hubble_dot_G_rad_ij;

    if (visc_du_term <= 0.) {
      visc_acc_term = 0.;
      visc_du_term = 0.;
    }


    /* Get the time derivative for h. */
    

    /* Velocity divergence is from the SPH estimator, not the G_ij vector */
    const hydro_real_t dv_dot_dx_ij = hydro_vec3_vec3_dot(dv_raw, dx_ij);
    pi->force.h_dt -= hydro_get_h_dt_sum(dv_dot_dx_ij, dv_dot_G_ij, 
                                         mj, rhoj_inv, r_inv, wi_dr);
    pj->force.h_dt -= hydro_get_h_dt_sum(dv_dot_dx_ij, dv_dot_G_ij, 
                                         mi, rhoi_inv, r_inv, wj_dr);


    /* Timestepping */


    /* Compute based on raw velocities */
    const hydro_real_t v_sig_visc =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = fmin(pi->h_min, h_ij);
    pj->h_min = fmin(pj->h_min, h_ij);

    /* New timestep estimate */
    const hydro_real_t dt_min_i = h_ij / v_sig_visc;
    const hydro_real_t dt_min_j = h_ij / v_sig_visc;
    pi->dt_min = fmin(pi->dt_min, dt_min_i);
    pj->dt_min = fmin(pj->dt_min, dt_min_j);

#ifdef MAGMA2_DEBUG_CHECKS
    pi->debug.N_force_high_order_grad++;
    pj->debug.N_force_high_order_grad++;
#endif
  }
  else { /* Gasoline-like SPH fallback */


    /* Compute dv dot dr. */
    const hydro_real_t dvdr = hydro_vec3_vec3_dot(dv_raw, G_ij);

    /* Includes the hubble flow term; not used for du/dt */
    const hydro_real_t dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const hydro_real_t omega_ij = fmin(dvdr_Hubble, 0.);
    
    hydro_real_t mu_full_ij = fac_mu * r_inv * dvdr_Hubble;
    const hydro_real_t mu_ij = fac_mu * r_inv * omega_ij;

    /* Construct the full viscosity term */
    const hydro_real_t b_ij = 
          0.5 * (pi->gradients.balsara + pj->gradients.balsara);
    const hydro_real_t rho_ij = rhoi + rhoj;
    const hydro_real_t cs_ij = pi->force.soundspeed + pj->force.soundspeed;
    const hydro_real_t c_ij = 0.5 * cs_ij;
    const hydro_real_t alpha = const_viscosity_alpha;
    const hydro_real_t visc =
        omega_ij < 0.
            ? (-0.25 * alpha * cs_ij *
                   mu_ij +
               const_viscosity_beta * mu_ij * mu_ij) * b_ij /
                  (0.5 * rho_ij)
            : 0.;

    visc_acc_term = const_fallback_reduction_factor * visc;
    visc_du_term = 0.5 * visc_acc_term;
    
    const hydro_real_t v_sig_alpha = const_viscosity_alpha_prefactor * c_ij;
    const hydro_real_t v_sig_beta = const_viscosity_beta_prefactor * mu_full_ij;

    if (v_sig_alpha > v_sig_beta) {
      const hydro_real_t rho_ij_inv = 2. / rho_ij;
      const hydro_real_t du = pi->u - pj->u;

      const hydro_real_t alpha_cond = const_conductivity_alpha;
      const hydro_real_t delta_P = fabs(pressurei - pressurej);
      const hydro_real_t P_lim = delta_P / (pressurei + pressurej);

      cond_du_term = alpha_cond * P_lim * fabs(mu_full_ij) * du * rho_ij_inv;
      cond_du_term *= const_fallback_reduction_factor;
    }

    /* New signal velocity */
    const hydro_real_t new_v_sig_visc =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = fmin(pi->h_min, h_ij);
    pj->h_min = fmin(pj->h_min, h_ij);
 
    /* New timestep estimate */
    const hydro_real_t dt_min_i = h_ij / new_v_sig_visc;
    const hydro_real_t dt_min_j = h_ij / new_v_sig_visc;
    pi->dt_min = fmin(pi->dt_min, dt_min_i);
    pj->dt_min = fmin(pj->dt_min, dt_min_j);

    const hydro_real_t kernel_gradient = 0.5 * (wi_dr + wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;

    sph_du_term_i *= dvdr * kernel_gradient;
    sph_du_term_j *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;
    cond_du_term *= kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= hydro_get_h_dt_sum(dvdr, 0., mj, rhoj_inv, r_inv, wi_dr);
    pj->force.h_dt -= hydro_get_h_dt_sum(dvdr, 0., mi, rhoi_inv, r_inv, wj_dr);

#ifdef MAGMA2_DEBUG_CHECKS
    pi->debug.N_force_low_order_grad++;
    pj->debug.N_force_low_order_grad++;
#endif
  }


  /* Get the time derivative for v. */


  /* Assemble the acceleration */
  const hydro_real_t acc[3] = {
      sph_acc_term * G_ij[0] + visc_acc_term * G_rad_ij[0],
      sph_acc_term * G_ij[1] + visc_acc_term * G_rad_ij[1],
      sph_acc_term * G_ij[2] + visc_acc_term * G_rad_ij[2]
  };

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc[0];
  pi->a_hydro[1] -= mj * acc[1];
  pi->a_hydro[2] -= mj * acc[2];

  pj->a_hydro[0] += mi * acc[0];
  pj->a_hydro[1] += mi * acc[1];
  pj->a_hydro[2] += mi * acc[2];


  /* Get the time derivative for u. */


  /* Assemble the energy equation term */
  const hydro_real_t du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;
  const hydro_real_t du_dt_j = sph_du_term_j + visc_du_term - cond_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  pi->u_dt_cond += cond_du_term * mj;
  pj->u_dt_cond -= cond_du_term * mi;
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
  const hydro_real_t a2_Hubble = a * a * H;

  const hydro_real_t r = sqrt(r2);
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
  const hydro_real_t hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
  const hydro_real_t wi_dr = hi_inv_dim_plus_one * wi_dx;
  const hydro_real_t wj_dr = hj_inv_dim_plus_one * wj_dx;

  /* Peculiar velocity difference vector */
  const hydro_real_t dv[3] = {pi->v[0] - pj->v[0],
                              pi->v[1] - pj->v[1],
                              pi->v[2] - pj->v[2]};
  const hydro_real_t dv_raw[3] = {dv[0], dv[1], dv[2]};

  /* For MAGMA2, this is the full anti-symmetric gradient vector. For the 
   * fallback Gasoline2-style SPH, this will just be the direction vector
   * between the two particles (r_i - r_j). */
  hydro_real_t G_ij[3] = {0., 0., 0.};
  hydro_real_t G_rad_ij[3] = {0., 0., 0.};
  hydro_real_t G_ij_norm = 0.;
  hydro_real_t G_rad_ij_norm = 0.;

  /* Separation vectors and swapped */
  const hydro_real_t dx_ij[3] = {dx[0], dx[1], dx[2]};
  const hydro_real_t dx_ji[3] = {-dx[0], -dx[1], -dx[2]};
  hydro_real_t dx_ij_hat[3] = {0., 0., 0.};
  hydro_vec3_unit(dx_ij, dx_ij_hat);

  /* These are set whether or not we fall back onto SPH gradients */
  hydro_real_t visc_acc_term = 0.;
  hydro_real_t sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  hydro_real_t sph_du_term_i = pressurei * rhoij_inv;
  hydro_real_t visc_du_term = 0.;
  hydro_real_t cond_du_term = 0.;

  /* Correction factor must be well-conditioned. */
  const unsigned char C_well_conditioned = 
      (pi->gradients.C_well_conditioned && pj->gradients.C_well_conditioned);
  /* First order density-independent velocity gradients */
  const unsigned char D_well_conditioned = 
      (pi->gradients.D_well_conditioned && pj->gradients.D_well_conditioned);
  const unsigned char u_well_conditioned = 
      (pi->gradients.u_well_conditioned && pj->gradients.u_well_conditioned);

  /* Flag to use second order gradients */
  unsigned char high_order_gradients_flag = 0;

  /* Always use SPH gradients between particles if one of them has an
   * ill-conditioned C matrix */
  if (C_well_conditioned && D_well_conditioned && u_well_conditioned) {

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

    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    hydro_get_average_kernel_gradient(pi, pj, G_i, G_j, G_ij);

        /* Compute from j perspective to ensure perfectly anti-symmetric */
    G_j[0] = 0.;
    G_j[1] = 0.;
    G_j[2] = 0.;
  
    G_i[0] = 0.;
    G_i[1] = 0.;
    G_i[2] = 0.;

    /* Swap G_i and G_j */
    hydro_mat3x3_vec3_dot(pj->gradients.correction_matrix, dx_ij, G_j);
    hydro_mat3x3_vec3_dot(pi->gradients.correction_matrix, dx_ij, G_i);

    G_j[0] *= wj * hj_inv_dim;
    G_j[1] *= wj * hj_inv_dim;
    G_j[2] *= wj * hj_inv_dim;

    G_i[0] *= wi * hi_inv_dim;
    G_i[1] *= wi * hi_inv_dim;
    G_i[2] *= wi * hi_inv_dim;

    hydro_real_t G_ji[3] = {0., 0., 0.};
    hydro_get_average_kernel_gradient(pj, pi, G_j, G_i, G_ji);

    /* Average the two estimators */
    G_ij[0] = 0.5 * (G_ij[0] - G_ji[0]);
    G_ij[1] = 0.5 * (G_ij[1] - G_ji[1]);
    G_ij[2] = 0.5 * (G_ij[2] - G_ji[2]);
  
    /* Check if G_ij is extremely misaligned with the radial direction */
    G_ij_norm = hydro_vec3_norm(G_ij);

    /* Get G_ij along the separation vector */
    hydro_real_t G_ij_dot_dx_ij_hat = hydro_vec3_vec3_dot(G_ij, dx_ij_hat);
    const hydro_real_t G_ij_dot_dx_ij_hat_abs = fabs(G_ij_dot_dx_ij_hat);

    /* Find the cos(theta) term between G and dx */
    hydro_real_t cosine_G_ij_dx_ij_hat = 
        G_ij_dot_dx_ij_hat_abs / (G_ij_norm + 1.e-10);

    /* Handle floating point errors */
    if (cosine_G_ij_dx_ij_hat > 1.) cosine_G_ij_dx_ij_hat = 1.;

    const unsigned char G_has_large_angle = 
        (cosine_G_ij_dx_ij_hat < const_viscosity_cosine_limit);
    const unsigned char G_in_wrong_direction = (G_ij_dot_dx_ij_hat > 0.);

    /* Good angle between separation and correct direction, good to go! */
    if (!G_has_large_angle && !G_in_wrong_direction) {

      /* Make sure we use the correct gradients below */
      high_order_gradients_flag = 1;

#ifdef hydro_props_use_radial_artificial_terms
      /* G along the separation vector */
      G_rad_ij[0] = G_ij_dot_dx_ij_hat * dx_ij_hat[0];
      G_rad_ij[1] = G_ij_dot_dx_ij_hat * dx_ij_hat[1];
      G_rad_ij[2] = G_ij_dot_dx_ij_hat * dx_ij_hat[2];
#else
      G_rad_ij[0] = G_ij[0];
      G_rad_ij[1] = G_ij[1];
      G_rad_ij[2] = G_ij[2];
#endif
    }
    else {
      /* Revert back to standard separation vector */
      G_ij[0] = dx[0];
      G_ij[1] = dx[1];
      G_ij[2] = dx[2];
      G_ij_norm = hydro_vec3_norm(G_ij);

      /* Use the separation vector for SPH */
      G_rad_ij[0] = dx[0];
      G_rad_ij[1] = dx[1];
      G_rad_ij[2] = dx[2];
    }

    G_rad_ij_norm = hydro_vec3_norm(G_rad_ij);
  }

  /* MAGMA2 style gradients (Matrix-Inversion-2 SPH) */
  if (high_order_gradients_flag) {

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
    const hydro_real_t phi_ji_vec = 
        hydro_vector_van_leer_phi(pj->gradients.velocity_tensor, 
                                  pi->gradients.velocity_tensor,
                                  dx_ji, xj, xi);

    /* Make sure no floating point problems */
    hydro_real_t phi_vec_sym = fmin(phi_ij_vec, phi_ji_vec);
    phi_vec_sym = fmin(1., phi_vec_sym);
    phi_vec_sym = fmax(0., phi_vec_sym);

    /* Need these recast in case of switching precision */
    const hydro_real_t v_i[3] = {pi->v[0], pi->v[1], pi->v[2]};
    const hydro_real_t v_j[3] = {pj->v[0], pj->v[1], pj->v[2]};

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_vector_second_order_reconstruction(phi_vec_sym, dx_ji, v_i, 
                                             pi->gradients.velocity_tensor,
                                             pi->gradients.velocity_hessian,
                                             vi_reconstructed);

    hydro_vector_second_order_reconstruction(phi_vec_sym, dx_ij, v_j, 
                                             pj->gradients.velocity_tensor,
                                             pj->gradients.velocity_hessian,
                                             vj_reconstructed);

    const hydro_real_t dv_reconstructed[3] = {
        vi_reconstructed[0] - vj_reconstructed[0],
        vi_reconstructed[1] - vj_reconstructed[1],
        vi_reconstructed[2] - vj_reconstructed[2]
    };

    /* Get velocity difference, but limit reconstructed values */
    hydro_real_t dv_ij[3] = {0., 0., 0.};
    hydro_vec_minmod_limiter(dv_reconstructed, dv_raw,
                             pi->gradients.dv_min, pi->gradients.dv_max,
                             pj->gradients.dv_min, pj->gradients.dv_max,
                             dv_ij);


    /* Artificial viscosity */


    /* Get the acceleration term, depends on the weighting scheme */
    visc_acc_term = 
        hydro_get_visc_acc_term(pi, pj, dv_ij, dx_ij, fac_mu, a2_Hubble);

    /* Split heating between the two particles */
    visc_du_term = 0.5 * visc_acc_term;


    /* Artificial conductivity */


    hydro_real_t ui_reconstructed = 0.;
    hydro_real_t uj_reconstructed = 0.;

    /* Compute global Van Leer limiter (scalar, not component-wise) */
    const hydro_real_t phi_ij_scalar = 
        hydro_scalar_van_leer_phi(pi->gradients.u, 
                                  pj->gradients.u,
                                  dx_ij, xi, xj);
    const hydro_real_t phi_ji_scalar = 
        hydro_scalar_van_leer_phi(pj->gradients.u, 
                                  pi->gradients.u,
                                  dx_ji, xj, xi);

    /* Make sure no floating point problems */
    hydro_real_t phi_scalar_sym = fmin(phi_ij_scalar, phi_ji_scalar);
    phi_scalar_sym = fmin(1., phi_scalar_sym);
    phi_scalar_sym = fmax(0., phi_scalar_sym);

    /* dx_ji for particle i and dx_ij for particle j */
    hydro_scalar_second_order_reconstruction(phi_scalar_sym, dx_ji, 
                                             (hydro_real_t)pi->u, 
                                             pi->gradients.u,
                                             pi->gradients.u_hessian,
                                             &ui_reconstructed);

    hydro_scalar_second_order_reconstruction(phi_scalar_sym, dx_ij, 
                                             (hydro_real_t)pj->u, 
                                             pj->gradients.u,
                                             pj->gradients.u_hessian,
                                             &uj_reconstructed);

    const hydro_real_t rho_ij = 0.5 * (rhoi + rhoj);
    const hydro_real_t rho_ij_inv = 1. / rho_ij;
    const hydro_real_t dv_Hubble[3] = {
        dv[0] + a2_Hubble * dx_ij[0],
        dv[1] + a2_Hubble * dx_ij[1],
        dv[2] + a2_Hubble * dx_ij[2]
    };

    const hydro_real_t dv_Hubble_dot_dx_ij_hat =
        hydro_vec3_vec3_dot(dv_Hubble, dx_ij_hat);

    /* Limit art. cond. to only when information is communicable */
    const hydro_real_t c_ij = 
        0.5 * (pi->force.soundspeed + pj->force.soundspeed);
    const hydro_real_t v_sig_alpha = const_viscosity_alpha_prefactor * c_ij;

    /* Must connect the particles along the LOS */
    hydro_real_t mu_ij = fac_mu * dv_Hubble_dot_dx_ij_hat;
    const hydro_real_t v_sig_beta = const_viscosity_beta_prefactor * mu_ij;

    /* Skip conduction if expansion beats sound speed along LOS */
    if (v_sig_alpha > v_sig_beta) {

      const hydro_real_t dv_ij_Hubble[3] = {
          dv_ij[0] + a2_Hubble * dx_ij[0],
          dv_ij[1] + a2_Hubble * dx_ij[1],
          dv_ij[2] + a2_Hubble * dx_ij[2]
      };

      /* Signal velocity from speed contributions */
#ifdef hydro_props_use_radial_artificial_terms
      const hydro_real_t dv_ij_Hubble_dot_dx_ij_hat =
          hydro_vec3_vec3_dot(dv_ij_Hubble, dx_ij_hat);
      const hydro_real_t v_sig_speed = 
          fac_mu * fabs(dv_ij_Hubble_dot_dx_ij_hat);
#else
      const hydro_real_t v_sig_speed = fac_mu * hydro_vec3_norm(dv_ij_Hubble);
#endif

      /* Get spec. energy difference, but limit reconstructed values */
      const hydro_real_t du_raw = pi->u - pj->u;
      const hydro_real_t du_reconstructed = ui_reconstructed - uj_reconstructed;
      const hydro_real_t du_ij = 
          hydro_scalar_minmod_limiter(du_reconstructed, du_raw,
                                      pi->gradients.du_min, pi->gradients.du_max,
                                      pj->gradients.du_min, pj->gradients.du_max);

      const hydro_real_t alpha_cond = const_conductivity_alpha;
      const hydro_real_t delta_P = fabs(pressurei - pressurej);
      const hydro_real_t P_lim = delta_P / (pressurei + pressurej);

      /* Add conductivity to the specific energy */
      cond_du_term = alpha_cond * P_lim * v_sig_speed * du_ij * rho_ij_inv;
    }
    else {
      mu_ij = 0.;
    }

    /* Finalize the viscosity and conductivity with correct normalizations. */


    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const hydro_real_t dv_dot_G_ij = hydro_vec3_vec3_dot(dv_raw, G_ij);

#ifdef hydro_props_use_radial_artificial_terms
    /* Get Hubble flow along dx */
    const hydro_real_t dv_Hubble_along_dx_ij[3] = {
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[0],
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[1],
        dv_Hubble_dot_dx_ij_hat * dx_ij_hat[2]
    };

    /* Get Hubble flow contribution */
    const hydro_real_t dv_Hubble_dot_G_rad_ij = 
        hydro_vec3_vec3_dot(dv_Hubble_along_dx_ij, G_rad_ij);
#else
    const hydro_real_t dv_Hubble_dot_G_rad_ij =
        hydro_vec3_vec3_dot(dv_Hubble, G_rad_ij);
#endif

    sph_du_term_i *= dv_dot_G_ij;
    cond_du_term *= -G_rad_ij_norm; /* Eq. 24 Rosswog 2020 */
    visc_du_term *= dv_Hubble_dot_G_rad_ij;

    if (visc_du_term <= 0.) {
      visc_acc_term = 0.;
      visc_du_term = 0.;
    }


    /* Get the time derivative for h. */


    const hydro_real_t dv_dot_dx_ij = hydro_vec3_vec3_dot(dv_raw, dx_ij);
    pi->force.h_dt -= hydro_get_h_dt_sum(dv_dot_dx_ij, dv_dot_G_ij, 
                                         mj, rhoj_inv, r_inv, wi_dr);


    /* Timestepping */


    /* Compute based on raw velocities */
    const hydro_real_t v_sig_visc =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = fmin(pi->h_min, h_ij);

    const hydro_real_t dt_min_i = h_ij / v_sig_visc;
    pi->dt_min = fmin(pi->dt_min, dt_min_i);

#ifdef MAGMA2_DEBUG_CHECKS
    pi->debug.N_force_high_order_grad++;
#endif
  }
  else { /* Gasoline2 SPH fallback*/


    /* Compute dv dot dr. */
    const hydro_real_t dvdr = hydro_vec3_vec3_dot(dv_raw, G_ij);

    /* Includes the hubble flow term; not used for du/dt */
    const hydro_real_t dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const hydro_real_t omega_ij = fmin(dvdr_Hubble, 0.);

    hydro_real_t mu_full_ij = fac_mu * r_inv * dvdr_Hubble;
    const hydro_real_t mu_ij = fac_mu * r_inv * omega_ij;

    /* Construct the full viscosity term */
    const hydro_real_t b_ij = 
        0.5 * (pi->gradients.balsara + pj->gradients.balsara);
    const hydro_real_t rho_ij = rhoi + rhoj;
    const hydro_real_t cs_ij = pi->force.soundspeed + pj->force.soundspeed;
    const hydro_real_t c_ij = 0.5 * cs_ij;
    const hydro_real_t alpha = const_viscosity_alpha;
    const hydro_real_t visc =
        omega_ij < 0.
            ? (-0.25 * alpha * cs_ij *
                  mu_ij +
              const_viscosity_beta * mu_ij * mu_ij) * b_ij /
                  (0.5 * rho_ij)
            : 0.;

    visc_acc_term = const_fallback_reduction_factor * visc;
    visc_du_term = 0.5 * visc_acc_term;
    
    const hydro_real_t v_sig_alpha = const_viscosity_alpha_prefactor * c_ij;
    const hydro_real_t v_sig_beta = const_viscosity_beta_prefactor * mu_full_ij;

    if (v_sig_alpha > v_sig_beta) {
      const hydro_real_t rho_ij_inv = 2. / rho_ij;
      const hydro_real_t du = pi->u - pj->u;

      const hydro_real_t alpha_cond = const_conductivity_alpha;
      const hydro_real_t delta_P = fabs(pressurei - pressurej);
      const hydro_real_t P_lim = delta_P / (pressurei + pressurej);

      cond_du_term = alpha_cond * P_lim * fabs(mu_full_ij) * du * rho_ij_inv;
      cond_du_term *= const_fallback_reduction_factor;
    }

    /* New signal velocity */
    const hydro_real_t new_v_sig_visc =
        signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

    const hydro_real_t h_ij = 0.5 * (hi + hj);
    pi->h_min = fmin(pi->h_min, h_ij);

    /* New time-step estimate */
    const hydro_real_t dt_min_i = h_ij / new_v_sig_visc;
    pi->dt_min = fmin(pi->dt_min, dt_min_i);

    /* Variable smoothing length term */
    const hydro_real_t kernel_gradient = 0.5 * (wi_dr + wj_dr) * r_inv;

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term *= dvdr_Hubble * kernel_gradient;
    cond_du_term *= kernel_gradient;

    /* Get the time derivative for h. */
    pi->force.h_dt -= hydro_get_h_dt_sum(dvdr, 0., mj, rhoj_inv, r_inv, wi_dr);

#ifdef MAGMA2_DEBUG_CHECKS
    pi->debug.N_force_low_order_grad++;
#endif
  }

  /* Assemble the acceleration */
  const hydro_real_t acc[3] = {
      sph_acc_term * G_ij[0] + visc_acc_term * G_rad_ij[0],
      sph_acc_term * G_ij[1] + visc_acc_term * G_rad_ij[1],
      sph_acc_term * G_ij[2] + visc_acc_term * G_rad_ij[2]
  };

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc[0];
  pi->a_hydro[1] -= mj * acc[1];
  pi->a_hydro[2] -= mj * acc[2];

  /* Assemble the energy equation term */
  const hydro_real_t du_dt_i = sph_du_term_i + visc_du_term + cond_du_term;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  pi->u_dt_cond += cond_du_term * mj;
}

#endif /* SWIFT_MAGMA2_HYDRO_IACT_H */
