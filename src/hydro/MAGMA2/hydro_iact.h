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

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Smooth pressure gradient */
  pi->gradients.pressure[0] += faci * pj->u * dx[0];
  pi->gradients.pressure[1] += faci * pj->u * dx[1];
  pi->gradients.pressure[2] += faci * pj->u * dx[2];

  pj->gradients.pressure[0] -= facj * pi->u * dx[0];
  pj->gradients.pressure[1] -= facj * pi->u * dx[1];
  pj->gradients.pressure[2] -= facj * pi->u * dx[2];

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] +=
          mj * dv[i] * dx[k] * wi_dx * r_inv;
      pj->gradients.velocity_tensor_aux[i][k] +=
          mi * dv[i] * dx[k] * wj_dx * r_inv;

      pi->gradients.velocity_tensor_aux_norm[i][k] +=
          mj * dx[i] * dx[k] * wi_dx * r_inv;
      pj->gradients.velocity_tensor_aux_norm[i][k] +=
          mi * dx[i] * dx[k] * wj_dx * r_inv;
    }
  }

  /* Number of neighbors */
  pi->num_ngb++;
  pj->num_ngb++;

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */

  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
  pj->weighted_wcount += mi * r2 * wj_dx * r_inv;
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

  /* Compute pressure gradient */
  pi->gradients.pressure[0] += faci * pj->u * dx[0];
  pi->gradients.pressure[1] += faci * pj->u * dx[1];
  pi->gradients.pressure[2] += faci * pj->u * dx[2];

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  /* Equations 19 & 20 in Rosswog 2020. Signs are all positive because
   * dv * dx always results in a positive sign. */
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++) {
      pi->gradients.velocity_tensor_aux[i][k] +=
          mj * dv[i] * dx[k] * wi_dx * r_inv;

      pi->gradients.velocity_tensor_aux_norm[i][k] +=
          mj * dx[i] * dx[k] * wi_dx * r_inv;
    }
  }

  pi->num_ngb++;

  /* Correction factors for kernel gradients, and norm for the velocity
   * gradient. */
  pi->weighted_wcount += mj * r2 * wi_dx * r_inv;
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

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);
  pj->viscosity.v_sig = max(pj->viscosity.v_sig, new_v_sig);

  /* Calculate Del^2 u for the thermal diffusion coefficient. */
  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx, wj, wj_dx;

  const float ui = r / hi;
  const float uj = r / hj;

  kernel_deval(ui, &wi, &wi_dx);
  kernel_deval(uj, &wj, &wj_dx);

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  const float shock_ratio_j =
      pi->viscosity.tensor_norm > 0.f
          ? pi->viscosity.shock_indicator / pi->viscosity.tensor_norm
          : 0.f;

  pi->viscosity.shock_limiter += mj * shock_ratio_i * wi;
  pj->viscosity.shock_limiter += mi * shock_ratio_j * wj;

  /* Correction factors for kernel gradients */
  const float rho_inv_i = 1.f / rhoi;
  const float rho_inv_j = 1.f / rhoj;

  pi->weighted_neighbour_wcount += mj * r2 * wi_dx * rho_inv_j * r_inv;
  pj->weighted_neighbour_wcount += mi * r2 * wj_dx * rho_inv_i * r_inv;

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.correction_matrix[k][i] += mj * rho_inv_j * dx[k] * dx[i] * wi;
      pj->gradients.correction_matrix[k][i] += mi * rho_inv_i * dx[k] * dx[i] * wj;

      /* Indices in Rosswog 2020 are i for dv and k for dx. In this loop,
       * they are swapped just because correction_matrix is computed with
       * the paper indices. */
      pi->gradients.velocity_tensor[k][i] += mj * rho_inv_j * dv[k] * dx[i] * wi;
      pj->gradients.velocity_tensor[k][i] += mi * rho_inv_i * dv[k] * dx[i] * wj;

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
        pi->gradients.velocity_hessian[k][i][j] +=
            mj * rho_inv_j * dv_grad_ki * dx[j] * wi;
        pj->gradients.velocity_hessian[k][i][j] +=
            mi * rho_inv_i * dv_grad_ki * dx[j] * wj;
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

  /* We need to construct the maximal signal velocity between our particle
   * and all of it's neighbours */

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Cosmology terms for the signal velocity */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */

  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float new_v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta);

  /* Update if we need to */
  pi->viscosity.v_sig = max(pi->viscosity.v_sig, new_v_sig);

  /* Need to get some kernel values F_ij = wi_dx */
  float wi, wi_dx;

  const float ui = r / hi;

  kernel_deval(ui, &wi, &wi_dx);

  /* Calculate the shock limiter component */
  const float shock_ratio_i =
      pj->viscosity.tensor_norm > 0.f
          ? pj->viscosity.shock_indicator / pj->viscosity.tensor_norm
          : 0.f;

  pi->viscosity.shock_limiter += mj * shock_ratio_i * wi;

  /* Correction factors for kernel gradients */

  const float rho_inv_j = 1.f / rhoj;

  pi->weighted_neighbour_wcount += mj * r2 * wi_dx * rho_inv_j * r_inv;

  const float dv[3] = {pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1],
                       pi->v[2] - pj->v[2]};

  for (int k = 0; k < 3; k++) {
    for (int i = 0; i < 3; i++) {
      /* dx is signed as (pi - pj), but it is symmetric so we add */
      pi->gradients.correction_matrix[k][i] += 
          mj * rho_inv_j * dx[k] * dx[i] * wi;

      /* Indices in Rosswog 2020 are i for dv and k for dx. In this loop,
       * they are swapped just because correction_matrix is computed with
       * the paper indices. */
      pi->gradients.velocity_tensor[k][i] += 
          mj * rho_inv_j * dv[k] * dx[i] * wi;

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
        pi->gradients.velocity_hessian[k][i][j] += 
            mj * rho_inv_j * dv_grad_ki * dx[j] * wi;
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

  /* The acceleration vector */
  float G_ij[3] = {0.f, 0.f, 0.f};

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = 0.f;
  float sph_acc_term = (pressurei + pressurej) * rhoij_inv;

  float sph_du_term_i = pressurei * rhoij_inv;
  float sph_du_term_j = pressurej * rhoij_inv;
  float visc_du_term_i = 0.f;
  float visc_du_term_j = 0.f;

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
    float vi_reconstruct[3] = {0.f, 0.f, 0.f};
    float vj_reconstruct[3] = {0.f, 0.f, 0.f};

    const float dx_ij[3] = {dx[0], dx[1], dx[2]};
    const float dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* TODO: make this constant */
    const float num_ngb_ij = 
        0.5f * ((float)pi->num_ngb + (float)pj->num_ngb);

    hydro_slope_limiter(dx_ji, xi, xj, num_ngb_ij, pi->v, pj->v,
                        pi->gradients.velocity_tensor,
                        pj->gradients.velocity_tensor,
                        pi->gradients.velocity_hessian,
                        vi_reconstruct);

    hydro_slope_limiter(dx_ij, xj, xi, num_ngb_ij, pj->v, pi->v,
                        pj->gradients.velocity_tensor,
                        pi->gradients.velocity_tensor,
                        pj->gradients.velocity_hessian,
                        vj_reconstruct);

    /* Artificial viscosity */
    const float dv_ij[3] = {vi_reconstruct[0] - vj_reconstruct[0],
                            vi_reconstruct[1] - vj_reconstruct[1],
                            vi_reconstruct[2] - vj_reconstruct[2]};
    const float dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const float eta_i[3] = {dx_ij[0] * hi_inv,
                            dx_ij[1] * hi_inv, 
                            dx_ij[2] * hi_inv};
    const float eta_j[3] = {dx_ji[0] * hj_inv,
                            dx_ji[1] * hj_inv,
                            dx_ji[2] * hj_inv};
    const float eta_i2 = eta_i[0] * eta_i[0] + 
                         eta_i[1] * eta_i[1] +
                         eta_i[2] * eta_i[2];
    const float eta_j2 = eta_j[0] * eta_j[0] +
                         eta_j[1] * eta_j[1] + 
                         eta_j[2] * eta_j[2];
    const float dv_dot_eta_i = dv_ij[0] * eta_i[0] + 
                               dv_ij[1] * eta_i[1] +
                               dv_ij[2] * eta_i[2];
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const float dv_dot_eta_i_phys = dv_dot_eta_i + a2_Hubble * r2 * hi_inv;
    const float dv_dot_eta_j = dv_ji[0] * eta_j[0] +
                               dv_ji[1] * eta_j[1] +
                               dv_ji[2] * eta_j[2];
    /* Scale Hubble flow by hj_inv so it is overall scaled */
    const float dv_dot_eta_j_phys = dv_dot_eta_j + a2_Hubble * r2 * hj_inv;

    /* mu_i and mu_j are physical, not comoving */
    const float mu_i = 
        fminf(0.f, dv_dot_eta_i_phys  / (eta_i2 + const_viscosity_epsilon2));
    const float mu_j = 
        fminf(0.f, dv_dot_eta_j_phys  / (eta_j2 + const_viscosity_epsilon2));
    
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
    visc_du_term_i = Q_i * rhoij_inv;
    visc_du_term_j = Q_j * rhoij_inv;

    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5f * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5f * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5f * (G_i[2] + G_j[2]);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                               (pi->v[1] - pj->v[1]) * G_ij[1] +
                               (pi->v[2] - pj->v[2]) * G_ij[2];

    sph_du_term_i *= dv_dot_G_ij;
    sph_du_term_j *= dv_dot_G_ij;
    visc_du_term_i *= dv_dot_G_ij + a2_Hubble * r2;
    visc_du_term_j *= dv_dot_G_ij + a2_Hubble * r2;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
    pj->force.h_dt -= mi * dv_dot_G_ij;
  }
  else {

    /* Dot product with the distance vector to align
     * with dW/dr (dW/dr = r_ij / |r_ij| dot grad W */
    G_ij[0] = dx[0];
    G_ij[1] = dx[1];
    G_ij[2] = dx[2];

    /* Compute dv dot dr. */
    const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                       (pi->v[1] - pj->v[1]) * dx[1] +
                       (pi->v[2] - pj->v[2]) * dx[2];

    /* Includes the hubble flow term; not used for du/dt */
    const float dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const float omega_ij = min(dvdr_Hubble, 0.f);
    const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

    /* Construct the full viscosity term */
    const float rho_ij = rhoi + rhoj;
    const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
    const float visc =
        omega_ij < 0.f
            ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                   mu_ij +
               const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5f * rho_ij)
            : 0.f;

    visc_acc_term = visc;
    visc_du_term_i = 0.5f * visc_acc_term;
    visc_du_term_j = visc_du_term_i;

    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient =
        0.5f * r_inv * (wi_dr * pi->force.f + wj_dr * pj->force.f);

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;

    sph_du_term_i *= dvdr * kernel_gradient;
    sph_du_term_j *= dvdr * kernel_gradient;
    visc_du_term_i *= dvdr_Hubble * kernel_gradient;
    visc_du_term_j = visc_du_term_i;

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
  const float du_dt_i = sph_du_term_i + visc_du_term_i;
  const float du_dt_j = sph_du_term_j + visc_du_term_j;

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

  float G_ij[3] = {0.f, 0.f, 0.f};

  /* These are set whether or not we fall back onto SPH gradients */
  float visc_acc_term = 0.f;
  float sph_acc_term = (pressurei + pressurej) * rhoij_inv;
  float sph_du_term_i = pressurei * rhoij_inv;
  float visc_du_term_i = 0.f;

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
    float vi_reconstruct[3] = {0.f, 0.f, 0.f};
    float vj_reconstruct[3] = {0.f, 0.f, 0.f};

    const float dx_ij[3] = {dx[0], dx[1], dx[2]};
    const float dx_ji[3] = {-dx[0], -dx[1], -dx[2]};

    /* TODO: make this constant */
    const float num_ngb_ij = 
        0.5f * ((float)pi->num_ngb + (float)pj->num_ngb);

    hydro_slope_limiter(dx_ji, xi, xj, num_ngb_ij, pi->v, pj->v,
                        pi->gradients.velocity_tensor,
                        pj->gradients.velocity_tensor,
                        pi->gradients.velocity_hessian,
                        vi_reconstruct);

    hydro_slope_limiter(dx_ij, xj, xi, num_ngb_ij, pj->v, pi->v,
                        pj->gradients.velocity_tensor,
                        pi->gradients.velocity_tensor,
                        pj->gradients.velocity_hessian,
                        vj_reconstruct);

    /* Artificial viscosity */
    const float dv_ij[3] = {vi_reconstruct[0] - vj_reconstruct[0],
                            vi_reconstruct[1] - vj_reconstruct[1],
                            vi_reconstruct[2] - vj_reconstruct[2]};
    const float dv_ji[3] = {-dv_ij[0], -dv_ij[1], -dv_ij[2]};

    const float eta_i[3] = {dx_ij[0] * hi_inv,
                            dx_ij[1] * hi_inv, 
                            dx_ij[2] * hi_inv};
    const float eta_j[3] = {dx_ji[0] * hj_inv,
                            dx_ji[1] * hj_inv,
                            dx_ji[2] * hj_inv};
    const float eta_i2 = eta_i[0] * eta_i[0] + 
                         eta_i[1] * eta_i[1] +
                         eta_i[2] * eta_i[2];
    const float eta_j2 = eta_j[0] * eta_j[0] +
                         eta_j[1] * eta_j[1] + 
                         eta_j[2] * eta_j[2];
    const float dv_dot_eta_i = dv_ij[0] * eta_i[0] + 
                               dv_ij[1] * eta_i[1] +
                               dv_ij[2] * eta_i[2];
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const float dv_dot_eta_i_phys = dv_dot_eta_i + a2_Hubble * r2 * hi_inv;
    const float dv_dot_eta_j = dv_ji[0] * eta_j[0] +
                               dv_ji[1] * eta_j[1] +
                               dv_ji[2] * eta_j[2];
    /* Scale Hubble flow by hi_inv so it is overall scaled */
    const float dv_dot_eta_j_phys = dv_dot_eta_j + a2_Hubble * r2 * hj_inv;

    /* mu_i and mu_j are physical, not comoving */
    const float mu_i = 
        fminf(0.f, dv_dot_eta_i_phys  / (eta_i2 + const_viscosity_epsilon2));
    const float mu_j = 
        fminf(0.f, dv_dot_eta_j_phys  / (eta_j2 + const_viscosity_epsilon2));
    
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
    visc_du_term_i = Q_i * rhoij_inv;

    /* Averaged correction gradient. Note: antisymmetric, so only need
     * a sign flip for pj */
    G_ij[0] = 0.5f * (G_i[0] + G_j[0]);
    G_ij[1] = 0.5f * (G_i[1] + G_j[1]);
    G_ij[2] = 0.5f * (G_i[2] + G_j[2]);

    /* Compute dv dot G_ij, reduces to dv dot dx in regular SPH. */
    const double dv_dot_G_ij = (pi->v[0] - pj->v[0]) * G_ij[0] +
                               (pi->v[1] - pj->v[1]) * G_ij[1] +
                               (pi->v[2] - pj->v[2]) * G_ij[2];

    sph_du_term_i *= dv_dot_G_ij;
    visc_du_term_i *= dv_dot_G_ij + a2_Hubble * r2;

    /* Get the time derivative for h. */
    pi->force.h_dt -= mj * dv_dot_G_ij;
  }
  else {
    
    /* Compute dv dot dr. */
    const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                       (pi->v[1] - pj->v[1]) * dx[1] +
                       (pi->v[2] - pj->v[2]) * dx[2];

    /* Includes the hubble flow term; not used for du/dt */
    const float dvdr_Hubble = dvdr + a2_Hubble * r2;

    /* Are the particles moving towards each others ? */
    const float omega_ij = min(dvdr_Hubble, 0.f);
    const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

    /* Construct the full viscosity term */
    const float rho_ij = rhoi + rhoj;
    const float alpha = pi->viscosity.alpha + pj->viscosity.alpha;
    const float visc =
        omega_ij < 0.f
            ? (-0.25f * alpha * (pi->force.soundspeed + pj->force.soundspeed) *
                  mu_ij +
              const_viscosity_beta * mu_ij * mu_ij) /
                  (0.5f * rho_ij)
            : 0.f;

    visc_acc_term = visc;
    visc_du_term_i = 0.5f * visc_acc_term;
    
    const float hi_inv_dim_plus_one = hi_inv * hi_inv_dim; /* 1/h^(d+1) */
    const float hj_inv_dim_plus_one = hj_inv * hj_inv_dim;
    const float wi_dr = hi_inv_dim_plus_one * wi_dx;
    const float wj_dr = hj_inv_dim_plus_one * wj_dx;

    /* Variable smoothing length term */
    const float kernel_gradient =
        0.5f * r_inv * (wi_dr * pi->force.f + wj_dr * pj->force.f);

    visc_acc_term *= kernel_gradient;
    sph_acc_term *= kernel_gradient;
    sph_du_term_i *= dvdr * kernel_gradient;
    visc_du_term_i *= dvdr_Hubble * kernel_gradient;

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
  const float du_dt_i = sph_du_term_i + visc_du_term_i;

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

}

#endif /* SWIFT_MAGMA2_HYDRO_IACT_H */
