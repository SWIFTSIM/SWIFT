/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Federico Stasyszyn (fstasyszyn@unc.edu.ar)
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
#ifndef SWIFT_VECTOR_POTENTIAL_MHD_IACT_H
#define SWIFT_VECTOR_POTENTIAL_MHD_IACT_H
#define VP_ADV_GAUGE

#include "periodic.h"

/**
 * @brief MHD-Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  double dA[3];
  for (int i = 0; i < 3; ++i)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];

  const double dAdr = dA[0] * dx[0] + dA[1] * dx[1] + dA[2] * dx[2];
  pi->mhd_data.divA -= faci * dAdr;
  pj->mhd_data.divA -= facj * dAdr;
  /////
  // bi = dj ak - dk aj
  // bj = dk ai - di ak
  // bk = di aj - dj ai
  //
  for (int i = 0; i < 3; ++i) {
    pi->mhd_data.BPred[i] += faci * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);
    pj->mhd_data.BPred[i] += facj * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);
  }
}

/**
 * @brief MHD-Density interaction between two particles. (non-symmetric)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mhd_density(const float r2, const float dx[3],
                               const float hi, const float hj,
                               struct part *restrict pi,
                               const struct part *restrict pj,
                               const double mu_0, const float a,
                               const float H) {
  float wi, wi_dx;

  const float r = sqrtf(r2);

  /* Get the mass. */
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;

  kernel_deval(ui, &wi, &wi_dx);

  /* Now we need to compute the div terms */
  const float r_inv = r ? 1.0f / r : 0.0f;
  const float faci = mj * wi_dx * r_inv;

  double dA[3];
  for (int i = 0; i < 3; ++i)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];

  const double dAdr = dA[0] * dx[0] + dA[1] * dx[1] + dA[2] * dx[2];
  pi->mhd_data.divA -= faci * dAdr;
  for (int i = 0; i < 3; ++i)
    pi->mhd_data.BPred[i] += faci * (dA[(i + 1) % 3] * dx[(i + 2) % 3] -
                                     dA[(i + 2) % 3] * dx[(i + 1) % 3]);
}

/**
 * @brief Calculate the MHD-gradient interaction between particle i and particle
 * j
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
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  /* Define kernel variables */
  float wi, wj, wi_dx, wj_dx;
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }
  float Ai[3], Aj[3];
  for (int i = 0; i < 3; ++i) {
    Ai[i] = pi->mhd_data.APred[i];
    Aj[i] = pj->mhd_data.APred[i];
  }

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr_tmp = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr_tmp = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  // const float f_ij = 1.f;
  // const float f_ji = 1.f;
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  const float wi_dr = (wi_dr_tmp * 0.5 + wj_dr_tmp * 0.5 * f_ji / f_ij);
  const float wj_dr = (wj_dr_tmp * 0.5 + wi_dr_tmp * 0.5 * f_ij / f_ji);

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];
  /* A dot r. */
  const float Ari = Ai[0] * dx[0] + Ai[1] * dx[1] + Ai[2] * dx[2];
  const float Arj = Aj[0] * dx[0] + Aj[1] * dx[1] + Aj[2] * dx[2];
  /* Compute gradient terms */
  const float over_rho_i = 1.0f / rhoi * f_ij;
  const float over_rho_j = 1.0f / rhoj * f_ji;

  /* Calculate divergence term */
  float B_mon_i = -over_rho_i * (Bri - Brj) * wi_dr * r_inv;
  float B_mon_j = -over_rho_j * (Bri - Brj) * wj_dr * r_inv;
  pi->mhd_data.divB += mj * B_mon_i;
  pj->mhd_data.divB += mi * B_mon_j;
  /* Calculate divergence term */
  float A_mon_i = -over_rho_i * (Ari - Arj) * wi_dr * r_inv;
  float A_mon_j = -over_rho_j * (Ari - Arj) * wj_dr * r_inv;
  pi->mhd_data.divA += mj * A_mon_i;
  pj->mhd_data.divA += mi * A_mon_j;

  /* dB */
  float dB[3];
  for (int k = 0; k < 3; ++k) dB[k] = Bi[k] - Bj[k];
  /* dA */
  float dA[3];
  for (int k = 0; k < 3; ++k) dA[k] = Ai[k] - Aj[k];

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];
  /* dB cross r */
  float dA_cross_dx[3];
  dA_cross_dx[0] = dA[1] * dx[2] - dA[2] * dx[1];
  dA_cross_dx[1] = dA[2] * dx[0] - dA[0] * dx[2];
  dA_cross_dx[2] = dA[0] * dx[1] - dA[1] * dx[0];

  /* Calculate Curl */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.curl_B[k] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[k];
    pj->mhd_data.curl_B[k] += mi * over_rho_j * wj_dr * r_inv * dB_cross_dx[k];
    pi->mhd_data.BSmooth[k] += mj * over_rho_i * wi_dr * r_inv * dA_cross_dx[k];
    pj->mhd_data.BSmooth[k] += mi * over_rho_j * wj_dr * r_inv * dA_cross_dx[k];
    //    pi->mhd_data.BSmooth[k] += mj * wi * Bi[k];
    //    pj->mhd_data.BSmooth[k] += mi * wj * Bj[k];
  }
  /* calculate the weights */
  pi->mhd_data.Q0 += pj->mass * wi;
  pj->mhd_data.Q0 += pi->mass * wj;
  /* Calculate gradient of B tensor */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.grad_B_tensor[i][j] -=
          mj * over_rho_i * wi_dr * r_inv * dB[i] * dx[j];
      pj->mhd_data.grad_B_tensor[i][j] -=
          mi * over_rho_j * wj_dr * r_inv * dB[i] * dx[j];
    }
  }

  /* Calculate SPH error */
  pi->mhd_data.mean_SPH_err += mj * wi;
  pj->mhd_data.mean_SPH_err += mi * wj;
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.mean_grad_SPH_err[k] +=
        mj * over_rho_i * wi_dr * r_inv * dx[k];
    pj->mhd_data.mean_grad_SPH_err[k] -=
        mi * over_rho_j * wj_dr * r_inv * dx[k];
  }

  /* Vect Potential difference */
  const float Aij[3] = {pj->mhd_data.APred[0] - pi->mhd_data.APred[0],
                        pj->mhd_data.APred[1] - pi->mhd_data.APred[1],
                        pj->mhd_data.APred[2] - pi->mhd_data.APred[2]};
  const float DAi[3] = {pi->mhd_data.APred[0] * (pj->v[0] - pi->v[0]),
                        pi->mhd_data.APred[1] * (pj->v[1] - pi->v[2]),
                        pi->mhd_data.APred[2] * (pj->v[2] - pi->v[2])};
  const float DAj[3] = {pj->mhd_data.APred[0] * (pj->v[0] - pi->v[0]),
                        pj->mhd_data.APred[1] * (pj->v[1] - pi->v[2]),
                        pj->mhd_data.APred[2] * (pj->v[2] - pi->v[2])};

  const float common_term_i = wi * mj / rhoj;

  /* The inverse of the C-matrix. eq. 6
   * It's symmetric so recall we only store the 6 useful terms. */
  pi->mhd_data.grad.c_matrix_inv.xx += common_term_i * dx[0] * dx[0];
  pi->mhd_data.grad.c_matrix_inv.yy += common_term_i * dx[1] * dx[1];
  pi->mhd_data.grad.c_matrix_inv.zz += common_term_i * dx[2] * dx[2];
  pi->mhd_data.grad.c_matrix_inv.xy += common_term_i * dx[0] * dx[1];
  pi->mhd_data.grad.c_matrix_inv.xz += common_term_i * dx[0] * dx[2];
  pi->mhd_data.grad.c_matrix_inv.yz += common_term_i * dx[1] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pi->mhd_data.grad.Mat_bx[0] -= common_term_i * Aij[0] * dx[0];
  pi->mhd_data.grad.Mat_bx[1] -= common_term_i * Aij[0] * dx[1];
  pi->mhd_data.grad.Mat_bx[2] -= common_term_i * Aij[0] * dx[2];

  pi->mhd_data.grad.Mat_by[0] -= common_term_i * Aij[1] * dx[0];
  pi->mhd_data.grad.Mat_by[1] -= common_term_i * Aij[1] * dx[1];
  pi->mhd_data.grad.Mat_by[2] -= common_term_i * Aij[1] * dx[2];

  pi->mhd_data.grad.Mat_bz[0] -= common_term_i * Aij[2] * dx[0];
  pi->mhd_data.grad.Mat_bz[1] -= common_term_i * Aij[2] * dx[1];
  pi->mhd_data.grad.Mat_bz[2] -= common_term_i * Aij[2] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pi->mhd_data.grad.Mat_dax[0] -= common_term_i * DAi[0] * dx[0];
  pi->mhd_data.grad.Mat_dax[1] -= common_term_i * DAi[0] * dx[1];
  pi->mhd_data.grad.Mat_dax[2] -= common_term_i * DAi[0] * dx[2];

  pi->mhd_data.grad.Mat_day[0] -= common_term_i * DAi[1] * dx[0];
  pi->mhd_data.grad.Mat_day[1] -= common_term_i * DAi[1] * dx[1];
  pi->mhd_data.grad.Mat_day[2] -= common_term_i * DAi[1] * dx[2];

  pi->mhd_data.grad.Mat_daz[0] -= common_term_i * DAi[2] * dx[0];
  pi->mhd_data.grad.Mat_daz[1] -= common_term_i * DAi[2] * dx[1];
  pi->mhd_data.grad.Mat_daz[2] -= common_term_i * DAi[2] * dx[2];

  const float common_term_j = wj * mi / rhoi;

  /* The inverse of the C-matrix. eq. 6
   * It's symmetric so recall we only store the 6 useful terms. */
  pj->mhd_data.grad.c_matrix_inv.xx += common_term_j * dx[0] * dx[0];
  pj->mhd_data.grad.c_matrix_inv.yy += common_term_j * dx[1] * dx[1];
  pj->mhd_data.grad.c_matrix_inv.zz += common_term_j * dx[2] * dx[2];
  pj->mhd_data.grad.c_matrix_inv.xy += common_term_j * dx[0] * dx[1];
  pj->mhd_data.grad.c_matrix_inv.xz += common_term_j * dx[0] * dx[2];
  pj->mhd_data.grad.c_matrix_inv.yz += common_term_j * dx[1] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pj->mhd_data.grad.Mat_bx[0] -= common_term_j * Aij[0] * dx[0];
  pj->mhd_data.grad.Mat_bx[1] -= common_term_j * Aij[0] * dx[1];
  pj->mhd_data.grad.Mat_bx[2] -= common_term_j * Aij[0] * dx[2];

  pj->mhd_data.grad.Mat_by[0] -= common_term_j * Aij[1] * dx[0];
  pj->mhd_data.grad.Mat_by[1] -= common_term_j * Aij[1] * dx[1];
  pj->mhd_data.grad.Mat_by[2] -= common_term_j * Aij[1] * dx[2];

  pj->mhd_data.grad.Mat_bz[0] -= common_term_j * Aij[2] * dx[0];
  pj->mhd_data.grad.Mat_bz[1] -= common_term_j * Aij[2] * dx[1];
  pj->mhd_data.grad.Mat_bz[2] -= common_term_j * Aij[2] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pj->mhd_data.grad.Mat_dax[0] -= common_term_j * DAj[0] * dx[0];
  pj->mhd_data.grad.Mat_dax[1] -= common_term_j * DAj[0] * dx[1];
  pj->mhd_data.grad.Mat_dax[2] -= common_term_j * DAj[0] * dx[2];

  pj->mhd_data.grad.Mat_day[0] -= common_term_j * DAj[1] * dx[0];
  pj->mhd_data.grad.Mat_day[1] -= common_term_j * DAj[1] * dx[1];
  pj->mhd_data.grad.Mat_day[2] -= common_term_j * DAj[1] * dx[2];

  pj->mhd_data.grad.Mat_daz[0] -= common_term_j * DAj[2] * dx[0];
  pj->mhd_data.grad.Mat_daz[1] -= common_term_j * DAj[2] * dx[1];
  pj->mhd_data.grad.Mat_daz[2] -= common_term_j * DAj[2] * dx[2];
}

/**
 * @brief Calculate the MHDgradient interaction between particle i and particle
 * j (non-symmetric)
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
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_mhd_gradient(const float r2, const float dx[3],
                                const float hi, const float hj,
                                struct part *restrict pi,
                                const struct part *restrict pj,
                                const float mu_0, const float a,
                                const float H) {

  /* Define kernel variables */
  float wi, wi_dx;
  float wj, wj_dx;
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get the mass. */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }
  float Ai[3], Aj[3];
  for (int i = 0; i < 3; ++i) {
    Ai[i] = pi->mhd_data.APred[i];
    Aj[i] = pj->mhd_data.APred[i];
  }

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr_tmp = hid_inv * wi_dx;

  /* Get the kernel for hi. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr_tmp = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  const float wi_dr = (wi_dr_tmp * 0.5 + wj_dr_tmp * 0.5 * f_ji / f_ij);

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]);
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]);
  /* A dot r. */
  const float Ari = (Ai[0] * dx[0] + Ai[1] * dx[1] + Ai[2] * dx[2]);
  const float Arj = (Aj[0] * dx[0] + Aj[1] * dx[1] + Aj[2] * dx[2]);

  /* Compute gradient terms */
  const float over_rho_i = 1.0f / rhoi * f_ij;

  /* Calculate divergence term */
  float B_mon_i = -over_rho_i * (Bri - Brj) * wi_dr * r_inv;
  pi->mhd_data.divB += mj * B_mon_i;
  /* Calculate divergence term */
  float A_mon_i = -over_rho_i * (Ari - Arj) * wi_dr * r_inv;
  pi->mhd_data.divA += mj * A_mon_i;

  /* dB */
  float dB[3];
  for (int k = 0; k < 3; ++k) dB[k] = Bi[k] - Bj[k];
  /* dA */
  float dA[3];
  for (int k = 0; k < 3; ++k) dA[k] = Ai[k] - Aj[k];

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];
  float dA_cross_dx[3];
  dA_cross_dx[0] = dA[1] * dx[2] - dA[2] * dx[1];
  dA_cross_dx[1] = dA[2] * dx[0] - dA[0] * dx[2];
  dA_cross_dx[2] = dA[0] * dx[1] - dA[1] * dx[0];

  /* Calculate Curl */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.curl_B[k] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[k];
    pi->mhd_data.BSmooth[k] += mj * over_rho_i * wi_dr * r_inv * dA_cross_dx[k];
    //    pi->mhd_data.BSmooth[k] += mj * wi * Bi[k];
  }
  /* calculate the weights */
  pi->mhd_data.Q0 += pj->mass * wi;
  /* Calculate gradient of B tensor */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.grad_B_tensor[i][j] -=
          mj * over_rho_i * wi_dr * r_inv * dB[i] * dx[j];
    }
  }

  /* Calculate SPH error */
  pi->mhd_data.mean_SPH_err += mj * wi;
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.mean_grad_SPH_err[k] +=
        mj * over_rho_i * wi_dr * r_inv * dx[k];
  }

  /* Vect Potential difference */
  const float Aij[3] = {pj->mhd_data.APred[0] - pi->mhd_data.APred[0],
                        pj->mhd_data.APred[1] - pi->mhd_data.APred[1],
                        pj->mhd_data.APred[2] - pi->mhd_data.APred[2]};
  const float DAi[3] = {pi->mhd_data.APred[0] * (pj->v[0] - pi->v[0]),
                        pi->mhd_data.APred[1] * (pj->v[1] - pi->v[2]),
                        pi->mhd_data.APred[2] * (pj->v[2] - pi->v[2])};

  const float common_term = wi * mj / rhoj;

  /* The inverse of the C-matrix. eq. 6
   * It's symmetric so recall we only store the 6 useful terms. */
  pi->mhd_data.grad.c_matrix_inv.xx += common_term * dx[0] * dx[0];
  pi->mhd_data.grad.c_matrix_inv.yy += common_term * dx[1] * dx[1];
  pi->mhd_data.grad.c_matrix_inv.zz += common_term * dx[2] * dx[2];
  pi->mhd_data.grad.c_matrix_inv.xy += common_term * dx[0] * dx[1];
  pi->mhd_data.grad.c_matrix_inv.xz += common_term * dx[0] * dx[2];
  pi->mhd_data.grad.c_matrix_inv.yz += common_term * dx[1] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pi->mhd_data.grad.Mat_bx[0] -= common_term * Aij[0] * dx[0];
  pi->mhd_data.grad.Mat_bx[1] -= common_term * Aij[0] * dx[1];
  pi->mhd_data.grad.Mat_bx[2] -= common_term * Aij[0] * dx[2];

  pi->mhd_data.grad.Mat_by[0] -= common_term * Aij[1] * dx[0];
  pi->mhd_data.grad.Mat_by[1] -= common_term * Aij[1] * dx[1];
  pi->mhd_data.grad.Mat_by[2] -= common_term * Aij[1] * dx[2];

  pi->mhd_data.grad.Mat_bz[0] -= common_term * Aij[2] * dx[0];
  pi->mhd_data.grad.Mat_bz[1] -= common_term * Aij[2] * dx[1];
  pi->mhd_data.grad.Mat_bz[2] -= common_term * Aij[2] * dx[2];

  /* Gradient of v (recall dx is pi - pj), eq. 18 */
  pi->mhd_data.grad.Mat_dax[0] -= common_term * DAi[0] * dx[0];
  pi->mhd_data.grad.Mat_dax[1] -= common_term * DAi[0] * dx[1];
  pi->mhd_data.grad.Mat_dax[2] -= common_term * DAi[0] * dx[2];

  pi->mhd_data.grad.Mat_day[0] -= common_term * DAi[1] * dx[0];
  pi->mhd_data.grad.Mat_day[1] -= common_term * DAi[1] * dx[1];
  pi->mhd_data.grad.Mat_day[2] -= common_term * DAi[1] * dx[2];

  pi->mhd_data.grad.Mat_daz[0] -= common_term * DAi[2] * dx[0];
  pi->mhd_data.grad.Mat_daz[1] -= common_term * DAi[2] * dx[1];
  pi->mhd_data.grad.Mat_daz[2] -= common_term * DAi[2] * dx[2];
}

/**
 * @brief MHD-Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_mhd_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float Pi = pi->force.pressure;
  const float Pj = pj->force.pressure;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr_tmp = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr_tmp = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  // const float f_ij = 1.f;
  // const float f_ji = 1.f;
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
  const float rho_ij = rhoi + rhoj;

  /* Construct the gradient functions (eq. 4 and 5) */
  /*  float G_i[3], G_j[3];
    sym_matrix_multiply_by_vector(G_i, &pi->mhd_data.force.c_matrix, dx);
    sym_matrix_multiply_by_vector(G_j, &pj->mhd_data.force.c_matrix, dx);
  */
  /* Note we multiply by -1 as dx is (pi - pj) and not (pj - pi) */
  /*  G_i[0] *= -wi * hid_inv;
    G_i[1] *= -wi * hid_inv;
    G_i[2] *= -wi * hid_inv;
    G_j[0] *= -wj * hjd_inv;
    G_j[1] *= -wj * hjd_inv;
    G_j[2] *= -wj * hjd_inv;

    const float G_ij[3] = {0.5f * (G_i[0] + G_j[0]), 0.5f * (G_i[1] + G_j[1]),
                           0.5f * (G_i[2] + G_j[2])};

   */
  const float wi_dr = (wi_dr_tmp * 0.5 + wj_dr_tmp * 0.5 * f_ji / f_ij);
  const float wj_dr = (wj_dr_tmp * 0.5 + wi_dr_tmp * 0.5 * f_ij / f_ji);

  const float a_fac =
      pow(a, 2.f * mhd_comoving_factor + 3.f * (hydro_gamma - 1.f));

  const float mag_faci = f_ij * wi_dr * r_inv / (rhoi * rhoi) / mu_0 * a_fac;
  const float mag_facj = f_ji * wj_dr * r_inv / (rhoj * rhoj) / mu_0 * a_fac;

  // const float mag_faci =  1.f / mu_0 * a_fac;
  // const float mag_facj =  1.f / mu_0 * a_fac;
  float Bi[3], Bj[3];
  float mm_i[3][3], mm_j[3][3];

  for (int i = 0; i < 3; i++) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  ///////////////////////////// FORCE MAXWELL TENSOR
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      mm_i[i][j] = Bi[i] * Bi[j];
      mm_j[i][j] = Bj[i] * Bj[j];
    }
  for (int j = 0; j < 3; j++) {
    mm_i[j][j] -= 0.5 * (Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2]);
    mm_j[j][j] -= 0.5 * (Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2]);
  }

  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  const float plasma_beta_j = B2j != 0.0f ? 2.0f * mu_0 * Pj / B2j : FLT_MAX;

  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float scale_j = 0.125f * (10.0f - plasma_beta_j);

  const float tensile_correction_scale_i = fmaxf(0.0f, fminf(scale_i, 1.0f));
  const float tensile_correction_scale_j = fmaxf(0.0f, fminf(scale_j, 1.0f));

  //////////////////////////// Apply to the Force and DIVB TERM SUBTRACTION
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      pi->a_hydro[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pj->a_hydro[i] -=
          mi * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->a_hydro[i] -= mj * Bi[i] * tensile_correction_scale_i *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
      pj->a_hydro[i] += mi * Bj[i] * tensile_correction_scale_j *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
    }

  /* Save forces*/
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.tot_mag_F[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pj->mhd_data.tot_mag_F[i] -=
          mi * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->mhd_data.tot_mag_F[i] -=
          mj * Bi[i] * (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
      pj->mhd_data.tot_mag_F[i] +=
          mi * Bj[i] * (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
    }
  }
  /////////////////////////// VP evolution
  const float mag_VPIndi = f_ij * wi_dr * r_inv / rhoi;
  const float mag_VPIndj = f_ji * wj_dr * r_inv / rhoj;
  // Normal Gauge
  double dA[3];
  for (int i = 0; i < 3; i++)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];
#ifdef VP_ADV_GAUGE
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float sourceAi = dv[0] * pi->mhd_data.APred[0] +
                         dv[1] * pi->mhd_data.APred[1] +
                         dv[2] * pi->mhd_data.APred[2];
  const float sourceAj = dv[0] * pj->mhd_data.APred[0] +
                         dv[1] * pj->mhd_data.APred[1] +
                         dv[2] * pj->mhd_data.APred[2];
#else
  const float sourceAi =
      -(dA[0] * pi->v[0] + dA[1] * pi->v[1] + dA[2] * pi->v[2]);
  const float sourceAj =
      -(dA[0] * pj->v[0] + dA[1] * pj->v[1] + dA[2] * pj->v[2]);
#endif
  float SAi = sourceAi + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  float SAj = sourceAj + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  SAi = a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  SAj = a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);

  for (int i = 0; i < 3; i++) {
    pi->mhd_data.dAdt[i] += mj * mag_VPIndi * SAi * dx[i];
    pj->mhd_data.dAdt[i] += mi * mag_VPIndj * SAj * dx[i];
    // pi->mhd_data.dAdt[i] += mj * SAi * G_ij[i] / rhoj *dx[i];
    // pj->mhd_data.dAdt[i] += mi * SAj * G_ij[i] / rhoi *dx[i];
  }
  /// DISSSIPATION
  const float mag_Disi =
      (f_ij * wi_dr + f_ji * wj_dr) / 2.f * r_inv * rhoi / (rho_ij * rho_ij);
  const float mag_Disj =
      (f_ji * wj_dr + f_ij * wi_dr) / 2.f * r_inv * rhoj / (rho_ij * rho_ij);
  for (int i = 0; i < 3; i++) {
    pi->mhd_data.dAdt[i] +=
        mj * 8.0 * pi->mhd_data.resistive_eta * mag_Disi * dA[i];
    pj->mhd_data.dAdt[i] -=
        mi * 8.0 * pj->mhd_data.resistive_eta * mag_Disj * dA[i];
  }

  /* Artificial resistivity */
  const float alpha_AR = 0.5f * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const float vsig_AR =
      0.5f * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const float rhoij = 0.5f * (rhoi + rhoj);
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1.0f / rhoij2;

  const float grad_term = 0.5f * (f_ij * wi_dr + f_ji * wj_dr);

  const float art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.dAdt[0] += mj * art_diff_pref * dA[0];
  pi->mhd_data.dAdt[1] += mj * art_diff_pref * dA[1];
  pi->mhd_data.dAdt[2] += mj * art_diff_pref * dA[2];

  pj->mhd_data.dAdt[0] -= mi * art_diff_pref * dA[0];
  pj->mhd_data.dAdt[1] -= mi * art_diff_pref * dA[1];
  pj->mhd_data.dAdt[2] -= mi * art_diff_pref * dA[2];

  /* Save induction sources */
  for (int i = 0; i < 3; i++) {
    pi->mhd_data.Adv_A_source[i] += mj * mag_VPIndi * SAi * dx[i];
    pj->mhd_data.Adv_A_source[i] += mi * mag_VPIndj * SAj * dx[i];
    pi->mhd_data.Diff_A_source[i] +=
        mj * 8.0 * pi->mhd_data.resistive_eta * mag_Disi * dA[i];
    pj->mhd_data.Diff_A_source[i] -=
        mi * 8.0 * pj->mhd_data.resistive_eta * mag_Disj * dA[i];
    pi->mhd_data.Delta_A[i] += mj * 8.0 * mag_Disi * dA[i];
    pj->mhd_data.Delta_A[i] -= mi * 8.0 * mag_Disj * dA[i];
  }
}

/**
 * @brief MHD-Force interaction between two particles. non-symmetric version.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 The vaccuum permeability constant in internal units.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_mhd_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  /* Cosmological factors entering the EoMs */
  // const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  // const float a2_Hubble = a * a * H;

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float Pi = pi->force.pressure;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr_tmp = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr_tmp = hjd_inv * wj_dx;

  /* Construct the gradient functions (eq. 4 and 5) */
  /*  float G_i[3], G_j[3];
    sym_matrix_multiply_by_vector(G_i, &pi->mhd_data.force.c_matrix, dx);
    sym_matrix_multiply_by_vector(G_j, &pj->mhd_data.force.c_matrix, dx);
  */
  /* Note we multiply by -1 as dx is (pi - pj) and not (pj - pi) */
  /*  G_i[0] *= -wi * hid_inv;
    G_i[1] *= -wi * hid_inv;
    G_i[2] *= -wi * hid_inv;
    G_j[0] *= -wj * hjd_inv;
    G_j[1] *= -wj * hjd_inv;
    G_j[2] *= -wj * hjd_inv;

    const float G_ij[3] = {0.5f * (G_i[0] + G_j[0]), 0.5f * (G_i[1] + G_j[1]),
                           0.5f * (G_i[2] + G_j[2])};

  */
  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
  const float rho_ij = rhoi + rhoj;

  const float wi_dr = (wi_dr_tmp * 0.5 + wj_dr_tmp * 0.5 * f_ji / f_ij);
  const float wj_dr = (wj_dr_tmp * 0.5 + wi_dr_tmp * 0.5 * f_ij / f_ji);

  const float a_fac =
      pow(a, 2.f * mhd_comoving_factor + 3.f * (hydro_gamma - 1.f));

  const float mag_faci = f_ij * wi_dr * r_inv / (rhoi * rhoi) / mu_0 * a_fac;
  const float mag_facj = f_ji * wj_dr * r_inv / (rhoj * rhoj) / mu_0 * a_fac;
  // const float mag_faci =  1.f / mu_0 * a_fac;
  // const float mag_facj =  1.f / mu_0 * a_fac;
  float Bi[3], Bj[3];
  float mm_i[3][3], mm_j[3][3];

  for (int i = 0; i < 3; i++) {
    Bi[i] = pi->mhd_data.BPred[i];
    Bj[i] = pj->mhd_data.BPred[i];
  }
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];

  ///////////////////////////// FORCE MAXWELL TENSOR
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      mm_i[i][j] = Bi[i] * Bi[j];
      mm_j[i][j] = Bj[i] * Bj[j];
    }
  for (int j = 0; j < 3; j++) {
    mm_i[j][j] -= 0.5 * (Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2]);
    mm_j[j][j] -= 0.5 * (Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2]);
  }

  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float tensile_correction_scale_i = fmaxf(0.0f, fminf(scale_i, 1.0f));

  //////////////////////////// Apply to the Force and DIVB TERM SUBTRACTION
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      pi->a_hydro[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->a_hydro[i] -= mj * Bi[i] * tensile_correction_scale_i *
                        (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
      // pi->a_hydro[i] +=
      //     mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * G_ij[j] /
      //     rhoj * dx[j];
      // pi->a_hydro[i] -= mj * Bi[i] * tensile_correction_scale_i *
      //                   (Bi[j] * mag_faci + Bj[j] * mag_facj) * G_ij[j] /
      //                   rhoi * dx[j];
    }

  /* Save forces*/
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.tot_mag_F[i] +=
          mj * (mm_i[i][j] * mag_faci + mm_j[i][j] * mag_facj) * dx[j];
      pi->mhd_data.tot_mag_F[i] -=
          mj * Bi[i] * (Bi[j] * mag_faci + Bj[j] * mag_facj) * dx[j];
    }
  }

  /////////////////////////// VP INDUCTION
  const float mag_VPIndi = f_ij * wi_dr * r_inv / rhoi;
  // Normal Gauge
  double dA[3];
  for (int i = 0; i < 3; i++)
    dA[i] = pi->mhd_data.APred[i] - pj->mhd_data.APred[i];
#ifdef VP_ADV_GAUGE
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float sourceAi = dv[0] * pi->mhd_data.APred[0] +
                         dv[1] * pi->mhd_data.APred[1] +
                         dv[2] * pi->mhd_data.APred[2];
#else
  const float sourceAi =
      -(dA[0] * pi->v[0] + dA[1] * pi->v[1] + dA[2] * pi->v[2]);
#endif
  float SAi = sourceAi + a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  SAi = a * a * (pi->mhd_data.Gau - pj->mhd_data.Gau);
  for (int i = 0; i < 3; i++) {
    pi->mhd_data.dAdt[i] += mj * mag_VPIndi * SAi * dx[i];
    //  pi->mhd_data.dAdt[i] += mj * SAi * G_ij[i] / rhoj * dx[i];
  }
  /// DISSSIPATION
  const float mag_Disi =
      (f_ij * wi_dr + f_ji * wj_dr) / 2.f * r_inv * rhoi / (rho_ij * rho_ij);
  for (int i = 0; i < 3; i++)
    pi->mhd_data.dAdt[i] +=
        mj * 8.0 * pi->mhd_data.resistive_eta * mag_Disi * dA[i];

  /* Artificial resistivity */
  const float alpha_AR = 0.5f * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const float vsig_AR =
      0.5f * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const float rhoij = 0.5f * (rhoi + rhoj);
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1.0f / rhoij2;

  const float grad_term = 0.5f * (f_ij * wi_dr + f_ji * wj_dr);

  const float art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.dAdt[0] += mj * art_diff_pref * dA[0];
  pi->mhd_data.dAdt[1] += mj * art_diff_pref * dA[1];
  pi->mhd_data.dAdt[2] += mj * art_diff_pref * dA[2];

  /* Save induction sources */
  for (int i = 0; i < 3; i++) {
    pi->mhd_data.Adv_A_source[i] += mj * mag_VPIndi * SAi * dx[i];
    pi->mhd_data.Diff_A_source[i] +=
        mj * 8.0 * pi->mhd_data.resistive_eta * mag_Disi * dA[i];
    pi->mhd_data.Delta_A[i] += mj * 8.0 * mag_Disi * dA[i];
  }
}
#endif /* SWIFT_VECTOR_POTENTIAL_MHD_H */
