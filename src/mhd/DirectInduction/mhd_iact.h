/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DIRECT_INDUCTION_MHD_IACT_H
#define SWIFT_DIRECT_INDUCTION_MHD_IACT_H

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

  /* Compute kernels */
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  const float hj_inv = 1.f / hj;
  const float xj = r * hj_inv;

  kernel_deval(xi, &wi, &wi_dx);
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute neighbour weight norm */
  pi->mhd_data.norm_Nw += 1;
  pi->mhd_data.norm_Kw += wi;
  pj->mhd_data.norm_Nw += 1;
  pj->mhd_data.norm_Kw += wj;

  pi->mhd_data.hb_over_ha_Nw += hj;
  pi->mhd_data.hb_over_ha_Kw += hj * wi;
  pj->mhd_data.hb_over_ha_Nw += hi;
  pj->mhd_data.hb_over_ha_Kw += hi * wj;

  /* Compute weighted distance sum */
  for (int k = 0; k < 3; k++) {
  
      pi->mhd_data.rcm_over_ha_Nw[k] -= dx[k];
      pi->mhd_data.rcm_over_ha_Kw[k] -= dx[k] * wi;

      pj->mhd_data.rcm_over_ha_Nw[k] += dx[k];
      pj->mhd_data.rcm_over_ha_Kw[k] += dx[k] * wj;

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
                               const struct part *restrict pj, const float mu_0,
                               const float a, const float H) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Compute kernels */
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  const float hj_inv = 1.f / hj;
  const float xj = r * hj_inv;

  kernel_deval(xi, &wi, &wi_dx);
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute neighbour weight norm */
  pi->mhd_data.norm_Nw += 1;
  pi->mhd_data.norm_Kw += wi;

  pi->mhd_data.hb_over_ha_Nw += hj;
  pi->mhd_data.hb_over_ha_Kw += hj * wi;

  /* Compute weighted distance sum */
  for (int k = 0; k < 3; k++) {
  
      pi->mhd_data.rcm_over_ha_Nw[k] -= dx[k] ;
      pi->mhd_data.rcm_over_ha_Kw[k] -= dx[k] * wi;

  }



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

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.B_over_rho[i] * rhoi;
    Bj[i] = pj->mhd_data.B_over_rho[i] * rhoj;
  }

  float dB[3];
  for (int i = 0; i < 3; ++i) dB[i] = Bi[i] - Bj[i];

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float over_rho_i = 1.0f / rhoi * f_ij;
  const float over_rho_j = 1.0f / rhoj * f_ji;


  /* Error corrections */
  float dx_corr_gradi[3];
  float dx_corr_gradj[3];
  for (int ki = 0; ki < 3; ki++) {
    dx_corr_gradi[ki] = 0.0f;
    dx_corr_gradj[ki] = 0.0f;
    for (int kj = 0; kj < 3; kj++) {
       dx_corr_gradi[ki] += pi->err_proj_tensor[ki][kj] * dx[kj]
       dx_corr_gradj[ki] += pj->err_proj_tensor[ki][kj] * dx[kj]
     }
  }
  float dB_cross_dx_corri[3];
  float dB_cross_dx_corrj[3];
  dB_cross_dx_corri[0] = dB[1] * dx_corr_gradi[2] - dB[2] * dx_corr_gradi[1];
  dB_cross_dx_corri[1] = dB[2] * dx_corr_gradi[0] - dB[0] * dx_corr_gradi[2];
  dB_cross_dx_corri[2] = dB[0] * dx_corr_gradi[1] - dB[1] * dx_corr_gradi[0];
  dB_cross_dx_corrj[0] = dB[1] * dx_corr_gradj[2] - dB[2] * dx_corr_gradj[1];
  dB_cross_dx_corrj[1] = dB[2] * dx_corr_gradj[0] - dB[0] * dx_corr_gradj[2];
  dB_cross_dx_corrj[2] = dB[0] * dx_corr_gradj[1] - dB[1] * dx_corr_gradj[0];


  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[0]-dB_cross_dx_corri[0]);
  pi->mhd_data.curl_B[1] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[1]-dB_cross_dx_corri[1]);
  pi->mhd_data.curl_B[2] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[2]-dB_cross_dx_corri[2]);
  pj->mhd_data.curl_B[0] += mi * over_rho_j * wj_dr * r_inv * (dB_cross_dx[0]-dB_cross_dx_corrj[0]);
  pj->mhd_data.curl_B[1] += mi * over_rho_j * wj_dr * r_inv * (dB_cross_dx[1]-dB_cross_dx_corrj[1]);
  pj->mhd_data.curl_B[2] += mi * over_rho_j * wj_dr * r_inv * (dB_cross_dx[2]-dB_cross_dx_corrj[2]);

  /* Calculate gradient of B tensor */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.grad_B_tensor[i][j] -=
          mj * over_rho_i * wi_dr * r_inv * dB[i] * (dx[j] - dx_corr_gradi[j]);
      pj->mhd_data.grad_B_tensor[i][j] -=
          mi * over_rho_j * wj_dr * r_inv * dB[i] * (dx[j] - dx_corr_gradj[j]);
    }
  }

  /* Calculate SPH error */
  pi->mhd_data.mean_SPH_err += mj * wi;
  pj->mhd_data.mean_SPH_err += mi * wj;

  pi->mhd_data.norm_SPHw += mj / rhoj * wi;
  pj->mhd_data.norm_SPHw += mi / rhoi * wj;

  pi->mhd_data.hb_over_ha_SPHw += hj * mj / rhoj * wi;
  pj->mhd_data.hb_over_ha_SPHw += hi * mi / rhoi * wj;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.mean_grad_SPH_err[k] +=
        mj * over_rho_i * wi_dr * r_inv * dx[k];
    pi->mhd_data.rcm_over_ha_SPHw[k] -= dx[k] * mj / rhoj * wi;
    pi->mhd_data.symmetric_gradient_err_fij[k] += mj * (wi_dr * f_ij / (rhoi * rhoi) + wj_dr * f_ji / (rhoj * rhoj)) * r_inv * dx[k];

    pj->mhd_data.mean_grad_SPH_err[k] -=
        mi * over_rho_j * wj_dr * r_inv * dx[k];
    pj->mhd_data.rcm_over_ha_SPHw[k] += dx[k] * mi / rhoi * wj;
    pj->mhd_data.symmetric_gradient_err_fij[k] -= mi * (wj_dr * f_ji / (rhoj * rhoj) + wi_dr * f_ij / (rhoi * rhoi)) * r_inv * dx[k];


  }

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
  float wi, wj, wi_dx, wj_dx;
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.B_over_rho[i] * rhoi;
    Bj[i] = pj->mhd_data.B_over_rho[i] * rhoj;
  }

  float dB[3];
  for (int i = 0; i < 3; ++i) dB[i] = Bi[i] - Bj[i];

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float over_rho_i = 1.0f / rhoi * f_ij;

  /* Error corrections */
  float dx_corr_gradi[3];
  for (int ki = 0; ki < 3; ki++) {
    dx_corr_gradi[ki] = 0.0f;
    for (int kj = 0; kj < 3; kj++) {
       dx_corr_gradi[ki] += pi->err_proj_tensor[ki][kj] * dx[kj]
     }
  }
  float dB_cross_dx_corri[3];
  dB_cross_dx_corri[0] = dB[1] * dx_corr_gradi[2] - dB[2] * dx_corr_gradi[1];
  dB_cross_dx_corri[1] = dB[2] * dx_corr_gradi[0] - dB[0] * dx_corr_gradi[2];
  dB_cross_dx_corri[2] = dB[0] * dx_corr_gradi[1] - dB[1] * dx_corr_gradi[0];

  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[0]-dB_cross_dx_corri[0]);
  pi->mhd_data.curl_B[1] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[1]-dB_cross_dx_corri[1]);
  pi->mhd_data.curl_B[2] += mj * over_rho_i * wi_dr * r_inv * (dB_cross_dx[2]-dB_cross_dx_corri[2]);

  /* Calculate gradient of B tensor */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->mhd_data.grad_B_tensor[i][j] -=
          mj * over_rho_i * wi_dr * r_inv * dB[i] * (dx[j] - dx_corr_gradi[j]);
    }
  }

  /* Calculate SPH error */
  pi->mhd_data.mean_SPH_err += mj * wi;

  pi->mhd_data.norm_SPHw += mj / rhoj * wi;

  pi->mhd_data.hb_over_ha_SPHw += hj * mj / rhoj * wi;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.mean_grad_SPH_err[k] +=
        mj * over_rho_i * wi_dr * r_inv * dx[k];
    pi->mhd_data.rcm_over_ha_SPHw[k] -= dx[k] * mj / rhoj * wi;
    pi->mhd_data.symmetric_gradient_err_fij[k] += mj * (wi_dr * f_ij / (rhoi * rhoi) + wj_dr * f_ji / (rhoj * rhoj)) * r_inv * dx[k];

  }


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

  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  float Bi[3];
  float Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const float permeability_inv = 1.0f / mu_0;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  float f_ij = 1.f - pi->force.f / mj;
  float f_ji = 1.f - pj->force.f / mi;

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];
  const float dBdr = Bri - Brj;

  /* Error corrections */
  float dx_corr_gradi[3];
  float dx_corr_gradj[3];
  for (int ki = 0; ki < 3; ki++) {
    dx_corr_gradi[ki] = 0.0f;
    dx_corr_gradj[ki] = 0.0f;
    for (int kj = 0; kj < 3; kj++) {
       dx_corr_gradi[ki] += pi->err_proj_tensor[ki][kj] * dx[kj]
       dx_corr_gradj[ki] += pj->err_proj_tensor[ki][kj] * dx[kj]
     }
  }
  const float Bri_corr_gradi = Bi[0] * dx_corr_gradi[0] + Bi[1] * dx_corr_gradi[1] + Bi[2] * dx_corr_gradi[2];
  const float Bri_corr_gradj = Bi[0] * dx_corr_gradj[0] + Bi[1] * dx_corr_gradj[1] + Bi[2] * dx_corr_gradj[2];
  const float Brj_corr_gradi = Bj[0] * dx_corr_gradi[0] + Bj[1] * dx_corr_gradi[1] + Bj[2] * dx_corr_gradi[2];
  const float Brj_corr_gradj = Bj[0] * dx_corr_gradj[0] + Bj[1] * dx_corr_gradj[1] + Bj[2] * dx_corr_gradj[2];
  const float dBdr_corr_gradi = Bri_corr_gradi - Brj_corr_gradi;
  const float dBdr_corr_gradj = Brj_corr_gradj - Bri_corr_gradj;

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Compute symmetric div(B) */
  const float grad_term_i = over_rho2_i * wi_dr * r_inv;
  const float grad_term_j = over_rho2_j * wj_dr * r_inv;

  const float asym_grad_term_i = f_ij * wi_dr * r_inv / rhoi;
  const float asym_grad_term_j = f_ji * wj_dr * r_inv / rhoj;

  const float divB_i = (dBdr-dBdr_corr_gradi) * asym_grad_term_i;
  const float divB_j = (dBdr-dBdr_corr_gradj) * asym_grad_term_j;

  pi->mhd_data.divB -= mj * divB_i;
  pj->mhd_data.divB -= mi * divB_j;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[0]-dx_corr_gradi[0]);
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[0]-dx_corr_gradj[0]);

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * (Bri - Bri_corr_gradi) * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * (Brj - Brj_corr_gradj) * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[1]-dx_corr_gradi[1]);
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[1]-dx_corr_gradj[1]);

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * (Bri - Bri_corr_gradi) * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * (Brj - Brj_corr_gradj) * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[2]-dx_corr_gradi[2]);
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[2]-dx_corr_gradj[2]);

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * (Bri - Bri_corr_gradi) * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * (Brj - Brj_corr_gradj) * permeability_inv * r_inv * Bj[2];

  /* SPH acceleration term in x direction, j_th particle */
  float sph_acc_term_j[3];
  sph_acc_term_j[0] = -sph_acc_term_i[0];
  sph_acc_term_j[1] = -sph_acc_term_i[1];
  sph_acc_term_j[2] = -sph_acc_term_i[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  const float monopole_beta = pi->mhd_data.monopole_beta;

  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  const float plasma_beta_j = B2j != 0.0f ? 2.0f * mu_0 * Pj / B2j : FLT_MAX;

  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float scale_j = 0.125f * (10.0f - plasma_beta_j);

  float tensile_correction_scale_i = fmaxf(0.0f, fminf(scale_i, 1.0f));
  float tensile_correction_scale_j = fmaxf(0.0f, fminf(scale_j, 1.0f));

  /* Mitigation: switching tensile correction to 0.5 */
  //tensile_correction_scale_i = fminf(tensile_correction_scale_i, 0.5f * ( 1.0f + pi->mhd_data.mhdsw) );
  //tensile_correction_scale_j = fminf(tensile_correction_scale_j, 0.5f * ( 1.0f + pj->mhd_data.mhdsw) );


  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bi[0] * tensile_correction_scale_i;
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bi[0] * tensile_correction_scale_i;

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bi[1] * tensile_correction_scale_i;
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bi[1] * tensile_correction_scale_i;

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bi[2] * tensile_correction_scale_i;
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bi[2] * tensile_correction_scale_i;

  sph_acc_term_j[0] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bj[0] * tensile_correction_scale_j;
  sph_acc_term_j[0] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bj[0] * tensile_correction_scale_j;

  sph_acc_term_j[1] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bj[1] * tensile_correction_scale_j;
  sph_acc_term_j[1] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bj[1] * tensile_correction_scale_j;

  sph_acc_term_j[2] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri - Bri_corr_gradi) * r_inv * Bj[2] * tensile_correction_scale_j;
  sph_acc_term_j[2] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj - Brj_corr_gradj) * r_inv * Bj[2] * tensile_correction_scale_j;


  /* Mitigation: switching off force component parallel to the chosen error variable */
/*
  float err_vec_i[3];
  float err_vec_j[3];
  for (int k = 0; k < 3; k++) {  
    err_vec_i[k] = pi->mhd_data.rcm_over_ha_Kw[k];
    err_vec_j[k] = pj->mhd_data.rcm_over_ha_Kw[k];
    // err_vec_i[k] = pi->mhd_data.symmetric_gradient_err_fij[k];
    // err_vec_j[k] = pi->mhd_data.symmetric_gradient_err_fij[k];
  }
  const float sph_acc_term_i_abs_pr = ( err_vec_i[0] * sph_acc_term_i[0] + err_vec_i[1] * sph_acc_term_i[1] + err_vec_i[2] * sph_acc_term_i[2] ) / (err_vec_i[0]*err_vec_i[0] + err_vec_i[1]*err_vec_i[1] + err_vec_i[2]*err_vec_i[2]);
  const float sph_acc_term_j_abs_pr = ( err_vec_j[0] * sph_acc_term_j[0] + err_vec_j[1] * sph_acc_term_j[1] + err_vec_j[2] * sph_acc_term_j[2] ) / (err_vec_j[0]*err_vec_j[0] + err_vec_j[1]*err_vec_j[1] + err_vec_j[2]*err_vec_j[2]);
  for (int k = 0; k < 3; k++) {  
    sph_acc_term_i[k] -= sph_acc_term_i_abs_pr * err_vec_i[k] * (1.0f -  pi->mhd_data.mhdsw);
    sph_acc_term_j[k] -= sph_acc_term_j_abs_pr * err_vec_j[k] * (1.0f -  pj->mhd_data.mhdsw);
  }
*/
  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2];

  pj->a_hydro[0] -= mi * sph_acc_term_j[0];
  pj->a_hydro[1] -= mi * sph_acc_term_j[1];
  pj->a_hydro[2] -= mi * sph_acc_term_j[2];

  /* Save forces */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
    pj->mhd_data.tot_mag_F[k] -= mi * sph_acc_term_j[k];
  }

  /* Direct Induction */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;
  const float dB_dt_pref_j = over_rho2_j * wj_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = - (Bri - Bri_corr_gradi) * dv[0];
  dB_dt_i[1] = - (Bri - Bri_corr_gradi) * dv[1];
  dB_dt_i[2] = - (Bri - Bri_corr_gradi) * dv[2];

  float dB_dt_j[3];
  dB_dt_j[0] = - (Brj - Brj_corr_gradj) * dv[0];
  dB_dt_j[1] = - (Brj - Brj_corr_gradj)* dv[1];
  dB_dt_j[2] = - (Brj - Brj_corr_gradj) * dv[2];

  /* */
  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  pj->mhd_data.B_over_rho_dt[0] += mi * dB_dt_pref_j * dB_dt_j[0];
  pj->mhd_data.B_over_rho_dt[1] += mi * dB_dt_pref_j * dB_dt_j[1];
  pj->mhd_data.B_over_rho_dt[2] += mi * dB_dt_pref_j * dB_dt_j[2];

  /* Physical resistivity */
  const float resistive_eta_i = pi->mhd_data.resistive_eta;
  const float resistive_eta_j = pj->mhd_data.resistive_eta;

  const float rho_term_PR = 1.0f / (rhoi * rhoj);
  const float grad_term_PR = f_ij * wi_dr + f_ji * wj_dr;

  const float dB_dt_pref_PR = rho_term_PR * grad_term_PR * r_inv;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_over_rho_dt[k] +=
        resistive_eta_i * mj * dB_dt_pref_PR * dB[k];
    pj->mhd_data.B_over_rho_dt[k] -=
        resistive_eta_j * mi * dB_dt_pref_PR * dB[k];
  }

  pi->u_dt -=
      0.5f * permeability_inv * resistive_eta_i * mj * dB_dt_pref_PR * dB_2;
  pj->u_dt -=
      0.5f * permeability_inv * resistive_eta_j * mi * dB_dt_pref_PR * dB_2;

  /* Artificial resistivity */
  const float alpha_AR = 0.5f * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const float vsig_AR =
      0.5f * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const float rhoij = 0.5f * (rhoi + rhoj);
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1.0f / rhoij2;

  const float grad_term = 0.5f * (f_ij * wi_dr + f_ji * wj_dr);

  const float art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref * dB[2];

  pj->mhd_data.B_over_rho_dt[0] -= mi * art_diff_pref * dB[0];
  pj->mhd_data.B_over_rho_dt[1] -= mi * art_diff_pref * dB[1];
  pj->mhd_data.B_over_rho_dt[2] -= mi * art_diff_pref * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_diff_pref * dB_2;
  pj->u_dt -= 0.5f * mi * permeability_inv * art_diff_pref * dB_2;

  /* Store AR terms */
  pi->mhd_data.B_over_rho_dt_AR[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt_AR[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt_AR[2] += mj * art_diff_pref * dB[2];

  pj->mhd_data.B_over_rho_dt_AR[0] -= mi * art_diff_pref * dB[0];
  pj->mhd_data.B_over_rho_dt_AR[1] -= mi * art_diff_pref * dB[1];
  pj->mhd_data.B_over_rho_dt_AR[2] -= mi * art_diff_pref * dB[2];

  pi->mhd_data.u_dt_AR -= 0.5f * mj * permeability_inv * art_diff_pref * dB_2;
  pj->mhd_data.u_dt_AR -= 0.5f * mi * permeability_inv * art_diff_pref * dB_2;

  /* Divergence diffusion */
  const float vsig_Dedner_i = 0.5f * pi->viscosity.v_sig;
  const float vsig_Dedner_j = 0.5f * pj->viscosity.v_sig;

  float grad_psi = grad_term_i * psi_over_ch_i * vsig_Dedner_i;
  grad_psi += grad_term_j * psi_over_ch_j * vsig_Dedner_j;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi * (dx[0]-dx_corr_gradi[0]);
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi * (dx[1]-dx_corr_gradi[1]);
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi * (dx[2]-dx_corr_gradi[2]);

  pj->mhd_data.B_over_rho_dt[0] += mi * grad_psi * (dx[0]-dx_corr_gradi[0]);
  pj->mhd_data.B_over_rho_dt[1] += mi * grad_psi * (dx[1]-dx_corr_gradi[1]);
  pj->mhd_data.B_over_rho_dt[2] += mi * grad_psi * (dx[2]-dx_corr_gradi[2]);

  /* Save induction sources */
  const float dB_dt_pref_Lap_i = 2.0f * r_inv / rhoj;
  const float dB_dt_pref_Lap_j = 2.0f * r_inv / rhoi;

  for (int i = 0; i < 3; i++) {
    pi->mhd_data.Adv_B_source[i] += mj * dB_dt_pref_i * dB_dt_i[i];
    pj->mhd_data.Adv_B_source[i] += mi * dB_dt_pref_j * dB_dt_j[i];
    pi->mhd_data.Adv_B_source[i] -= mj * grad_psi * dx[i];
    pj->mhd_data.Adv_B_source[i] += mi * grad_psi * dx[i];
    pi->mhd_data.Diff_B_source[i] +=
        mj * resistive_eta_i * mj * dB_dt_pref_PR * dB[i];
    pj->mhd_data.Diff_B_source[i] -=
        mi * resistive_eta_j * mi * dB_dt_pref_PR * dB[i];
    pi->mhd_data.Diff_B_source[i] += mj * art_diff_pref * dB[i];
    pj->mhd_data.Diff_B_source[i] -= mi * art_diff_pref * dB[i];
    pi->mhd_data.Delta_B[i] += mj * dB_dt_pref_Lap_i * wi_dr * dB[i];
    pj->mhd_data.Delta_B[i] -= mi * dB_dt_pref_Lap_j * wj_dr * dB[i];
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float Pi = pi->force.pressure;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  float Bi[3];
  float Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const float permeability_inv = 1.0f / mu_0;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  float f_ij = 1.f - pi->force.f / mj;
  float f_ji = 1.f - pj->force.f / mi;

  /* Mitigation: switching fij to 0 */
  //f_ij *= pi->mhd_data.mhdsw;
  //f_ji *= pj->mhd_data.mhdsw;

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];
  const float dBdr = Bri - Brj;

  /* Error corrections */
  float dx_corr_gradi[3];
  float dx_corr_gradj[3];
  for (int ki = 0; ki < 3; ki++) {
    dx_corr_gradi[ki] = 0.0f;
    dx_corr_gradj[ki] = 0.0f;
    for (int kj = 0; kj < 3; kj++) {
       dx_corr_gradi[ki] += pi->err_proj_tensor[ki][kj] * dx[kj]
       dx_corr_gradj[ki] += pj->err_proj_tensor[ki][kj] * dx[kj]
     }
  }
  const float Bri_corr_gradi = Bi[0] * dx_corr_gradi[0] + Bi[1] * dx_corr_gradi[1] + Bi[2] * dx_corr_gradi[2];
  const float Bri_corr_gradj = Bi[0] * dx_corr_gradj[0] + Bi[1] * dx_corr_gradj[1] + Bi[2] * dx_corr_gradj[2];
  const float Brj_corr_gradi = Bj[0] * dx_corr_gradi[0] + Bj[1] * dx_corr_gradi[1] + Bj[2] * dx_corr_gradi[2];
  const float Brj_corr_gradj = Bj[0] * dx_corr_gradj[0] + Bj[1] * dx_corr_gradj[1] + Bj[2] * dx_corr_gradj[2];
  const float dBdr_corr_gradi = Bri_corr_gradi - Brj_corr_gradi;


  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Compute symmetric div(B) */
  const float grad_term_i = over_rho2_i * wi_dr * r_inv;
  const float grad_term_j = over_rho2_j * wj_dr * r_inv;

  const float asym_grad_term_i = f_ij * wi_dr * r_inv / rhoi;

  const float divB_i = dBdr * asym_grad_term_i;

  pi->mhd_data.divB -= mj * divB_i;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[0]-dx_corr_gradi[0]);
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[0]-dx_corr_gradj[0]);

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * (Bri-Bri_corr_gradi) * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * (Brj-Brj_corr_gradj) * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[1]-dx_corr_gradi[1]);
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[1]-dx_corr_gradj[1]);

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * (Bri-Bri_corr_gradi) * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * (Brj-Brj_corr_gradj) * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * (dx[2]-dx_corr_gradi[2]);
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * (dx[2]-dx_corr_gradj[2]);

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * (Bri-Bri_corr_gradi) * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * (Brj-Brj_corr_gradj) * permeability_inv * r_inv * Bj[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  const float monopole_beta = pi->mhd_data.monopole_beta;

  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  float tensile_correction_scale_i = fmaxf(0.0f, fminf(scale_i, 1.0f));

  /* Mitigation: switching tensile correction to 0.5 */
  //tensile_correction_scale_i = fminf(tensile_correction_scale_i, 0.5f * ( 1.0f + pi->mhd_data.mhdsw) );

  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri-Bri_corr_gradi) * r_inv * Bi[0] * tensile_correction_scale_i;
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj-Brj_corr_gradj) * r_inv * Bi[0] * tensile_correction_scale_i;

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri-Bri_corr_gradi) * r_inv * Bi[1] * tensile_correction_scale_i;
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj-Brj_corr_gradj) * r_inv * Bi[1] * tensile_correction_scale_i;

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       (Bri-Bri_corr_gradi) * r_inv * Bi[2] * tensile_correction_scale_i;
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       (Brj-Brj_corr_gradj) * r_inv * Bi[2] * tensile_correction_scale_i;

  /* Mitigation: switching off force component parallel to the chosen error variable */
/*
  float err_vec_i[3];
  for (int k = 0; k < 3; k++) {  
    err_vec_i[k] = pi->mhd_data.rcm_over_ha_Kw[k];
    // err_vec_i[k] = pi->mhd_data.symmetric_gradient_err_fij[k];
  }
  const float sph_acc_term_i_abs_pr = ( err_vec_i[0] * sph_acc_term_i[0] + err_vec_i[1] * sph_acc_term_i[1] + err_vec_i[2] * sph_acc_term_i[2] ) / (err_vec_i[0]*err_vec_i[0] + err_vec_i[1]*err_vec_i[1] + err_vec_i[2]*err_vec_i[2]);
  for (int k = 0; k < 3; k++) {  
    sph_acc_term_i[k] -= sph_acc_term_i_abs_pr * err_vec_i[k] * (1.0f -  pi->mhd_data.mhdsw);
  }
*/


  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2];

  /* Save forces */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
  }

  /* */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = - (Bri-Bri_corr_gradi) * dv[0];
  dB_dt_i[1] = - (Bri-Bri_corr_gradi) * dv[1];
  dB_dt_i[2] = - (Bri-Bri_corr_gradi) * dv[2];

  /* */
  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  /* Physical resistivity */
  const float resistive_eta_i = pi->mhd_data.resistive_eta;

  const float rho_term_PR = 1.0f / (rhoi * rhoj);
  const float grad_term_PR = f_ij * wi_dr + f_ji * wj_dr;

  const float dB_dt_pref_PR = rho_term_PR * grad_term_PR * r_inv;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_over_rho_dt[k] +=
        resistive_eta_i * mj * dB_dt_pref_PR * dB[k];
  }

  pi->u_dt -=
      0.5f * permeability_inv * resistive_eta_i * mj * dB_dt_pref_PR * dB_2;

  /* Artificial resistivity */
  const float alpha_AR = 0.5f * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const float vsig_AR =
      0.5f * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const float rhoij = 0.5f * (rhoi + rhoj);
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1.0f / rhoij2;

  const float grad_term = 0.5f * (f_ij * wi_dr + f_ji * wj_dr);

  const float art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_diff_pref * dB_2;

  /* Store AR terms */
  pi->mhd_data.B_over_rho_dt_AR[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt_AR[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt_AR[2] += mj * art_diff_pref * dB[2];

  pi->mhd_data.u_dt_AR -= 0.5f * mj * permeability_inv * art_diff_pref * dB_2;

  /* Divergence diffusion */
  const float vsig_Dedner_i = 0.5f * pi->viscosity.v_sig;
  const float vsig_Dedner_j = 0.5f * pj->viscosity.v_sig;

  float grad_psi = grad_term_i * psi_over_ch_i * vsig_Dedner_i;
  grad_psi += grad_term_j * psi_over_ch_j * vsig_Dedner_j;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi * (dx[0]-dx_corr_gradi[0]);
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi * (dx[1]-dx_corr_gradi[1]);
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi * (dx[2]-dx_corr_gradi[2]);

  /* Save induction sources */
  const float dB_dt_pref_Lap = 2.0f * r_inv / rhoj;

  for (int i = 0; i < 3; i++) {
    pi->mhd_data.Adv_B_source[i] += mj * dB_dt_pref_i * dB_dt_i[i];
    pi->mhd_data.Adv_B_source[i] -= mj * grad_psi * dx[i];
    pi->mhd_data.Diff_B_source[i] +=
        resistive_eta_i * mj * dB_dt_pref_PR * dB[i];
    pi->mhd_data.Diff_B_source[i] += mj * art_diff_pref * dB[i];
    pi->mhd_data.Delta_B[i] += mj * dB_dt_pref_Lap * wi_dr * dB[i];
  }
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
