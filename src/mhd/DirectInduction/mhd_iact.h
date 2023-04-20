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

extern float monopole_beta;
extern float resistivity_beta;

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
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const float a, const float H) {}

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
                               const float H) {}

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
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const float a, const float H) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3];
  float Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

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
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* dB cross r */
  float dB_cross_dx[3]; 
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Calculate monopole term */
  float B_mon_i = -over_rho2_i * rhoi * (Bri - Brj) * wi_dr * r_inv;
  float B_mon_j = -over_rho2_j * rhoj * (Bri - Brj) * wj_dr * r_inv;
  pi->mhd_data.B_mon += mj * B_mon_i;
  pj->mhd_data.B_mon += mi * B_mon_j;

  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[0];
  pi->mhd_data.curl_B[1] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[1];
  pi->mhd_data.curl_B[2] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[2];
  pj->mhd_data.curl_B[0] += mi * over_rho2_j * rhoj * wj_dr * r_inv * dB_cross_dx[0];
  pj->mhd_data.curl_B[1] += mi * over_rho2_j * rhoj * wj_dr * r_inv * dB_cross_dx[1];
  pj->mhd_data.curl_B[2] += mi * over_rho2_j * rhoj * wj_dr * r_inv * dB_cross_dx[2];
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
                                const double mu_0, const float a,
                                const float H) {

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3];
  float Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  // const float hj_inv = 1.0f / hj;
  // const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  // const float xj = r * hj_inv;
  // float wj, wj_dx;
  // kernel_deval(xj, &wj, &wj_dx);
  // const float wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  // const float f_ji = 1.f - pj->force.f / mi;

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]);
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]);

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  // const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Calculate monopole term */
  float B_mon_i = -over_rho2_i * rhoi * (Bri - Brj) * wi_dr * r_inv;
  pi->mhd_data.B_mon += mj * B_mon_i;

  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[0];
  pi->mhd_data.curl_B[1] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[1];
  pi->mhd_data.curl_B[2] += mj * over_rho2_i * rhoi * wi_dr * r_inv * dB_cross_dx[2];
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
    struct part *restrict pi, struct part *restrict pj, 
    const double mu_0, const float a, const float H) {

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

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const float permeability_inv = 1.f / mu_0;

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
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Isotropic pressure */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* SPH acceleration term in x direction, j_th particle */
  float sph_acc_term_j[3];
  sph_acc_term_j[0] = -sph_acc_term_i[0];
  sph_acc_term_j[1] = -sph_acc_term_i[1];
  sph_acc_term_j[2] = -sph_acc_term_i[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  // const float monopole_beta = hydro_props->mhd.monopole_subtraction;

  const float plasma_beta_i = 2.0f * mu_0 * Pi / B2i;
  const float plasma_beta_j = 2.0f * mu_0 * Pj / B2j;
  
  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float scale_j = 0.125f * (10.0f - plasma_beta_j);

  const float tensile_correction_scale_i = fmax(0.0f, fmin(scale_i, 1.0f));
  const float tensile_correction_scale_j = fmax(0.0f, fmin(scale_j, 1.0f));

  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[0] * tensile_correction_scale_i;
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[0] * tensile_correction_scale_i;

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[1] * tensile_correction_scale_i;
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[1] * tensile_correction_scale_i;

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[2] * tensile_correction_scale_i;
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[2] * tensile_correction_scale_i;

  sph_acc_term_j[0] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[0] * tensile_correction_scale_j;
  sph_acc_term_j[0] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[0] * tensile_correction_scale_j;

  sph_acc_term_j[1] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[1] * tensile_correction_scale_j;
  sph_acc_term_j[1] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[1] * tensile_correction_scale_j;

  sph_acc_term_j[2] -= monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bj[2] * tensile_correction_scale_j;
  sph_acc_term_j[2] -= monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bj[2] * tensile_correction_scale_j;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2];

  pj->a_hydro[0] -= mi * sph_acc_term_j[0];
  pj->a_hydro[1] -= mi * sph_acc_term_j[1];
  pj->a_hydro[2] -= mi * sph_acc_term_j[2];

  /* Calculate monopole term */
  // float B_mon_i = -over_rho2_i * rhoi * (Bri - Brj) * wi_dr * r_inv;
  // float B_mon_j = -over_rho2_j * rhoj * (Bri - Brj) * wj_dr * r_inv;
  // pi->mhd_data.B_mon += mj * B_mon_i;
  // pj->mhd_data.B_mon += mi * B_mon_j;

  /* */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;
  const float dB_dt_pref_j = over_rho2_j * wj_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = -Bri * dv[0];
  dB_dt_i[1] = -Bri * dv[1];
  dB_dt_i[2] = -Bri * dv[2];

  float dB_dt_j[3];
  dB_dt_j[0] = -Brj * dv[0];
  dB_dt_j[1] = -Brj * dv[1];
  dB_dt_j[2] = -Brj * dv[2];

  /* */
  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  pj->mhd_data.B_over_rho_dt[0] += mi * dB_dt_pref_j * dB_dt_j[0];
  pj->mhd_data.B_over_rho_dt[1] += mi * dB_dt_pref_j * dB_dt_j[1];
  pj->mhd_data.B_over_rho_dt[2] += mi * dB_dt_pref_j * dB_dt_j[2];

  /*Artificial resistivity*/

  // const float resistivity_beta = hydro_props->mhd.art_resistivity;

  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                          dv_cross_dx[1] * dv_cross_dx[1] +
                          dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B = sqrtf(v_sig_B_2) * r_inv;

  const float art_res_pref = 0.5f * resistivity_beta * v_sig_B *
                             (wi_dr * over_rho2_i + wj_dr * over_rho2_j);

  /*
  const float div_err_i = hi * fabs(pi->mhd_data.B_mon) / (sqrtf(B2i) + 1.e-18);
  const float alphares_i = min(div_err_i, 1.0f);

  const float div_err_j = hj * fabs(pj->mhd_data.B_mon) / (sqrtf(B2j) + 1.e-18);
  const float alphares_j = min(div_err_j, 1.0f);
  */
 
  /*
  const float cs2_i = pi->force.soundspeed * pi->force.soundspeed;
  const float cs2_j = pj->force.soundspeed * pj->force.soundspeed;

  const float vA2_i = B2i * permeability_inv / rhoi;
  const float vA2_j = B2j * permeability_inv / rhoj;

  const float vsig_AR2_i = cs2_i + vA2_i;
  const float vsig_AR2_j = cs2_j + vA2_j;

  const float vsig_AR_i = sqrtf(vsig_AR2_i);
  const float vsig_AR_j = sqrtf(vsig_AR2_j);
  */  

  /*
  const float corr_AR_i = -4.0f * cs2_i * Bri * Bri * r_inv * r_inv *
  permeability_inv /rhoi; const float corr_AR_j = -4.0f * cs2_j * Brj * Brj *
  r_inv * r_inv * permeability_inv /rhoj;

  const float vsig_B_i = sqrtf(0.5f * vsig_AR2_i + 0.5f * sqrtf(vsig_AR2_i +
  corr_AR_i)); const float vsig_B_j = sqrtf(0.5f * vsig_AR2_j + 0.5f *
  sqrtf(vsig_AR2_j + corr_AR_j));

  const float art_res_pref_i = alphares_i * vsig_B_i * dB_dt_pref_i;
  const float art_res_pref_j = alphares_j * vsig_B_j * dB_dt_pref_j;

  const float art_res_pref = resistivity_beta * 0.5f * (art_res_pref_i +
  art_res_pref_j);
  */

  pi->mhd_data.B_over_rho_dt[0] += mj * art_res_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_res_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_res_pref * dB[2];

  pj->mhd_data.B_over_rho_dt[0] -= mi * art_res_pref * dB[0];
  pj->mhd_data.B_over_rho_dt[1] -= mi * art_res_pref * dB[1];
  pj->mhd_data.B_over_rho_dt[2] -= mi * art_res_pref * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_res_pref * dB_2;
  pj->u_dt -= 0.5f * mi * permeability_inv * art_res_pref * dB_2;

  /*Divergence diffusion */
  
  const float vsig_Dedner_i = pi->viscosity.v_sig;
  const float vsig_Dedner_j = pj->viscosity.v_sig;
   
  float grad_psi_i = over_rho2_i * psi_over_ch_i * vsig_Dedner_i * wi_dr * r_inv;
  grad_psi_i += over_rho2_j * psi_over_ch_j * vsig_Dedner_j * wj_dr * r_inv;
  float grad_psi_j = grad_psi_i;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi_i * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi_i * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi_i * dx[2];

  pj->mhd_data.B_over_rho_dt[0] += mi * grad_psi_j * dx[0];
  pj->mhd_data.B_over_rho_dt[1] += mi * grad_psi_j * dx[1];
  pj->mhd_data.B_over_rho_dt[2] += mi * grad_psi_j * dx[2];
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
    struct part *restrict pi, const struct part *restrict pj,  
    const double mu_0, const float a, const float H) {

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
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;

  /* Isotropic pressure */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5f * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5f * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      -1.f * over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      -1.f * over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  // const float monopole_beta = hydro_props->mhd.monopole_subtraction;

  const float plasma_beta_i = 2.0f * mu_0 * Pi / B2i;

  const float scale_i = 0.125f * (10.0f - plasma_beta_i);

  const float tensile_correction_scale_i = fmax(0.0f, fmin(scale_i, 1.0f));

  sph_acc_term_i[0] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[0] * tensile_correction_scale_i;
  sph_acc_term_i[0] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[0] * tensile_correction_scale_i;

  sph_acc_term_i[1] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[1] * tensile_correction_scale_i;
  sph_acc_term_i[1] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[1] * tensile_correction_scale_i;

  sph_acc_term_i[2] += monopole_beta * over_rho2_i * wi_dr * permeability_inv *
                       Bri * r_inv * Bi[2] * tensile_correction_scale_i;
  sph_acc_term_i[2] += monopole_beta * over_rho2_j * wj_dr * permeability_inv *
                       Brj * r_inv * Bi[2] * tensile_correction_scale_i;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2];

  /* Calculate monopole term */
  // float B_mon_i = -over_rho2_i * rhoi * (Bri - Brj) * wi_dr * r_inv;
  // pi->mhd_data.B_mon += mj * B_mon_i;

  /* */
  const float dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;
  // const float dB_dt_pref_j = over_rho2_j * wj_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = -Bri * dv[0];
  dB_dt_i[1] = -Bri * dv[1];
  dB_dt_i[2] = -Bri * dv[2];

  /* */
  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  /*Artificial resistivity*/
  
  // const float resistivity_beta = hydro_props->mhd.art_resistivity;

  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                          dv_cross_dx[1] * dv_cross_dx[1] +
                          dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B = sqrtf(v_sig_B_2) * r_inv;

  const float art_res_pref = 0.5f * resistivity_beta * v_sig_B *
                             (wi_dr * over_rho2_i + wj_dr * over_rho2_j);

  /*
  const float div_err_i = hi * fabs(pi->mhd_data.B_mon) / (sqrtf(B2i) + 1.e-18);
  const float alphares_i = min(div_err_i, 1.0f);

  const float div_err_j = hj * fabs(pj->mhd_data.B_mon) / (sqrtf(B2j) + 1.e-18);
  const float alphares_j = min(div_err_j, 1.0f);
  */
  
  /*
  const float cs2_i = pi->force.soundspeed * pi->force.soundspeed;
  const float cs2_j = pj->force.soundspeed * pj->force.soundspeed;

  const float vA2_i = B2i * permeability_inv / rhoi;
  const float vA2_j = B2j * permeability_inv / rhoj;

  const float vsig_AR2_i = cs2_i + vA2_i;
  const float vsig_AR2_j = cs2_j + vA2_j;

  const float vsig_AR_i = sqrtf(vsig_AR2_i);
  const float vsig_AR_j = sqrtf(vsig_AR2_j);
  */   

  /*
  const float corr_AR_i = -4.0f * cs2_i * Bri * Bri * r_inv * r_inv *
  permeability_inv /rhoi; const float corr_AR_j = -4.0f * cs2_j * Brj * Brj *
  r_inv * r_inv * permeability_inv /rhoj;

  const float vsig_B_i = sqrtf(0.5f * vsig_AR2_i + 0.5f * sqrtf(vsig_AR2_i +
  corr_AR_i)); const float vsig_B_j = sqrtf(0.5f * vsig_AR2_j + 0.5f *
  sqrtf(vsig_AR2_j + corr_AR_j));

  const float art_res_pref_i = alphares_i * vsig_B_i * dB_dt_pref_i;
  const float art_res_pref_j = alphares_j * vsig_B_j * dB_dt_pref_j;

  const float art_res_pref = resistivity_beta * 0.5f * (art_res_pref_i +
  art_res_pref_j);
  */

  pi->mhd_data.B_over_rho_dt[0] += mj * art_res_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_res_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_res_pref * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_res_pref * dB_2;

  /*Divergence diffusion */

  const float vsig_Dedner_i = pi->viscosity.v_sig;
  const float vsig_Dedner_j = pj->viscosity.v_sig;
      
  float grad_psi_i = over_rho2_i * psi_over_ch_i * vsig_Dedner_i * wi_dr * r_inv;
  grad_psi_i += over_rho2_j * psi_over_ch_j * vsig_Dedner_j * wj_dr * r_inv;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi_i * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi_i * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi_i * dx[2];
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
