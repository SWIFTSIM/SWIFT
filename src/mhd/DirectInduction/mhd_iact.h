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

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float over_rho_i = 1.0f / rhoi * f_ij;
  const float over_rho_j = 1.0f / rhoj * f_ji;

  /* Calculate monopole term */
  float divB_i = -over_rho_i * (Bri - Brj) * wi_dr * r_inv;
  float divB_j = -over_rho_j * (Bri - Brj) * wj_dr * r_inv;
  pi->mhd_data.divB += mj * divB_i;
  pj->mhd_data.divB += mi * divB_j;

  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[0];
  pi->mhd_data.curl_B[1] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[1];
  pi->mhd_data.curl_B[2] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[2];
  pj->mhd_data.curl_B[0] += mi * over_rho_j * wj_dr * r_inv * dB_cross_dx[0];
  pj->mhd_data.curl_B[1] += mi * over_rho_j * wj_dr * r_inv * dB_cross_dx[1];
  pj->mhd_data.curl_B[2] += mi * over_rho_j * wj_dr * r_inv * dB_cross_dx[2];

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

  /* Define kernel variables */
  float wi, wi_dx;
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  // const float mi = pi->mass;
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
  const float over_rho_i = 1.0f / rhoi * f_ij;

  /* Calculate monopole term */
  float divB_i = -over_rho_i * (Bri - Brj) * wi_dr * r_inv;
  pi->mhd_data.divB += mj * divB_i;

  /* Calculate curl */
  pi->mhd_data.curl_B[0] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[0];
  pi->mhd_data.curl_B[1] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[1];
  pi->mhd_data.curl_B[2] += mj * over_rho_i * wi_dr * r_inv * dB_cross_dx[2];

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
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
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

  const float normBi = sqrtf(B2i);
  const float normBj = sqrtf(B2j);

  /*
  float curlBi[3];
  float curlBj[3];
  curlBi[0] = pi->mhd_data.curl_B[0];
  curlBi[1] = pi->mhd_data.curl_B[1];
  curlBi[2] = pi->mhd_data.curl_B[2];
  curlBj[0] = pj->mhd_data.curl_B[0];
  curlBj[1] = pj->mhd_data.curl_B[1];
  curlBj[2] = pj->mhd_data.curl_B[2];
  */

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const double permeability_inv = 1.f / mu_0;

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
  const float monopole_beta = pi->mhd_data.monopole_beta;

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

  /* Save forces */
  for (int k = 1; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
    pj->mhd_data.tot_mag_F[k] -= mi * sph_acc_term_j[k];
  }

  /* Direct Induction */
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

  /* Physical resistivity */
  const float resistive_eta_i = pi->mhd_data.resistive_eta;
  const float resistive_eta_j = pj->mhd_data.resistive_eta;
  const float dB_dt_pref_PR_i = 2.0f * resistive_eta_i * r_inv / (rhoi * rhoj);
  const float dB_dt_pref_PR_j = 2.0f * resistive_eta_j * r_inv / (rhoi * rhoj);

  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_PR_i * wi_dr * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_PR_i * wi_dr * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_PR_i * wi_dr * dB[2];

  pj->mhd_data.B_over_rho_dt[0] -= mi * dB_dt_pref_PR_j * wj_dr * dB[0];
  pj->mhd_data.B_over_rho_dt[1] -= mi * dB_dt_pref_PR_j * wj_dr * dB[1];
  pj->mhd_data.B_over_rho_dt[2] -= mi * dB_dt_pref_PR_j * wj_dr * dB[2];

  /*
  float curlB_cross_dxi[3];
  float curlB_cross_dxj[3];

  curlB_cross_dxi[0] = curlBi[1] * dx[2] - curlBi[2] * dx[1];
  curlB_cross_dxi[1] = curlBi[2] * dx[0] - curlBi[0] * dx[2];
  curlB_cross_dxi[2] = curlBi[0] * dx[1] - curlBi[1] * dx[0];

  curlB_cross_dxj[0] = curlBj[1] * dx[2] - curlBj[2] * dx[1];
  curlB_cross_dxj[1] = curlBj[2] * dx[0] - curlBj[0] * dx[2];
  curlB_cross_dxj[2] = curlBj[0] * dx[1] - curlBj[1] * dx[0];

  pi->mhd_data.B_over_rho_dt[0] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[0];
  pi->mhd_data.B_over_rho_dt[0] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[0];
  pi->mhd_data.B_over_rho_dt[1] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[1];
  pi->mhd_data.B_over_rho_dt[1] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[1];
  pi->mhd_data.B_over_rho_dt[2] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[2];
  pi->mhd_data.B_over_rho_dt[2] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[2];

  pj->mhd_data.B_over_rho_dt[0] -=
      mi * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[0];
  pj->mhd_data.B_over_rho_dt[0] -=
      mi * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[0];
  pj->mhd_data.B_over_rho_dt[1] -=
      mi * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[1];
  pj->mhd_data.B_over_rho_dt[1] -=
      mi * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[1];
  pj->mhd_data.B_over_rho_dt[2] -=
      mi * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[2];
  pj->mhd_data.B_over_rho_dt[2] -=
      mi * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[2];
  */

  /*Artificial resistivity*/

  // const float resistivity_beta = hydro_props->mhd.art_resistivity;
  const float art_diff_beta_i = pi->mhd_data.art_diff_beta;
  const float art_diff_beta_j = pj->mhd_data.art_diff_beta;

  /*
  const float rhoij = rhoi + rhoj;
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1 / rhoij2;

  const float alpha_ARij = pi->mhd_data.alpha_AR  + pj->mhd_data.alpha_AR;
  const float v_sig_Bij  = mhd_get_fast_magnetosonic_wave_speed(dx, pi, a, mu_0)
  + mhd_get_fast_magnetosonic_wave_speed(dx, pj, a, mu_0);

  const float art_res_prefi = resistivity_beta * alpha_ARij * v_sig_Bij *
  rhoij2_inv * wi_dr; const float art_res_prefj = resistivity_beta * alpha_ARij
  * v_sig_Bij * rhoij2_inv * wj_dr;
  */

  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                          dv_cross_dx[1] * dv_cross_dx[1] +
                          dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B = sqrtf(v_sig_B_2) * r_inv;

  const float art_diff_pref_i = 0.5f * art_diff_beta_i * v_sig_B *
                                (wi_dr * over_rho2_i + wj_dr * over_rho2_j);
  const float art_diff_pref_j = 0.5f * art_diff_beta_j * v_sig_B *
                                (wi_dr * over_rho2_i + wj_dr * over_rho2_j);

  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref_i * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref_i * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref_i * dB[2];

  pj->mhd_data.B_over_rho_dt[0] -= mi * art_diff_pref_j * dB[0];
  pj->mhd_data.B_over_rho_dt[1] -= mi * art_diff_pref_j * dB[1];
  pj->mhd_data.B_over_rho_dt[2] -= mi * art_diff_pref_j * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_diff_pref_i * dB_2;
  pj->u_dt -= 0.5f * mi * permeability_inv * art_diff_pref_j * dB_2;

  /*Saving AR contribution*/
  pi->mhd_data.ARIS[0] += mj * art_diff_pref_i * dB[0];
  pi->mhd_data.ARIS[1] += mj * art_diff_pref_i * dB[1];
  pi->mhd_data.ARIS[2] += mj * art_diff_pref_i * dB[2];

  pj->mhd_data.ARIS[0] -= mi * art_diff_pref_j * dB[0];
  pj->mhd_data.ARIS[1] -= mi * art_diff_pref_j * dB[1];
  pj->mhd_data.ARIS[2] -= mi * art_diff_pref_j * dB[2];

  /*Divergence diffusion */

  // const float vsig_Dedner_i = pi->viscosity.v_sig;
  // const float vsig_Dedner_j = pj->viscosity.v_sig;

  const float vsig_Dedner_i = mhd_get_magnetosonic_speed(pi, a, mu_0);
  const float vsig_Dedner_j = mhd_get_magnetosonic_speed(pj, a, mu_0);

  float grad_psi_i =
      over_rho2_i * psi_over_ch_i * vsig_Dedner_i * wi_dr * r_inv;
  grad_psi_i += over_rho2_j * psi_over_ch_j * vsig_Dedner_j * wj_dr * r_inv;
  float grad_psi_j = grad_psi_i;

  const float psi_over_ch_i_inv =
      psi_over_ch_i != 0.f ? 1.f / psi_over_ch_i : 0.;
  const float psi_over_ch_j_inv =
      psi_over_ch_j != 0.f ? 1.f / psi_over_ch_j : 0.;

  const float corr_ratio_i = fabsf(normBi * psi_over_ch_i_inv);
  const float corr_ratio_j = fabsf(normBj * psi_over_ch_j_inv);

  const float Qi = corr_ratio_i < 2 ? 0.5 * corr_ratio_i : 1.0f;
  const float Qj = corr_ratio_j < 2 ? 0.5 * corr_ratio_j : 1.0f;

  pi->mhd_data.B_over_rho_dt[0] -= mj * Qi * grad_psi_i * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * Qi * grad_psi_i * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * Qi * grad_psi_i * dx[2];

  pj->mhd_data.B_over_rho_dt[0] += mi * Qj * grad_psi_j * dx[0];
  pj->mhd_data.B_over_rho_dt[1] += mi * Qj * grad_psi_j * dx[1];
  pj->mhd_data.B_over_rho_dt[2] += mi * Qj * grad_psi_j * dx[2];
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
    struct part *restrict pi, const struct part *restrict pj, const double mu_0,
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

  const float normBi = sqrtf(B2i);

  /*
  float curlBi[3];
  float curlBj[3];
  curlBi[0] = pi->mhd_data.curl_B[0];
  curlBi[1] = pi->mhd_data.curl_B[1];
  curlBi[2] = pi->mhd_data.curl_B[2];
  curlBj[0] = pj->mhd_data.curl_B[0];
  curlBj[1] = pj->mhd_data.curl_B[1];
  curlBj[2] = pj->mhd_data.curl_B[2];
  */

  float dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const double permeability_inv = 1.0f / mu_0;

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
  const float monopole_beta = pi->mhd_data.monopole_beta;

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

  /* Save forces */
  for (int k = 1; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
  }

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

  /* Physical resistivity */
  const float resistive_eta = pi->mhd_data.resistive_eta;

  const float dB_dt_pref_PR = 2.0f * resistive_eta * r_inv / (rhoi * rhoj);

  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_PR * wi_dr * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_PR * wi_dr * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_PR * wi_dr * dB[2];

  /*
  float curlB_cross_dxi[3];
  float curlB_cross_dxj[3];

  curlB_cross_dxi[0] = curlBi[1] * dx[2] - curlBi[2] * dx[1];
  curlB_cross_dxi[1] = curlBi[2] * dx[0] - curlBi[0] * dx[2];
  curlB_cross_dxi[2] = curlBi[0] * dx[1] - curlBi[1] * dx[0];

  curlB_cross_dxj[0] = curlBj[1] * dx[2] - curlBj[2] * dx[1];
  curlB_cross_dxj[1] = curlBj[2] * dx[0] - curlBj[0] * dx[2];
  curlB_cross_dxj[2] = curlBj[0] * dx[1] - curlBj[1] * dx[0];

  pi->mhd_data.B_over_rho_dt[0] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[0];
  pi->mhd_data.B_over_rho_dt[0] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[0];
  pi->mhd_data.B_over_rho_dt[1] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[1];
  pi->mhd_data.B_over_rho_dt[1] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[1];
  pi->mhd_data.B_over_rho_dt[2] +=
      mj * resistive_eta * over_rho2_i * wi_dr * r_inv * curlB_cross_dxi[2];
  pi->mhd_data.B_over_rho_dt[2] +=
      mj * resistive_eta * over_rho2_j * wj_dr * r_inv * curlB_cross_dxj[2];
  */

  /*Artificial resistivity*/

  // const float resistivity_beta = hydro_props->mhd.art_resistivity;

  const float art_diff_beta = pi->mhd_data.art_diff_beta;

  /*
  const float rhoij = rhoi + rhoj;
  const float rhoij2 = rhoij * rhoij;
  const float rhoij2_inv = 1 / rhoij2;

  const float alpha_ARij = pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR;
  const float v_sig_Bij  = mhd_get_fast_magnetosonic_wave_speed(dx, pi, a, mu_0)
  + mhd_get_fast_magnetosonic_wave_speed(dx, pj, a, mu_0);

  const float art_res_pref = resistivity_beta * alpha_ARij * v_sig_Bij *
  rhoij2_inv * wi_dr;
  */

  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                          dv_cross_dx[1] * dv_cross_dx[1] +
                          dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B = sqrtf(v_sig_B_2) * r_inv;

  const float art_diff_pref = 0.5f * art_diff_beta * v_sig_B *
                              (wi_dr * over_rho2_i + wj_dr * over_rho2_j);


  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref * dB[2];

  pi->u_dt -= 0.5f * mj * permeability_inv * art_diff_pref * dB_2;

  /*Saving AR contribution*/
  pi->mhd_data.ARIS[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.ARIS[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.ARIS[2] += mj * art_diff_pref * dB[2];


  /*Divergence diffusion */

  // const float vsig_Dedner_i = pi->viscosity.v_sig;
  // const float vsig_Dedner_j = pj->viscosity.v_sig;

  const float vsig_Dedner_i = mhd_get_magnetosonic_speed(pi, a, mu_0);
  const float vsig_Dedner_j = mhd_get_magnetosonic_speed(pj, a, mu_0);

  float grad_psi_i =
      over_rho2_i * psi_over_ch_i * vsig_Dedner_i * wi_dr * r_inv;
  grad_psi_i += over_rho2_j * psi_over_ch_j * vsig_Dedner_j * wj_dr * r_inv;

  const float psi_over_ch_i_inv =
      psi_over_ch_i != 0.f ? 1.f / psi_over_ch_i : 0.;

  const float corr_ratio_i = fabsf(normBi * psi_over_ch_i_inv);

  const float Qi = corr_ratio_i < 2 ? 0.5 * corr_ratio_i : 1.0f;

  pi->mhd_data.B_over_rho_dt[0] -= mj * Qi * grad_psi_i * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * Qi * grad_psi_i * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * Qi * grad_psi_i * dx[2];
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
