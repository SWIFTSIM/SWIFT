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
    const double r2, const float dx[3], const double hi, const double hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const double a, const double H) {}

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
runner_iact_nonsym_mhd_density(const double r2, const float dx[3],
                               const double hi, const double hj,
                               struct part *restrict pi,
                               const struct part *restrict pj, const double mu_0,
                               const double a, const double H) {}

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
    const double r2, const float dx[3], const double hi, const double hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const double a, const double H) {

  /* Define kernel variables */
  float wi, wj, wi_dx, wj_dx;
  /* Get r and 1/r. */
  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mi = pi->mass;
  const double mj = pj->mass;
  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  double Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.B_over_rho[i] * rhoi;
    Bj[i] = pj->mhd_data.B_over_rho[i] * rhoj;
  }

  double dB[3];
  for (int i = 0; i < 3; ++i) dB[i] = Bi[i] - Bj[i];

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const double hj_inv = 1.0 / hj;
  const double hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const double xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const double wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* dB cross r */
  double dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const double over_rho_i = 1.0 / rhoi * f_ij;
  const double over_rho_j = 1.0 / rhoj * f_ji;

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
runner_iact_nonsym_mhd_gradient(const double r2, const float dx[3],
                                const double hi, const double hj,
                                struct part *restrict pi,
                                const struct part *restrict pj,
                                const double mu_0, const double a,
                                const double H) {

  /* Define kernel variables */
  float wi, wi_dx;
  /* Get r and 1/r. */
  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mj = pj->mass;
  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  double Bi[3], Bj[3];
  for (int i = 0; i < 3; ++i) {
    Bi[i] = pi->mhd_data.B_over_rho[i] * rhoi;
    Bj[i] = pj->mhd_data.B_over_rho[i] * rhoj;
  }

  double dB[3];
  for (int i = 0; i < 3; ++i) dB[i] = Bi[i] - Bj[i];

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;

  /* dB cross r */
  double dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const double over_rho_i = 1.0 / rhoi * f_ij;

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
    const double r2, const float dx[3], const double hi, const double hj,
    struct part *restrict pi, struct part *restrict pj, const double mu_0,
    const double a, const double H) {

  /* Get r and 1/r. */
  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mi = pi->mass;
  const double mj = pj->mass;
  const double Pi = pi->force.pressure;
  const double Pj = pj->force.pressure;
  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  double dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  double Bi[3];
  double Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  const double B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const double B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  double dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const double dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const double psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const double psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const double permeability_inv = 1.0 / mu_0;

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const double hj_inv = 1.0 / hj;
  const double hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const double xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const double wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* B dot r. */
  const double Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const double Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];
  const double dBdr = Bri - Brj;

  /* Compute gradient terms */
  const double over_rho2_i = 1.0 / (rhoi * rhoi) * f_ij;
  const double over_rho2_j = 1.0 / (rhoj * rhoj) * f_ji;

  /* Compute symmetric div(B) */
  const double grad_term_i = over_rho2_i * wi_dr * r_inv;
  const double grad_term_j = over_rho2_j * wj_dr * r_inv;

  const double asym_grad_term_i = f_ij * wi_dr * r_inv / rhoi;
  const double asym_grad_term_j = f_ji * wj_dr * r_inv / rhoj;

  const double divB_i = dBdr * asym_grad_term_i;
  const double divB_j = dBdr * asym_grad_term_j;

  pi->mhd_data.divB -= mj * divB_i;
  pj->mhd_data.divB -= mi * divB_j;

  /* SPH acceleration term in x direction, i_th particle */
  double sph_acc_term_i[3] = {0.0, 0.0, 0.0};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* SPH acceleration term in x direction, j_th particle */
  double sph_acc_term_j[3];
  sph_acc_term_j[0] = -sph_acc_term_i[0];
  sph_acc_term_j[1] = -sph_acc_term_i[1];
  sph_acc_term_j[2] = -sph_acc_term_i[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  const double monopole_beta = pi->mhd_data.monopole_beta;

  const double plasma_beta_i = B2i != 0.0 ? 2.0 * mu_0 * Pi / B2i : 10.0;
  const double plasma_beta_j = B2j != 0.0 ? 2.0 * mu_0 * Pj / B2j : 10.0;

  const double scale_i = 0.125 * (10.0 - plasma_beta_i);
  const double scale_j = 0.125 * (10.0 - plasma_beta_j);

  const double tensile_correction_scale_i = fmax(0.0, fmin(scale_i, 1.0));
  const double tensile_correction_scale_j = fmax(0.0, fmin(scale_j, 1.0));

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
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
    pj->mhd_data.tot_mag_F[k] -= mi * sph_acc_term_j[k];
  }

  /* Direct Induction */
  const double dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;
  const double dB_dt_pref_j = over_rho2_j * wj_dr * r_inv;

  /* */
  double dB_dt_i[3];
  dB_dt_i[0] = -Bri * dv[0];
  dB_dt_i[1] = -Bri * dv[1];
  dB_dt_i[2] = -Bri * dv[2];

  double dB_dt_j[3];
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
  const double resistive_eta_i = pi->mhd_data.resistive_eta;
  const double resistive_eta_j = pj->mhd_data.resistive_eta;

  const double rho_term_PR = 1.0 / (rhoi * rhoj);
  const double grad_term_PR = f_ij * wi_dr + f_ji * wj_dr;

  const double dB_dt_pref_PR = rho_term_PR * grad_term_PR * r_inv;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_over_rho_dt[k] +=
        resistive_eta_i * mj * dB_dt_pref_PR * dB[k];
    pj->mhd_data.B_over_rho_dt[k] -=
        resistive_eta_j * mi * dB_dt_pref_PR * dB[k];
  }

  pi->u_dt -=
      0.5 * permeability_inv * resistive_eta_i * mj * dB_dt_pref_PR * dB_2;
  pj->u_dt -=
      0.5 * permeability_inv * resistive_eta_j * mi * dB_dt_pref_PR * dB_2;

  /* Artificial resistivity */
  const double alpha_AR = 0.5 * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const double vsig_AR =
      0.5 * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const double rhoij = 0.5 * (rhoi + rhoj);
  const double rhoij2 = rhoij * rhoij;
  const double rhoij2_inv = 1.0 / rhoij2;

  const double grad_term = 0.5 * (f_ij * wi_dr + f_ji * wj_dr);

  const double art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref * dB[2];

  pj->mhd_data.B_over_rho_dt[0] -= mi * art_diff_pref * dB[0];
  pj->mhd_data.B_over_rho_dt[1] -= mi * art_diff_pref * dB[1];
  pj->mhd_data.B_over_rho_dt[2] -= mi * art_diff_pref * dB[2];

  pi->u_dt -= 0.5 * mj * permeability_inv * art_diff_pref * dB_2;
  pj->u_dt -= 0.5 * mi * permeability_inv * art_diff_pref * dB_2;

  /* Store AR terms */
  pi->mhd_data.B_over_rho_dt_AR[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt_AR[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt_AR[2] += mj * art_diff_pref * dB[2];

  pj->mhd_data.B_over_rho_dt_AR[0] -= mi * art_diff_pref * dB[0];
  pj->mhd_data.B_over_rho_dt_AR[1] -= mi * art_diff_pref * dB[1];
  pj->mhd_data.B_over_rho_dt_AR[2] -= mi * art_diff_pref * dB[2];

  pi->mhd_data.u_dt_AR -= 0.5 * mj * permeability_inv * art_diff_pref * dB_2;
  pj->mhd_data.u_dt_AR -= 0.5 * mi * permeability_inv * art_diff_pref * dB_2;

  /* Divergence diffusion */
  const double vsig_Dedner_i = 0.5 * pi->viscosity.v_sig;
  const double vsig_Dedner_j = 0.5 * pj->viscosity.v_sig;

  double grad_psi = grad_term_i * psi_over_ch_i * vsig_Dedner_i;
  grad_psi += grad_term_j * psi_over_ch_j * vsig_Dedner_j;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi * dx[2];

  pj->mhd_data.B_over_rho_dt[0] += mi * grad_psi * dx[0];
  pj->mhd_data.B_over_rho_dt[1] += mi * grad_psi * dx[1];
  pj->mhd_data.B_over_rho_dt[2] += mi * grad_psi * dx[2];

  /* Save induction sources */
  const double dB_dt_pref_Lap_i = 2.0 * r_inv / rhoj;
  const double dB_dt_pref_Lap_j = 2.0 * r_inv / rhoi;

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
    const double r2, const float dx[3], const double hi, const double hj,
    struct part *restrict pi, const struct part *restrict pj, const double mu_0,
    const double a, const double H) {

  /* Get r and 1/r. */
  const double r = sqrt(r2);
  const double r_inv = r ? 1.0 / r : 0.0;

  /* Recover some data */
  const double mi = pi->mass;
  const double mj = pj->mass;
  const double Pi = pi->force.pressure;
  const double rhoi = pi->rho;
  const double rhoj = pj->rho;

  double dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  double Bi[3];
  double Bj[3];
  Bi[0] = pi->mhd_data.B_over_rho[0] * rhoi;
  Bi[1] = pi->mhd_data.B_over_rho[1] * rhoi;
  Bi[2] = pi->mhd_data.B_over_rho[2] * rhoi;
  Bj[0] = pj->mhd_data.B_over_rho[0] * rhoj;
  Bj[1] = pj->mhd_data.B_over_rho[1] * rhoj;
  Bj[2] = pj->mhd_data.B_over_rho[2] * rhoj;

  const double B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const double B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  double dB[3];
  dB[0] = Bi[0] - Bj[0];
  dB[1] = Bi[1] - Bj[1];
  dB[2] = Bi[2] - Bj[2];

  const double dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  const double psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const double psi_over_ch_j = pj->mhd_data.psi_over_ch;

  const double permeability_inv = 1.0 / mu_0;

  /* Get the kernel for hi. */
  const double hi_inv = 1.0 / hi;
  const double hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const double xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const double wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const double hj_inv = 1.0 / hj;
  const double hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const double xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const double wj_dr = hjd_inv * wj_dx;

  /* Variable smoothing length term */
  const double f_ij = 1.0 - pi->force.f / mj;
  const double f_ji = 1.0 - pj->force.f / mi;

  /* B dot r. */
  const double Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const double Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];
  const double dBdr = Bri - Brj;

  /* Compute gradient terms */
  const double over_rho2_i = 1.0 / (rhoi * rhoi) * f_ij;
  const double over_rho2_j = 1.0 / (rhoj * rhoj) * f_ji;

  /* Compute symmetric div(B) */
  const double grad_term_i = over_rho2_i * wi_dr * r_inv;
  const double grad_term_j = over_rho2_j * wj_dr * r_inv;

  const double asym_grad_term_i = f_ij * wi_dr * r_inv / rhoi;

  const double divB_i = dBdr * asym_grad_term_i;

  pi->mhd_data.divB -= mj * divB_i;

  /* SPH acceleration term in x direction, i_th particle */
  double sph_acc_term_i[3] = {0.0, 0.0, 0.0};

  /* Accelerations along X */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[0] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[0];
  sph_acc_term_i[0] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[0];

  /* Anisotropic MHD term */
  sph_acc_term_i[0] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[0];
  sph_acc_term_i[0] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[0];

  /* Accelerations along Y */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[1] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[1];
  sph_acc_term_i[1] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[1];

  /* Anisotropic MHD term */
  sph_acc_term_i[1] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[1];
  sph_acc_term_i[1] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[1];

  /* Accelerations along Z */

  /* Isotropic MHD pressure term */
  sph_acc_term_i[2] +=
      0.5 * B2i * permeability_inv * over_rho2_i * wi_dr * r_inv * dx[2];
  sph_acc_term_i[2] +=
      0.5 * B2j * permeability_inv * over_rho2_j * wj_dr * r_inv * dx[2];

  /* Anisotropic MHD term */
  sph_acc_term_i[2] +=
      - over_rho2_i * wi_dr * Bri * permeability_inv * r_inv * Bi[2];
  sph_acc_term_i[2] +=
      - over_rho2_j * wj_dr * Brj * permeability_inv * r_inv * Bj[2];

  /* Divergence cleaning term */
  /* Manifestly *NOT* symmetric in i <-> j */

  const double monopole_beta = pi->mhd_data.monopole_beta;

  const double plasma_beta_i = B2i != 0.0 ? 2.0 * mu_0 * Pi / B2i : 10.0;
  const double scale_i = 0.125 * (10.0 - plasma_beta_i);
  const double tensile_correction_scale_i = fmax(0.0, fmin(scale_i, 1.0));

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
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.tot_mag_F[k] -= mj * sph_acc_term_i[k];
  }

  /* */
  const double dB_dt_pref_i = over_rho2_i * wi_dr * r_inv;

  /* */
  double dB_dt_i[3];
  dB_dt_i[0] = -Bri * dv[0];
  dB_dt_i[1] = -Bri * dv[1];
  dB_dt_i[2] = -Bri * dv[2];

  /* */
  pi->mhd_data.B_over_rho_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  /* Physical resistivity */
  const double resistive_eta_i = pi->mhd_data.resistive_eta;

  const double rho_term_PR = 1.0 / (rhoi * rhoj);
  const double grad_term_PR = f_ij * wi_dr + f_ji * wj_dr;

  const double dB_dt_pref_PR = rho_term_PR * grad_term_PR * r_inv;

  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_over_rho_dt[k] +=
        resistive_eta_i * mj * dB_dt_pref_PR * dB[k];
  }

  pi->u_dt -=
      0.5 * permeability_inv * resistive_eta_i * mj * dB_dt_pref_PR * dB_2;

  /* Artificial resistivity */
  const double alpha_AR = 0.5 * (pi->mhd_data.alpha_AR + pj->mhd_data.alpha_AR);

  const double vsig_AR =
      0.5 * (pi->mhd_data.Alfven_speed + pj->mhd_data.Alfven_speed);

  const double rhoij = 0.5 * (rhoi + rhoj);
  const double rhoij2 = rhoij * rhoij;
  const double rhoij2_inv = 1.0 / rhoij2;

  const double grad_term = 0.5 * (f_ij * wi_dr + f_ji * wj_dr);

  const double art_diff_pref = alpha_AR * vsig_AR * rhoij2_inv * grad_term;

  pi->mhd_data.B_over_rho_dt[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt[2] += mj * art_diff_pref * dB[2];

  pi->u_dt -= 0.5 * mj * permeability_inv * art_diff_pref * dB_2;

  /* Store AR terms */
  pi->mhd_data.B_over_rho_dt_AR[0] += mj * art_diff_pref * dB[0];
  pi->mhd_data.B_over_rho_dt_AR[1] += mj * art_diff_pref * dB[1];
  pi->mhd_data.B_over_rho_dt_AR[2] += mj * art_diff_pref * dB[2];

  pi->mhd_data.u_dt_AR -= 0.5 * mj * permeability_inv * art_diff_pref * dB_2;

  /* Divergence diffusion */
  const double vsig_Dedner_i = 0.5 * pi->viscosity.v_sig;
  const double vsig_Dedner_j = 0.5 * pj->viscosity.v_sig;

  double grad_psi = grad_term_i * psi_over_ch_i * vsig_Dedner_i;
  grad_psi += grad_term_j * psi_over_ch_j * vsig_Dedner_j;

  pi->mhd_data.B_over_rho_dt[0] -= mj * grad_psi * dx[0];
  pi->mhd_data.B_over_rho_dt[1] -= mj * grad_psi * dx[1];
  pi->mhd_data.B_over_rho_dt[2] -= mj * grad_psi * dx[2];

  /* Save induction sources */
  const double dB_dt_pref_Lap = 2.0 * r_inv / rhoj;

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
