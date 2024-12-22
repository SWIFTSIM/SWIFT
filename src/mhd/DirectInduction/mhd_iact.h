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
                               const float mu_0, const float a,
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
  for (int k = 0; k < 3; k++) {
    Bi[k] = pi->mhd_data.B[k];
    Bj[k] = pj->mhd_data.B[k];
  }

  float dB[3];
  for (int k = 0; k < 3; k++) dB[k] = Bi[k] - Bj[k];

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
  const float grad_term_i = f_ij * wi_dr * r_inv / rhoi;
  const float grad_term_j = f_ji * wj_dr * r_inv / rhoj;

  /* Calculate divergence */
  float div_B_i = grad_term_i * (Bri - Brj);
  float div_B_j = grad_term_j * (Bri - Brj);
  pi->mhd_data.div_B -= mj * div_B_i;
  pj->mhd_data.div_B -= mi * div_B_j;

  /* Calculate curl */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.curl_B[k] += mj * grad_term_i * dB_cross_dx[k];
    pj->mhd_data.curl_B[k] += mi * grad_term_j * dB_cross_dx[k];
  }
  
  /* Calculate gradient of B tensor */
  for (int k = 0; k < 3; k++) {
    for (int m = 0; m < 3; m++) {
      pi->mhd_data.grad_B_tensor[k][m] -=
          mj * grad_term_i * dB[k] * dx[m];
      pj->mhd_data.grad_B_tensor[k][m] -=
          mi * grad_term_j * dB[k] * dx[m];
    }
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
  float wi, wi_dx;
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  
  float Bi[3], Bj[3];
  for (int k = 0; k < 3; k++) {
    Bi[k] = pi->mhd_data.B[k];
    Bj[k] = pj->mhd_data.B[k];
  }

  float dB[3];
  for (int k = 0; k < 3; k++) dB[k] = Bi[k] - Bj[k];

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;

  /* B dot r. */
  const float Bri = Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2];
  const float Brj = Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2];

  /* dB cross r */
  float dB_cross_dx[3];
  dB_cross_dx[0] = dB[1] * dx[2] - dB[2] * dx[1];
  dB_cross_dx[1] = dB[2] * dx[0] - dB[0] * dx[2];
  dB_cross_dx[2] = dB[0] * dx[1] - dB[1] * dx[0];

  /* Compute gradient terms */
  const float grad_term_i = f_ij * wi_dr * r_inv / rhoi;

  /* Calculate divergence */
  float div_B_i = grad_term_i * (Bri - Brj);
  pi->mhd_data.div_B -= mj * div_B_i;

  /* Calculate curl */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.curl_B[k] += mj * grad_term_i * dB_cross_dx[k];
  }
  
  /* Calculate gradient of B tensor */
  for (int k = 0; k < 3; k++) {
    for (int m = 0; m < 3; m++) {
      pi->mhd_data.grad_B_tensor[k][m] -=
          mj * grad_term_i * dB[k] * dx[m];
    }
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

  /* Force coupling between MHD and the rest */
  const float permeability_inv = 1.0f / mu_0;
  
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

  float Bi[3], Bj[3];
  for (int k = 0; k < 3; k++) {
    Bi[k] = pi->mhd_data.B[k];
    Bj[k] = pj->mhd_data.B[k];
  }

  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;
  
  /* Useful combinations of magnetic field components */ 
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  float dB[3];
  for (int k = 0; k < 3; k++) dB[k] = Bi[k] - Bj[k];
  
  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  /* Useful terms involving position and velocity differences */ 
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  
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
  const float grad_term_i = f_ij * wi_dr * r_inv / (rhoi * rhoi);
  const float grad_term_j = f_ji * wj_dr * r_inv / (rhoj * rhoj);

  /* SPH acceleration term, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  for (int k = 0; k < 3; k++) {

    /* Isotropic MHD pressure term */
    sph_acc_term_i[k] += 0.5f * B2i * grad_term_i * dx[k];
    sph_acc_term_i[k] += 0.5f * B2j * grad_term_j * dx[k];

    /* Anisotropic MHD tension term */
    sph_acc_term_i[k] -= grad_term_i * Bri * Bi[k];
    sph_acc_term_i[k] -= grad_term_j * Brj * Bj[k];

  }

  /* SPH acceleration term, j_th particle */
  float sph_acc_term_j[3];
  for (int k = 0; k < 3; k++) {
    sph_acc_term_j[k] = - sph_acc_term_i[k];
  }
  
  /* Tensile instability correction term */
  const float monopole_beta_i = pi->mhd_data.monopole_beta;
  const float monopole_beta_j = pj->mhd_data.monopole_beta;
  
  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  const float plasma_beta_j = B2j != 0.0f ? 2.0f * mu_0 * Pj / B2j : FLT_MAX;

  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  const float scale_j = 0.125f * (10.0f - plasma_beta_j);

  const float tensile_correction_scale_i = monopole_beta_i * fmaxf(0.0f, fminf(scale_i, 1.0f));
  const float tensile_correction_scale_j = monopole_beta_j * fmaxf(0.0f, fminf(scale_j, 1.0f));

  for (int k = 0; k < 3; k++) {

    sph_acc_term_i[k] += tensile_correction_scale_i * Bri * grad_term_i * Bi[k];
    sph_acc_term_i[k] += tensile_correction_scale_i * Brj * grad_term_j * Bi[k];

    sph_acc_term_j[k] -= tensile_correction_scale_j * Bri * grad_term_i * Bj[k];
    sph_acc_term_j[k] -= tensile_correction_scale_j * Brj * grad_term_j	* Bj[k];

  }
  
  /* Scale forces appropriately */
  for (int k = 0; k < 3; k++) {
    sph_acc_term_i[k] *= permeability_inv;
    sph_acc_term_j[k] *= permeability_inv;
  }
  
  /* Use the force Luke ! */
  for (int k = 0; k < 3; k++) {
    pi->a_hydro[k] -= mj * sph_acc_term_i[k];
    pj->a_hydro[k] -= mi * sph_acc_term_j[k];
  }
  
  /* Physical induction */
  const float induction_grad_term_i = f_ij * wi_dr * r_inv / rhoi;
  const float induction_grad_term_j = f_ji * wj_dr * r_inv / rhoj;
  
  float dB_dt_i[3], dB_dt_j[3];
  for (int k = 0; k < 3; k++) {
    dB_dt_i[k] = induction_grad_term_i * (Bri * dv[k] - dvdr * Bi[k]);
    dB_dt_j[k] = induction_grad_term_j * (Brj * dv[k] - dvdr * Bj[k]);
  }

  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] -= mj * dB_dt_i[k];
    pj->mhd_data.B_dt[k] -= mi * dB_dt_j[k];
  }

  /* Physical resistivity */
  const float resistive_eta_i = pi->mhd_data.resistive_eta;
  const float resistive_eta_j = pj->mhd_data.resistive_eta;

  const float dB_dt_PR_prefactor_i = 2.0f * f_ij * wi_dr * r_inv / rhoj;
  const float dB_dt_PR_prefactor_j = 2.0f * f_ji * wj_dr * r_inv / rhoi;
  
  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] += resistive_eta_i * mj * dB_dt_PR_prefactor_i * dB[k];
    pj->mhd_data.B_dt[k] -= resistive_eta_j * mi * dB_dt_PR_prefactor_j * dB[k];
  }
  
  /* Artificial resistivity */
  const float art_diff_beta_i = pi->mhd_data.art_diff_beta;
  const float art_diff_beta_j = pj->mhd_data.art_diff_beta;

  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_i_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                            dv_cross_dx[1] * dv_cross_dx[1] +
                            dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B_i = sqrtf(v_sig_B_i_2) * r_inv;
  const float v_sig_B_j = v_sig_B_i;
  
  const float dB_dt_AR_grad_term_i = f_ij * wi_dr / (rhoi * rhoi);
  const float dB_dt_AR_grad_term_j = f_ji * wj_dr / (rhoj * rhoj);

  float dB_dt_AR_pref_i = art_diff_beta_i * v_sig_B_i * dB_dt_AR_grad_term_i;
  dB_dt_AR_pref_i += art_diff_beta_j * v_sig_B_j * dB_dt_AR_grad_term_j;
  const float dB_dt_AR_pref_j = dB_dt_AR_pref_i;

  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] += 0.5f * rhoi * mj * dB_dt_AR_pref_i * dB[k];
    pj->mhd_data.B_dt[k] -= 0.5f * rhoj * mi * dB_dt_AR_pref_j * dB[k];
  }

  /* Put back diffused magnetic energy as heat */ 
  pi->u_dt -= 0.25f * permeability_inv * mj * dB_dt_AR_pref_i * dB_2;
  pj->u_dt -= 0.25f * permeability_inv * mi * dB_dt_AR_pref_j * dB_2;

  /* Divergence cleaning */
  const float ch_i = mhd_get_magnetosonic_speed(pi, mu_0);
  const float ch_j = mhd_get_magnetosonic_speed(pj, mu_0);

  float dB_dt_Dedner_pref_i = psi_over_ch_i * ch_i * grad_term_i;
  dB_dt_Dedner_pref_i += psi_over_ch_j * ch_j * grad_term_j;
  float dB_dt_Dedner_pref_j = dB_dt_Dedner_pref_i;

  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] -= rhoi * mj * dB_dt_Dedner_pref_i * dx[k];
    pj->mhd_data.B_dt[k] += rhoj * mi * dB_dt_Dedner_pref_j * dx[k];
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

  /* Force coupling between MHD and the rest */  
  const float permeability_inv = 1.0f / mu_0;
  
  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float Pi = pi->force.pressure;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  float Bi[3], Bj[3];
  for (int k = 0; k < 3; k++) {
    Bi[k] = pi->mhd_data.B[k];
    Bj[k] = pj->mhd_data.B[k];
  }
 
  const float psi_over_ch_i = pi->mhd_data.psi_over_ch;
  const float psi_over_ch_j = pj->mhd_data.psi_over_ch;

  /* Useful combinations of magnetic field components */
  const float B2i = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float B2j = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];
  
  float dB[3];
  for (int k = 0; k < 3; k++) dB[k] = Bi[k] - Bj[k];
  
  const float dB_2 = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  /* Useful terms involving position and velocity differences */
  float dv[3];
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  
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
  const float grad_term_i = f_ij * wi_dr * r_inv / (rhoi * rhoi);
  const float grad_term_j = f_ji * wj_dr * r_inv / (rhoj * rhoj);
  
  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3] = {0.f, 0.f, 0.f};

  for (int k = 0; k < 3; k++) {

    /* Isotropic MHD pressure term */
    sph_acc_term_i[k] += 0.5f * B2i * grad_term_i * dx[k];
    sph_acc_term_i[k] += 0.5f * B2j * grad_term_j * dx[k];

    /* Anisotropic MHD tension term */
    sph_acc_term_i[k] -= grad_term_i * Bri * Bi[k];
    sph_acc_term_i[k] -= grad_term_j * Brj * Bj[k];

  }
  
  /* Tensile instability corection term */
  const float monopole_beta_i = pi->mhd_data.monopole_beta;
  
  const float plasma_beta_i = B2i != 0.0f ? 2.0f * mu_0 * Pi / B2i : FLT_MAX;
  
  const float scale_i = 0.125f * (10.0f - plasma_beta_i);
  
  const float tensile_correction_scale_i = monopole_beta_i * fmaxf(0.0f, fminf(scale_i, 1.0f));

  for (int k = 0; k < 3; k++) {

    sph_acc_term_i[k] += tensile_correction_scale_i * Bri * grad_term_i * Bi[k];
    sph_acc_term_i[k] += tensile_correction_scale_i * Brj * grad_term_j	* Bi[k];

  }
  
  /* Scale forces appropriately */
  for (int k = 0; k < 3; k++) {
    sph_acc_term_i[k] *= permeability_inv;
  }
  
  /* Use the force Luke ! */
  for (int k = 0; k < 3; k++) {
    pi->a_hydro[k] -= mj * sph_acc_term_i[k];
  }
  
  /* Physical induction */
  const float induction_grad_term_i = f_ij * wi_dr * r_inv / rhoi;

  float dB_dt_i[3];
  for (int k = 0; k < 3; k++) {
    dB_dt_i[k] = induction_grad_term_i * (Bri * dv[k] - dvdr * Bi[k]);
  }

  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] -= mj * dB_dt_i[k];
  }
  
  /* Physical resistivity */
  const float resistive_eta_i = pi->mhd_data.resistive_eta;

  const float dB_dt_PR_prefactor_i = 2.0f * f_ij * wi_dr * r_inv / rhoj;
  
  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] += resistive_eta_i * mj * dB_dt_PR_prefactor_i * dB[k];
  }
  
  /* Artificial resistivity */
  const float art_diff_beta_i = pi->mhd_data.art_diff_beta;
  const float art_diff_beta_j = pj->mhd_data.art_diff_beta;
  
  float dv_cross_dx[3];
  dv_cross_dx[0] = dv[1] * dx[2] - dv[2] * dx[1];
  dv_cross_dx[1] = dv[2] * dx[0] - dv[0] * dx[2];
  dv_cross_dx[2] = dv[0] * dx[1] - dv[1] * dx[0];

  const float v_sig_B_i_2 = dv_cross_dx[0] * dv_cross_dx[0] +
                            dv_cross_dx[1] * dv_cross_dx[1] +
                            dv_cross_dx[2] * dv_cross_dx[2];
  const float v_sig_B_i = sqrtf(v_sig_B_i_2) * r_inv;
  const float v_sig_B_j = v_sig_B_i;
  
  const float dB_dt_AR_grad_term_i = f_ij * wi_dr / (rhoi * rhoi);
  const float dB_dt_AR_grad_term_j = f_ji * wj_dr / (rhoj * rhoj);

  float dB_dt_AR_pref_i = art_diff_beta_i * v_sig_B_i * dB_dt_AR_grad_term_i;
  dB_dt_AR_pref_i += art_diff_beta_j * v_sig_B_j * dB_dt_AR_grad_term_j;
  
  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] += 0.5f * rhoi * mj * dB_dt_AR_pref_i * dB[k];
  }

  /* Put back diffused magnetic energy as heat */ 
  pi->u_dt -= 0.25f * permeability_inv * mj * dB_dt_AR_pref_i * dB_2;
  
  /* Divergence cleaning */
  const float ch_i = mhd_get_magnetosonic_speed(pi, mu_0);
  const float ch_j = mhd_get_magnetosonic_speed(pj, mu_0);

  float dB_dt_Dedner_pref_i = psi_over_ch_i * ch_i * grad_term_i;
  dB_dt_Dedner_pref_i += psi_over_ch_j * ch_j * grad_term_j;

  /* Update time derivative of B */
  for (int k = 0; k < 3; k++) {
    pi->mhd_data.B_dt[k] -= rhoi * mj * dB_dt_Dedner_pref_i * dx[k];
  }
}

#endif /* SWIFT_DIRECT_INDUCTION_MHD_H */
