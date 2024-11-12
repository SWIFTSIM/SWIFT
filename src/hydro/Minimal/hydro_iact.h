/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MINIMAL_HYDRO_IACT_H
#define SWIFT_MINIMAL_HYDRO_IACT_H

/**
 * @file Minimal/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adaptive_softening_iact.h"
#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"
#include "signal_velocity.h"

#include <float.h>

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  float wi, wi_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
  adaptive_softening_add_correction_term(pi, ui, h_inv, mj);

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];
}

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  runner_iact_nonsym_density(r2, dx, hi, hj, pi, pj, mu_0, a, H);
  const float dx_inv[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_density(r2, dx_inv, hj, hi, pj, pi, mu_0, a, H);

#if 0
  
  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

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

  /* Compute dv dot r */
  float dv[3], curlvr[3];

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->density.div_v -= faci * dvdr;
  pj->density.div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += faci * curlvr[0];
  pi->density.rot_v[1] += faci * curlvr[1];
  pi->density.rot_v[2] += faci * curlvr[2];

  pj->density.rot_v[0] += facj * curlvr[0];
  pj->density.rot_v[1] += facj * curlvr[1];
  pj->density.rot_v[2] += facj * curlvr[2];

#endif
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * Nothing to do here in this scheme.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  /* MATTHIEU START --------------------------------------- */

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get masses and densities */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get B for both particles */
  const float Bi[3] = {pi->B_over_rho[0] * rhoi,   // x
                       pi->B_over_rho[1] * rhoi,   // y
                       pi->B_over_rho[2] * rhoi};  // z
  const float Bj[3] = {pj->B_over_rho[0] * rhoj,   // x
                       pj->B_over_rho[1] * rhoj,   // y
                       pj->B_over_rho[2] * rhoj};  // z

  /* Get the kernel for hi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Difference in B between particles */
  const float dB[3] = {Bi[0] - Bj[0],   // x
                       Bi[1] - Bj[1],   // y
                       Bi[2] - Bj[2]};  // z

  /* Compute dB dot r */
  const float faci = mj * wi_dx * r_inv;
  const float dBdr = dB[0] * dx[0] + dB[1] * dx[1] + dB[2] * dx[2];

  pi->div_B -= faci * dBdr;

  /* Compute dB cross r */
  const float curlBr[3] = {dB[1] * dx[2] - dB[2] * dx[1],   // x
                           dB[2] * dx[0] - dB[0] * dx[2],   // y
                           dB[0] * dx[1] - dB[1] * dx[0]};  // z

  pi->curl_B[0] += faci * curlBr[0];
  pi->curl_B[1] += faci * curlBr[1];
  pi->curl_B[2] += faci * curlBr[2];

  /* MATTHIEU END ----------------------------------------- */
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * Nothing to do here in this scheme.
 *
 * @param r2 Comoving squared distance between particle i and particle j.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  runner_iact_nonsym_gradient(r2, dx, hi, hj, pi, pj, mu_0, a, H);
  const float dx_inv[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_gradient(r2, dx_inv, hj, hi, pj, pi, mu_0, a, H);
}

/**
 * @brief Force interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float mu_0,
    const float a, const float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute signal velocity */
  const float v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Adaptive softening acceleration term */
  const float adapt_soft_acc_term =
      adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv);

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr * f_ij;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

  /* MATTHIEU START --------------------------------------- */

  /* Thermal diffusion ------------------------------------ */

  /* Diffusion signal velocity
   * Price 2018, eq. 43 */
  const float v_diff =
      sqrtf(2.0 * fabsf(pressurei - pressurej) / (rhoi + rhoj));
  // TODO: cosmo terms

  /* Take average of both diffusion constants */
  const float alpha_diff = 0.5f * (pi->alpha_u + pj->alpha_u);

  /* Price 2018, eq. 42 second term
   * Note that wi_dx + wj_dx / 2 is F_ij */
  const float diff_du_term = alpha_diff * v_diff * (pi->u - pj->u) * 0.5f *
                             (wi_dr * f_ij / rhoi + wj_dr * f_ji / rhoj);

  pi->u_dt += mj * diff_du_term;

  /* MHD equations ---------------------------------------- */

  const float one_over_mu0 = 1.f / mu_0;

  /* Get B for both particles */
  const float Bi[3] = {pi->B_over_rho[0] * rhoi,   // x
                       pi->B_over_rho[1] * rhoi,   // y
                       pi->B_over_rho[2] * rhoi};  // z
  const float Bj[3] = {pj->B_over_rho[0] * rhoj,   // x
                       pj->B_over_rho[1] * rhoj,   // y
                       pj->B_over_rho[2] * rhoj};  // z

  /* Difference in B between particles */
  const float dB[3] = {Bi[0] - Bj[0],   // x
                       Bi[1] - Bj[1],   // y
                       Bi[2] - Bj[2]};  // z

  /* Square of the norm of the difference in B between particles */
  const float square_norm_dB = dB[0] * dB[0] + dB[1] * dB[1] + dB[2] * dB[2];

  /* Square norm of B fields */
  const float Bi_2 = Bi[0] * Bi[0] + Bi[1] * Bi[1] + Bi[2] * Bi[2];
  const float Bj_2 = Bj[0] * Bj[0] + Bj[1] * Bj[1] + Bj[2] * Bj[2];

  /* Velocity difference */
  const float dv[3] = {pi->v[0] - pj->v[0],   // x
                       pi->v[1] - pj->v[1],   // y
                       pi->v[2] - pj->v[2]};  // z

  const float weight_i = f_ij * wi_dr * r_inv / (rhoi * rhoi);
  const float weight_j = f_ji * wj_dr * r_inv / (rhoj * rhoj);

  /* Induction equation ---------------------------------- */

  /* Induction equation (Price 2012, eq. 110) */
  float dB_over_rho_dt_i[3] = {0.f, 0.f, 0.f};
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      dB_over_rho_dt_i[k] += dv[k] * weight_i * Bi[l] * dx[l];
    }
  }

  pi->B_over_rho_dt[0] -= mj * dB_over_rho_dt_i[0];
  pi->B_over_rho_dt[1] -= mj * dB_over_rho_dt_i[1];
  pi->B_over_rho_dt[2] -= mj * dB_over_rho_dt_i[2];

  /* Accelerations  -------------------------------------- */

  /* Stress tensor (Price 2012, eq. 116) not including pressure
   * as we did pressure above already. */
  float S_i[3][3];
  float S_j[3][3];
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      S_i[k][l] = -Bi[k] * Bi[l];
      S_j[k][l] = -Bj[k] * Bj[l];
    }
  }
  for (int k = 0; k < 3; ++k) {
    S_i[k][k] += 0.5f * Bi_2;
    S_j[k][k] += 0.5f * Bj_2;
  }
  for (int k = 0; k < 3; ++k) {
    for (int l = 0; l < 3; ++l) {
      S_i[k][l] *= one_over_mu0;
      S_j[k][l] *= one_over_mu0;
    }
  }

  /* (Price 2012, eq. 115) not including P in S */
  float mhd_acc_i[3] = {0.f, 0.f, 0.f};
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      mhd_acc_i[k] += (S_i[k][l] * weight_i + S_j[k][l] * weight_j) * dx[l];
    }
  }

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * mhd_acc_i[0];
  pi->a_hydro[1] -= mj * mhd_acc_i[1];
  pi->a_hydro[2] -= mj * mhd_acc_i[2];

  /* Tensile instability correction -------------------------- */

  /* Plasma beta */
  const float beta_i = Bi_2 > 0.f ? 2.f * mu_0 * pressurei / Bi_2 : FLT_MAX;

  /* Magnitude of correction (Price 2018, eq. 178) */
  float B_hat_corr_i;
  if (beta_i < 2.f)
    B_hat_corr_i = 1.f;
  else if (beta_i < 10.f)
    B_hat_corr_i = (10.f - beta_i) * 0.125f;
  else
    B_hat_corr_i = 0.f;

  const float B_hat_i[3] = {B_hat_corr_i * Bi[0],   // x
                            B_hat_corr_i * Bi[1],   // y
                            B_hat_corr_i * Bi[2]};  // z

  /* Tensile correction (Price 2012, eq. 129 second term)
   * (Note that a 1/mu0 factor is missing in the paper) */
  float mhd_acc_corr_i[3] = {0.f, 0.f, 0.f};
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      mhd_acc_corr_i[k] +=
          B_hat_i[k] * (Bi[l] * weight_i + Bj[l] * weight_j) * dx[l];
    }
  }

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * one_over_mu0 * mhd_acc_corr_i[0];
  pi->a_hydro[1] -= mj * one_over_mu0 * mhd_acc_corr_i[1];
  pi->a_hydro[2] -= mj * one_over_mu0 * mhd_acc_corr_i[2];

  /* Dedner correction ----------------------------------- */

  /* Magnetosonic speed (Price 2018, eq. 180) */
  const float c_h_i = sqrtf(pi->force.soundspeed * pi->force.soundspeed +
                            pi->Alfven_speed * pi->Alfven_speed);
  const float c_h_j = sqrtf(pj->force.soundspeed * pj->force.soundspeed +
                            pj->Alfven_speed * pj->Alfven_speed);

  const float Psi_i = pi->Dedner_Psi_over_c * c_h_i;
  const float Psi_j = pj->Dedner_Psi_over_c * c_h_j;

  /* Induction equation (Price 2018, eq. 172 second term) */
  float dB_over_rho_dt_Dedner_i[3] = {0.f, 0.f, 0.f};
  for (int k = 0; k < 3; k++) {
    dB_over_rho_dt_Dedner_i[k] += (weight_i * Psi_i + weight_j * Psi_j) * dx[k];
  }

  pi->B_over_rho_dt[0] -= mj * dB_over_rho_dt_Dedner_i[0];
  pi->B_over_rho_dt[1] -= mj * dB_over_rho_dt_Dedner_i[1];
  pi->B_over_rho_dt[2] -= mj * dB_over_rho_dt_Dedner_i[2];

  /* Evolution of Dedner scalar ------------------------- */

  /* Div B for Dedner (Price 2018, eq. 172 first term)
   * (multiplication by c_h takes place later)
   * Note we have a negative sign overall as we want div B here
   * for checks and not -div B as the paper assumes */
  float Dedner_div_B = 0.f;
  for (int l = 0; l < 3; l++) {
    Dedner_div_B += weight_i * rhoi * dB[l] * dx[l];
  }

  /* Div B for Dedner (Price 2018, eq. 172 second term)
   * (multiplication by Psi/2c_h takes place later)
   * Note we have a negative sign overall as we want div B here
   * for checks and not -div v as the paper assumes */
  float Dedner_div_v = 0.f;
  for (int l = 0; l < 3; l++) {
    Dedner_div_v += weight_i * rhoi * dv[l] * dx[l];
  }

  pi->Dedner_div_B -= mj * Dedner_div_B;
  pi->Dedner_div_v -= mj * Dedner_div_v;

  /* Shock capturing (artifical resistivity) ------------ */

  /* Artificial resistivity speed (Price 2018, eq. 184) */
  const float alpha_B_i = 1.f;
  const float alpha_B_j = 1.f;

  const float v_sig_AR_i = c_h_i * alpha_B_i;
  const float v_sig_AR_j = c_h_j * alpha_B_j;

  /* Artifical resistivity (Price 2018, eq. 181 second term)
   * without the rho_a in front
   * Note: wi_dx + wj_dx / 2 is F_ij */
  float dB_over_rho_dt_AR_i[3] = {0.f, 0.f, 0.f};
  for (int k = 0; k < 3; k++) {
    dB_over_rho_dt_AR_i[k] += v_sig_AR_i * wi_dr * f_ij / (rhoi * rhoi);
    dB_over_rho_dt_AR_i[k] += v_sig_AR_j * wj_dr * f_ji / (rhoj * rhoj);
    dB_over_rho_dt_AR_i[k] *= 0.5f * dB[k];
  }

  pi->B_over_rho_dt[0] += mj * dB_over_rho_dt_AR_i[0];
  pi->B_over_rho_dt[1] += mj * dB_over_rho_dt_AR_i[1];
  pi->B_over_rho_dt[2] += mj * dB_over_rho_dt_AR_i[2];

  /* Artificial resistivity energy change (Price 2018, eq. 182) */
  float du_dt_AR_i = 0.f;
  for (int k = 0; k < 3; k++) {
    dB_over_rho_dt_AR_i[k] += v_sig_AR_i * wi_dr * f_ij / (rhoi * rhoi);
    dB_over_rho_dt_AR_i[k] += v_sig_AR_j * wj_dr * f_ji / (rhoj * rhoj);
    dB_over_rho_dt_AR_i[k] *= -0.25f * square_norm_dB;
  }

  pi->u_dt -= mj * du_dt_AR_i;

  /* MATTHIEU END --------------------------------------- */
}

/**
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param mu_0 Vaccum permeability in internal units (for the v_sig in the MHD
 * case).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float mu_0,
    const float a, const float H) {

  runner_iact_nonsym_force(r2, dx, hi, hj, pi, pj, mu_0, a, H);
  const float dx_inv[3] = {-dx[0], -dx[1], -dx[2]};
  runner_iact_nonsym_force(r2, dx_inv, hj, hi, pj, pi, mu_0, a, H);

#if 0

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute signal velocity */
  const float v_sig =
      signal_velocity(dx, pi, pj, mu_ij, const_viscosity_beta, a, mu_0);

  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Adaptive softening acceleration term */
  const float adapt_soft_acc_term =
      adaptive_softening_get_acc_term(pi, pj, wi_dr, wj_dr, f_ij, f_ji, r_inv);

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term + adapt_soft_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = P_over_rho2_j * dvdr * r_inv * wj_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr * f_ij;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr * f_ji;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

#endif
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
