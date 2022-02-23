/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include "adiabatic_index.h"
#include "hydro_parameters.h"
#include "minmax.h"

/**
 * @brief Density interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

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

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

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
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

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
 * @brief Force interaction between two particles.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

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
  float Bi[3];
  float Bj[3];
  Bi[0] = pi->B[0];
  Bi[1] = pi->B[1];
  Bi[2] = pi->B[2];
  Bj[0] = pj->B[0];
  Bj[1] = pj->B[1];
  Bj[2] = pj->B[2];

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
  const float isoPi = pressurei + 0.5f * B2i / const_vacuum_permeability;
  const float isoPj = pressurej + 0.5f * B2j / const_vacuum_permeability;

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]) /
                    const_vacuum_permeability;
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]) /
                    const_vacuum_permeability;

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float c2i = ci * ci;
  const float c2j = cj * cj;
  const float v_A2i = B2i / (rhoi * const_vacuum_permeability);
  const float v_A2j = B2j / (rhoj * const_vacuum_permeability);
  const float c2effi = c2i + v_A2i;
  const float c2effj = c2j + v_A2j;
  const float v_sig2i =
      0.5f * (c2effi + sqrtf(c2effi * c2effi -
                             4.0f * c2i * (Bri * r_inv) * (Bri * r_inv) *
                                 const_vacuum_permeability / rhoi));
  const float v_sig2j =
      0.5f * (c2effj + sqrtf(c2effj * c2effj -
                             4.0f * c2j * (Brj * r_inv) * (Brj * r_inv) *
                                 const_vacuum_permeability / rhoj));
  const float v_sig =
      sqrtf(v_sig2i) + sqrtf(v_sig2j) - const_viscosity_beta * mu_ij;

  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* SPH acceleration term in x direction, i_th particle */
  float sph_acc_term_i[3];
  sph_acc_term_i[0] =
      (over_rho2_i * wi_dr * (isoPi * dx[0] - Bri * Bi[0] + Bri * Bi[0]) +
       over_rho2_j * wj_dr * (isoPj * dx[0] - Brj * (Bj[0] - Bi[0]))) *
      r_inv;

  /* SPH acceleration term in y direction */
  sph_acc_term_i[1] =
      (over_rho2_i * wi_dr * (isoPi * dx[1] - Bri * Bi[1] + Bri * Bi[1]) +
       over_rho2_j * wj_dr * (isoPj * dx[1] - Brj * (Bj[1] - Bi[1]))) *
      r_inv;

  /* SPH acceleration term in z direction */
  sph_acc_term_i[2] =
      (over_rho2_i * wi_dr * (isoPi * dx[2] - Bri * Bi[2] + Bri * Bi[2]) +
       over_rho2_j * wj_dr * (isoPj * dx[2] - Brj * (Bj[2] - Bi[2]))) *
      r_inv;

  /* SPH acceleration term in x direction, j_th particle */
  float sph_acc_term_j[3];
  sph_acc_term_j[0] =
      (over_rho2_j * wj_dr * (-isoPj * dx[0] + Brj * Bj[0] - Brj * Bj[0]) +
       over_rho2_i * wi_dr * (-isoPi * dx[0] + Bri * (Bi[0] - Bj[0]))) *
      r_inv;

  /* SPH acceleration term in y direction */
  sph_acc_term_j[1] =
      (over_rho2_j * wj_dr * (-isoPj * dx[1] + Brj * Bj[1] - Brj * Bj[1]) +
       over_rho2_i * wi_dr * (-isoPi * dx[1] + Bri * (Bi[1] - Bj[1]))) *
      r_inv;

  /* SPH acceleration term in z direction */
  sph_acc_term_j[2] =
      (over_rho2_j * wj_dr * (-isoPj * dx[2] + Brj * Bj[2] - Brj * Bj[2]) +
       over_rho2_i * wi_dr * (-isoPi * dx[2] + Bri * (Bi[2] - Bj[2]))) *
      r_inv;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term_i[0] + mj * visc_acc_term * dx[0];
  pi->a_hydro[1] -= mj * sph_acc_term_i[1] + mj * visc_acc_term * dx[1];
  pi->a_hydro[2] -= mj * sph_acc_term_i[2] + mj * visc_acc_term * dx[2];

  pj->a_hydro[0] -= mi * sph_acc_term_j[0] - mi * visc_acc_term * dx[0];
  pj->a_hydro[1] -= mi * sph_acc_term_j[1] - mi * visc_acc_term * dx[1];
  pj->a_hydro[2] -= mi * sph_acc_term_j[2] - mi * visc_acc_term * dx[2];

  /* Calculate monopole term */
  float B_mon_i =
      (over_rho2_i * wi_dr * Bri + over_rho2_j * wj_dr * Brj) * r_inv;
  float B_mon_j =
      (over_rho2_i * wi_dr * Bri + over_rho2_j * wj_dr * Brj) * r_inv;
  pi->B_mon += mj * B_mon_i;
  pj->B_mon += mi * B_mon_j;

  /* Get the time derivative for u. */
  const float sph_du_term_i = pressurei * over_rho2_i * dvdr * r_inv * wi_dr;
  const float sph_du_term_j = pressurej * over_rho2_j * dvdr * r_inv * wj_dr;

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

  /* Update the density squared estimate */
  pi->rhosq += mj * rhoj * wi;
  pj->rhosq += mi * rhoi * wj;

  /* */
  const float dB_dt_pref_i = over_rho2_i * rhoi * wi_dr * r_inv;
  const float dB_dt_pref_j = over_rho2_j * rhoj * wj_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = Bi[0] * dvdr - Bri * (pi->v[0] - pj->v[0]);
  dB_dt_i[1] = Bi[1] * dvdr - Bri * (pi->v[1] - pj->v[1]);
  dB_dt_i[2] = Bi[2] * dvdr - Bri * (pi->v[2] - pj->v[2]);

  float dB_dt_j[3];
  dB_dt_j[0] = Bj[0] * dvdr - Brj * (pi->v[0] - pj->v[0]);
  dB_dt_j[1] = Bj[1] * dvdr - Brj * (pi->v[1] - pj->v[1]);
  dB_dt_j[2] = Bj[2] * dvdr - Brj * (pi->v[2] - pj->v[2]);

  /* */
  pi->B_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->B_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->B_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];

  pj->B_dt[0] += mi * dB_dt_pref_j * dB_dt_j[0];
  pj->B_dt[1] += mi * dB_dt_pref_j * dB_dt_j[1];
  pj->B_dt[2] += mi * dB_dt_pref_j * dB_dt_j[2];
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
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

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
  float Bi[3];
  float Bj[3];
  Bi[0] = pi->B[0];
  Bi[1] = pi->B[1];
  Bi[2] = pi->B[2];
  Bj[0] = pj->B[0];
  Bj[1] = pj->B[1];
  Bj[2] = pj->B[2];

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
  const float isoPi = pressurei + 0.5f * B2i / const_vacuum_permeability;
  const float isoPj = pressurej + 0.5f * B2j / const_vacuum_permeability;

  /* B dot r. */
  const float Bri = (Bi[0] * dx[0] + Bi[1] * dx[1] + Bi[2] * dx[2]) /
                    const_vacuum_permeability;
  const float Brj = (Bj[0] * dx[0] + Bj[1] * dx[1] + Bj[2] * dx[2]) /
                    const_vacuum_permeability;

  /* Compute gradient terms */
  const float over_rho2_i = 1.0f / (rhoi * rhoi) * f_ij;
  const float over_rho2_j = 1.0f / (rhoj * rhoj) * f_ji;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float c2i = ci * ci;
  const float c2j = cj * cj;
  const float v_A2i = B2i / (rhoi * const_vacuum_permeability);
  const float v_A2j = B2j / (rhoj * const_vacuum_permeability);
  const float c2effi = c2i + v_A2i;
  const float c2effj = c2j + v_A2j;
  const float v_sig2i =
      0.5f * (c2effi + sqrtf(c2effi * c2effi -
                             4.0f * c2i * (Bri * r_inv) * (Bri * r_inv) *
                                 const_vacuum_permeability / rhoi));
  const float v_sig2j =
      0.5f * (c2effj + sqrtf(c2effj * c2effj -
                             4.0f * c2j * (Brj * r_inv) * (Brj * r_inv) *
                                 const_vacuum_permeability / rhoj));
  const float v_sig =
      sqrtf(v_sig2i) + sqrtf(v_sig2j) - const_viscosity_beta * mu_ij;

  /* Grab balsara switches */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * (balsara_i + balsara_j) * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term =
      0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;

  /* SPH acceleration term in x direction */
  float sph_acc_term[3];
  sph_acc_term[0] =
      (over_rho2_i * wi_dr * (isoPi * dx[0] - Bri * Bi[0] + Bri * Bi[0]) +
       over_rho2_j * wj_dr * (isoPj * dx[0] - Brj * (Bj[0] - Bi[0]))) *
      r_inv;

  /* SPH acceleration term in y direction */
  sph_acc_term[1] =
      (over_rho2_i * wi_dr * (isoPi * dx[1] - Bri * Bi[1] + Bri * Bi[1]) +
       over_rho2_j * wj_dr * (isoPj * dx[1] - Brj * (Bj[1] - Bi[1]))) *
      r_inv;

  /* SPH acceleration term in z direction */
  sph_acc_term[2] =
      (over_rho2_i * wi_dr * (isoPi * dx[2] - Bri * Bi[2] + Bri * Bi[2]) +
       over_rho2_j * wj_dr * (isoPj * dx[2] - Brj * (Bj[2] - Bi[2]))) *
      r_inv;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_acc_term[0] + mj * visc_acc_term * dx[0];
  pi->a_hydro[1] -= mj * sph_acc_term[1] + mj * visc_acc_term * dx[1];
  pi->a_hydro[2] -= mj * sph_acc_term[2] + mj * visc_acc_term * dx[2];

  /* Calculate monopole term */
  float B_mon_i =
      (over_rho2_i * wi_dr * Bri + over_rho2_j * wj_dr * Brj) * r_inv;
  pi->B_mon += mj * B_mon_i;

  /* Get the time derivative for u. */
  const float sph_du_term_i = pressurei * over_rho2_i * dvdr * r_inv * wi_dr;

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

  pi->rhosq += mj * rhoj * wi;

  /* */
  const float dB_dt_pref_i = over_rho2_i * rhoi * wi_dr * r_inv;

  /* */
  float dB_dt_i[3];
  dB_dt_i[0] = Bi[0] * dvdr - Bri * (pi->v[0] - pj->v[0]);
  dB_dt_i[1] = Bi[1] * dvdr - Bri * (pi->v[1] - pj->v[1]);
  dB_dt_i[2] = Bi[2] * dvdr - Bri * (pi->v[2] - pj->v[2]);

  /* */
  pi->B_dt[0] += mj * dB_dt_pref_i * dB_dt_i[0];
  pi->B_dt[1] += mj * dB_dt_pref_i * dB_dt_i[1];
  pi->B_dt[2] += mj * dB_dt_pref_i * dB_dt_i[2];
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
