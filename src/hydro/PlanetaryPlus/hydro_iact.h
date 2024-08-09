/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk).
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
#ifndef SWIFT_PLANETARY_HYDRO_IACT_H
#define SWIFT_PLANETARY_HYDRO_IACT_H

/**
 * @file Planetary/hydro_iact.h
 * @brief Minimal conservative implementation of SPH (Neighbour loop equations)
 *
 * The thermal variable is the internal energy (u). Simple constant
 * viscosity term with the Balsara (1995) switch (optional).
 * No thermal conduction term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "const.h"
#include "hydro_density_estimate.h"
#include "hydro_kernels_etc.h"
#include "hydro_misc_utils.h"
#include "hydro_parameters.h"
#include "hydro_viscosity.h"
#include "math.h"
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pj->n_density += wj;
  pi->N_density++;
  pj->N_density++;
#endif

  hydro_runner_iact_density_extra_density_estimate(pi, pj, dx, wi, wj, wi_dx,
                                                   wj_dx);
  hydro_runner_iact_density_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_density_extra_viscosity(pi, pj, dx, wi, wj, wi_dx, wj_dx);
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

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_density += wi;
  pi->N_density++;
#endif

  hydro_runner_iact_nonsym_density_extra_density_estimate(pi, pj, dx, wi,
                                                          wi_dx);
  hydro_runner_iact_nonsym_density_extra_kernel(pi, pj, dx, wi, wi_dx);
  hydro_runner_iact_nonsym_density_extra_viscosity(pi, pj, dx, wi, wi_dx);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j
 *
 * This method wraps around hydro_gradients_collect, which can be an empty
 * method, in which case no gradients are used.
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
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Compute kernel of pi. */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute kernel of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Correction factors for kernel gradients */
  const float rho_inv_i = 1.f / pi->rho;
  const float rho_inv_j = 1.f / pj->rho;

  pi->weighted_wcount += pj->mass * r2 * wi_dx * r_inv;
  pj->weighted_wcount += pi->mass * r2 * wj_dx * r_inv;

  pi->weighted_neighbour_wcount += pj->mass * r2 * wi_dx * rho_inv_j * r_inv;
  pj->weighted_neighbour_wcount += pi->mass * r2 * wj_dx * rho_inv_i * r_inv;

  hydro_runner_iact_gradient_extra_density_estimate(pi, pj, dx, wi, wj, wi_dx,
                                                    wj_dx);
  hydro_runner_iact_gradient_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_gradient_extra_viscosity(pi, pj, dx, wi, wj, wi_dx, wj_dx);
}

/**
 * @brief Calculate the gradient interaction between particle i and particle j:
 * non-symmetric version
 *
 * This method wraps around hydro_gradients_nonsym_collect, which can be an
 * empty method, in which case no gradients are used.
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
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float wi, wj, wi_dx, wj_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Compute kernel of pi. */
  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);
    
      /* Compute kernel of pj. */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Correction factors for kernel gradients */
  // Nominator and denominator for f_gdf factors, Wadsley+2017 Eqn. 8
  const float rho_inv_j = 1.f / pj->rho;
  pi->weighted_wcount += pj->mass * r2 * wi_dx * r_inv;
  pi->weighted_neighbour_wcount += pj->mass * r2 * wi_dx * rho_inv_j * r_inv;

  hydro_runner_iact_nonsym_gradient_extra_density_estimate(pi, pj, dx, wi,
                                                           wi_dx);
  hydro_runner_iact_nonsym_gradient_extra_kernel(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  hydro_runner_iact_nonsym_gradient_extra_viscosity(pi, pj, dx, wi, wi_dx);
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);


  /* Compute the G and Q gradient and viscosity factors, either using matrix
   * inversions or standard GDF SPH */

  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley+2017
  or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);

  /* Density factors for GDF or standard equations */
  float rho_factor_i, rho_factor_j;
  hydro_set_rho_factors(&rho_factor_i, &rho_factor_j, pi, pj);

  /* Calculate the viscous pressures Q */
  float Qi, Qj;
  float visc_signal_velocity, cond_signal_velocity;  
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &cond_signal_velocity, pi, pj, dx, a, H);

  /* set kernel gradient terms to be used in eolution equations */
  float kernel_gradient_i[3], kernel_gradient_j[3];
  float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];
  hydro_set_kernel_gradient_terms(dx, kernel_gradient_i, kernel_gradient_j,
                                  Q_kernel_gradient_i, Q_kernel_gradient_j, Gi,
                                  Gj);
     
  /* Pressure terms to be used in evolution equations */
  float P_i_term = pressurei / rho_factor_i;
  float P_j_term = pressurej / rho_factor_j;
  float Q_i_term = Qi / rho_factor_i;
  float Q_j_term = Qj / rho_factor_j;

  /* Use the force Luke! */
  pi->a_hydro[0] -=
      mj *
      (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] + 
       Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pi->a_hydro[1] -=
      mj *
      (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] +
       Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pi->a_hydro[2] -=
      mj *
      (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] +
       Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);

  pj->a_hydro[0] +=
      mi *
      (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] +
       Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pj->a_hydro[1] +=
      mi *
      (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] +
       Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pj->a_hydro[2] +=
      mi *
      (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] +
       Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);

  /* dx dot kernel gradient term needed for du/dt in e.g. eq 13 of Wadsley and
   * equiv*/
  const float dvdG_i = (pi->v[0] - pj->v[0]) * kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_i[2];
  const float dvdG_j = (pi->v[0] - pj->v[0]) * kernel_gradient_j[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_j[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_j[2];
  const float Q_dvdG_i = (pi->v[0] - pj->v[0]) * Q_kernel_gradient_i[0] +
                         (pi->v[1] - pj->v[1]) * Q_kernel_gradient_i[1] +
                         (pi->v[2] - pj->v[2]) * Q_kernel_gradient_i[2];
  const float Q_dvdG_j = (pi->v[0] - pj->v[0]) * Q_kernel_gradient_j[0] +
                         (pi->v[1] - pj->v[1]) * Q_kernel_gradient_j[1] +
                         (pi->v[2] - pj->v[2]) * Q_kernel_gradient_j[2];

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = P_i_term * dvdG_i + Q_i_term * Q_dvdG_i;
  float du_dt_j = P_j_term * dvdG_j + Q_j_term * Q_dvdG_j;

#ifdef PLANETARY_FIXED_ENTROPY
  du_dt_i = P_i_term * dvdG_i;
  du_dt_j = P_j_term * dvdG_j;
#endif

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdG_i / rhoj;
  pj->force.h_dt -= mi * dvdG_j / rhoi;

    const float v_sig = visc_signal_velocity;
    
  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
    
  pi->drho_dt += mj * (pi->rho / pj->rho)  * dvdG_i;
  pj->drho_dt += mi * (pj->rho / pi->rho)  * dvdG_j;    

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pj->n_force += wi + wj;
  pi->N_force++;
  pj->N_force++;
#endif
 if (!pi->is_h_max && !pj->is_h_max) {   
      if (pi->mat_id == pj->mat_id){  
          float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
          hydro_set_u_rho_cond(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi, pj, dx, a, H); 

          float mean_rho = 0.5f * (pi->rho + pj->rho);  
          float v_sig_cond = cond_signal_velocity;//sqrtf(fabs(pressurei - pressurej) / mean_rho);  //vtilde_signal_velocity;// 0.5f * (ci + cj);//vtilde_signal_velocity;//sqrtf(fabs(pressurei - pressurej) / mean_rho);  
          float mean_G = 0.5f * sqrtf((kernel_gradient_i[0] + kernel_gradient_j[0]) * (kernel_gradient_i[0] + kernel_gradient_j[0]) +
                                 (kernel_gradient_i[1] + kernel_gradient_j[1]) * (kernel_gradient_i[1] + kernel_gradient_j[1]) +
                                 (kernel_gradient_i[2] + kernel_gradient_j[2]) * (kernel_gradient_i[2] + kernel_gradient_j[2]));  

          float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
          
          float alpha_u = mean_balsara;//1.f;//0.1f + mean_balsara * 0.9f;//0.05f;//0.5f;//1.f;//0.1f;//1.f 
          float du_dt_cond_i = -alpha_u * mj * v_sig_cond * (utilde_i - utilde_j) * mean_G / mean_rho;
          float du_dt_cond_j = -alpha_u * mi * v_sig_cond * (utilde_j - utilde_i) * mean_G / mean_rho;   

          pi->u_dt += du_dt_cond_i;
          pj->u_dt += du_dt_cond_j;  

          float alpha_rho = mean_balsara;//1.f;//0.1f + mean_balsara * 0.9f;//0.05f;//0.5f;//1.f;//0.1f;//1.f 
          float drho_dt_cond_i = -alpha_rho * mj * (pi->rho / pj->rho) * v_sig_cond * (rhotilde_i - rhotilde_j) * mean_G / mean_rho; 
          float drho_dt_cond_j = -alpha_rho * mi * (pj->rho / pi->rho) * v_sig_cond * (rhotilde_j - rhotilde_i) * mean_G / mean_rho;   

          pi->drho_dt += drho_dt_cond_i;
          pj->drho_dt += drho_dt_cond_j;
      }
 }
    
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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);

  /* Recover some data */
  const float mj = pj->mass;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);

  /* Compute the G and Q gradient and viscosity factors, either using matrix
   * inversions or standard GDF SPH */

  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley+2017
  or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  hydro_set_Gi_Gj_forceloop(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);

  /* Density factors for GDF or standard equations */
  float rho_factor_i, rho_factor_j;
  hydro_set_rho_factors(&rho_factor_i, &rho_factor_j, pi, pj);

  /* Calculate the viscous pressures Q */
  float Qi, Qj;
  float visc_signal_velocity, cond_signal_velocity;  
  hydro_set_Qi_Qj(&Qi, &Qj, &visc_signal_velocity, &cond_signal_velocity, pi, pj, dx, a, H);

  /* set kernel gradient terms to be used in eolution equations */
  float kernel_gradient_i[3], kernel_gradient_j[3];
  float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];
  hydro_set_kernel_gradient_terms(dx, kernel_gradient_i, kernel_gradient_j,
                                  Q_kernel_gradient_i, Q_kernel_gradient_j, Gi,
                                  Gj);

  /* Pressure terms to be used in evolution equations */
  float P_i_term = pressurei / rho_factor_i;
  float P_j_term = pressurej / rho_factor_j;
  float Q_i_term = Qi / rho_factor_i;
  float Q_j_term = Qj / rho_factor_j;

  /* Use the force Luke! */
  pi->a_hydro[0] -=
      mj *
      (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] +
       Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pi->a_hydro[1] -=
      mj *
      (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] +
       Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pi->a_hydro[2] -=
      mj *
      (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] +
       Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);

  /* dx dot kernel gradient term needed for du/dt in e.g. eq 13 of Wadsley and
   * equiv*/
  const float dvdG_i = (pi->v[0] - pj->v[0]) * kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_i[2];
  const float Q_dvdG_i = (pi->v[0] - pj->v[0]) * Q_kernel_gradient_i[0] +
                         (pi->v[1] - pj->v[1]) * Q_kernel_gradient_i[1] +
                         (pi->v[2] - pj->v[2]) * Q_kernel_gradient_i[2];

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = P_i_term * dvdG_i + Q_i_term * Q_dvdG_i;
#ifdef PLANETARY_FIXED_ENTROPY
  du_dt_i = P_i_term * dvdG_i;
#endif

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdG_i / rhoj;

  const float v_sig = visc_signal_velocity;
    
  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
    
  pi->drho_dt += mj * (pi->rho / pj->rho)  * dvdG_i;  

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
  if (!pi->is_h_max && !pj->is_h_max) {  
      if (pi->mat_id == pj->mat_id){  
          float utilde_i, utilde_j, rhotilde_i, rhotilde_j;
          hydro_set_u_rho_cond(&utilde_i, &utilde_j, &rhotilde_i, &rhotilde_j, pi, pj, dx, a, H); 

          float mean_rho = 0.5f * (pi->rho + pj->rho);  
          float v_sig_cond = cond_signal_velocity;//sqrtf(fabs(pressurei - pressurej) / mean_rho);//vtilde_signal_velocity;//0.5f * (ci + cj);//vtilde_signal_velocity;//sqrtf(fabs(pressurei - pressurej) / mean_rho);  
          float mean_G = 0.5f * sqrtf((kernel_gradient_i[0] + kernel_gradient_j[0]) * (kernel_gradient_i[0] + kernel_gradient_j[0]) +
                                 (kernel_gradient_i[1] + kernel_gradient_j[1]) * (kernel_gradient_i[1] + kernel_gradient_j[1]) +
                                 (kernel_gradient_i[2] + kernel_gradient_j[2]) * (kernel_gradient_i[2] + kernel_gradient_j[2]));  

          float mean_balsara = 0.5f * (pi->force.balsara + pj->force.balsara);
          
          float alpha_u = mean_balsara;//1.f;//0.1f + mean_balsara * 0.9f;//0.05f;//0.5f;//1.f;//0.1f;//1.f 
          float du_dt_cond_i = -alpha_u * mj * v_sig_cond * (utilde_i - utilde_j) * mean_G / mean_rho;

          pi->u_dt += du_dt_cond_i;

          float alpha_rho = mean_balsara;//1.f;//0.1f + mean_balsara * 0.9f;//0.05f;//0.5f;//1.f;//0.1f; 
          float drho_dt_cond_i = -alpha_rho * mj * (pi->rho / pj->rho) * v_sig_cond * (rhotilde_i - rhotilde_j) * mean_G / mean_rho; 

          pi->drho_dt += drho_dt_cond_i; 
      }
  }
}

#endif /* SWIFT_PLANETARY_HYDRO_IACT_H */
