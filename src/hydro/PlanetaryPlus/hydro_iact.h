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
#include "hydro_parameters.h"
#include "math.h"
#include "minmax.h"
#include "hydro_density_estimate.h"
#include "hydro_kernels_etc.h"
#include "hydro_viscosity.h"
#include "hydro_misc_utils.h"

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
    
  hydro_runner_iact_density_extra_density_estimate(pi, pj, dx, wi, wj, wi_dx, wj_dx);
  // hydro_runner_iact_density_extra_kernels(pi, pj, dx, wi, wj, wi_dx, wj_dx); // This will eventually need to be here with new methods
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

  hydro_runner_iact_nonsym_density_extra_density_estimate(pi, pj, dx, wi, wi_dx);
  // hydro_runner_iact_nonsym_density_extra_kernel(pi, pj, dx, wi, wi_dx); // This will eventually need to be here with new methods
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

  hydro_runner_iact_gradient_extra_density_estimate(pi, pj, dx, wi, wj, wi_dx, wj_dx);
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
    const struct part *restrict pj, float a, float H) {

  float wi, wi_dx;

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Compute kernel of pi. */
  const float h_inv = 1.f / hi;
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Correction factors for kernel gradients */
  // Nominator and denominator for f_gdf factors, Wadsley+2017 Eqn. 8
  const float rho_inv_j = 1.f / pj->rho;
  pi->weighted_wcount += pj->mass * r2 * wi_dx * r_inv;
  pi->weighted_neighbour_wcount += pj->mass * r2 * wi_dx * rho_inv_j * r_inv;

  hydro_runner_iact_nonsym_gradient_extra_density_estimate(pi, pj, dx, wi, wi_dx);
  hydro_runner_iact_nonsym_gradient_extra_kernel(pi, pj, dx, wi, wi_dx);
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

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;

  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij;
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Compute the G and Q gradient and viscosity factors, either using matrix
   * inversions or standard GDF SPH */
    
  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley+2017
  or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  hydro_set_Gi_Gj(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);
    

  /* Density factors for GDF or standard equations */
  float rho_factor_i, rho_factor_j;
  hydro_set_rho_factors(&rho_factor_i, &rho_factor_j, pi, pj);

  /* Calculate the viscous pressures Q */  
  float Qi, Qj;
  hydro_set_Qi_Qj(&Qi, &Qj, pi, pj, dx, a, H);

  /* set kernel gradient terms to be used in eolution equations */  
  float kernel_gradient_i[3], kernel_gradient_j[3];
  float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];  
  hydro_set_kernel_gradient_terms(kernel_gradient_i, kernel_gradient_j, Q_kernel_gradient_i, Q_kernel_gradient_j, Gi, Gj);
      
    
    
  float  sigma_dot_kernel_gradient_i[3], sigma_dot_kernel_gradient_j[3], Q_term_i[3], Q_term_j[3];
  for (int i = 0; i < 3; i++) {
      Q_term_i[i] = -Qi * Q_kernel_gradient_i[i];
      Q_term_j[i] = -Qj * Q_kernel_gradient_j[i];
      
      sigma_dot_kernel_gradient_i[i] = 0.f;
      sigma_dot_kernel_gradient_j[i]  = 0.f;   
      for (int j = 0; j < 3; j++) {
          sigma_dot_kernel_gradient_i[i] += pi->stress_tensor_sigma[i][j] * kernel_gradient_i[j];
          sigma_dot_kernel_gradient_j[i] += pj->stress_tensor_sigma[i][j] * kernel_gradient_j[j];
      }
      
  }  
    
  /* Use the force Luke! */  
  pi->a_hydro[0] += mj * ((sigma_dot_kernel_gradient_i[0]  + Q_term_i[0]) / rho_factor_i + (sigma_dot_kernel_gradient_j[0] + Q_term_j[0]) / rho_factor_j);
  pi->a_hydro[1] += mj * ((sigma_dot_kernel_gradient_i[1]  + Q_term_i[1]) / rho_factor_i + (sigma_dot_kernel_gradient_j[1] + Q_term_j[1]) / rho_factor_j);
  pi->a_hydro[2] += mj * ((sigma_dot_kernel_gradient_i[2]  + Q_term_i[2]) / rho_factor_i + (sigma_dot_kernel_gradient_j[2] + Q_term_j[2]) / rho_factor_j);

  pj->a_hydro[0] -= mi * ((sigma_dot_kernel_gradient_i[0]  + Q_term_i[0]) / rho_factor_i + (sigma_dot_kernel_gradient_j[0] + Q_term_j[0]) / rho_factor_j);
  pj->a_hydro[1] -= mi * ((sigma_dot_kernel_gradient_i[1]  + Q_term_i[1]) / rho_factor_i + (sigma_dot_kernel_gradient_j[1] + Q_term_j[1]) / rho_factor_j);
  pj->a_hydro[2] -= mi * ((sigma_dot_kernel_gradient_i[2]  + Q_term_i[2]) / rho_factor_i + (sigma_dot_kernel_gradient_j[2] + Q_term_j[2]) / rho_factor_j);
    
    

  /* dx dot kernel gradient term needed for du/dt in e.g. eq 13 of Wadsley and
   * equiv*/
  const float dv_dot_sigma_dot_kernel_gradient_i = (pi->v[0] - pj->v[0]) * sigma_dot_kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * sigma_dot_kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * sigma_dot_kernel_gradient_i[2];
  const float dv_dot_sigma_dot_kernel_gradient_j = (pi->v[0] - pj->v[0]) * sigma_dot_kernel_gradient_j[0] +
                       (pi->v[1] - pj->v[1]) * sigma_dot_kernel_gradient_j[1] +
                       (pi->v[2] - pj->v[2]) * sigma_dot_kernel_gradient_j[2];
  const float dv_dot_Q_term_i = (pi->v[0] - pj->v[0]) * Q_term_i[0] +
                       (pi->v[1] - pj->v[1]) * Q_term_i[1] +
                       (pi->v[2] - pj->v[2]) * Q_term_i[2];
  const float dv_dot_Q_term_j = (pi->v[0] - pj->v[0]) * Q_term_j[0] +
                       (pi->v[1] - pj->v[1]) * Q_term_j[1] +
                       (pi->v[2] - pj->v[2]) * Q_term_j[2];
    

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = (dv_dot_sigma_dot_kernel_gradient_i + dv_dot_Q_term_i) / rho_factor_i;
  float du_dt_j = (dv_dot_sigma_dot_kernel_gradient_j + dv_dot_Q_term_j) / rho_factor_j;

#ifdef PLANETARY_FIXED_ENTROPY
  du_dt_i = dv_dot_sigma_dot_kernel_gradient_i / rho_factor_i;
  du_dt_j = dv_dot_sigma_dot_kernel_gradient_j / rho_factor_j;
#endif

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;
    
  const float dvdG_i = (pi->v[0] - pj->v[0]) * kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_i[2];
  const float dvdG_j = (pi->v[0] - pj->v[0]) * kernel_gradient_j[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_j[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_j[2];  

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdG_i / rhoj;
  pj->force.h_dt -= mi * dvdG_j / rhoi;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pj->n_force += wi + wj;
  pi->N_force++;
  pj->N_force++;
#endif
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
  const float mj = pj->mass;
  const float rhoj = pj->rho;

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
    
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;

  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij;
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Compute the G and Q gradient and viscosity factors, either using matrix
   * inversions or standard GDF SPH */
  
  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley+2017
  or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  hydro_set_Gi_Gj(Gi, Gj, pi, pj, dx, wi, wj, wi_dx, wj_dx);
    

  /* Density factors for GDF or standard equations */
  float rho_factor_i, rho_factor_j;
  hydro_set_rho_factors(&rho_factor_i, &rho_factor_j, pi, pj);

  /* Calculate the viscous pressures Q */  
  float Qi, Qj;
  hydro_set_Qi_Qj(&Qi, &Qj, pi, pj, dx, a, H);

  /* set kernel gradient terms to be used in eolution equations */  
  float kernel_gradient_i[3], kernel_gradient_j[3];
  float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];  
  hydro_set_kernel_gradient_terms(kernel_gradient_i, kernel_gradient_j, Q_kernel_gradient_i, Q_kernel_gradient_j, Gi, Gj);


  float  sigma_dot_kernel_gradient_i[3], sigma_dot_kernel_gradient_j[3], Q_term_i[3], Q_term_j[3];
  for (int i = 0; i < 3; i++) {
      Q_term_i[i] = -Qi * Q_kernel_gradient_i[i];
      Q_term_j[i] = -Qj * Q_kernel_gradient_j[i];
      
      sigma_dot_kernel_gradient_i[i] = 0.f;
      sigma_dot_kernel_gradient_j[i]  = 0.f;   
      for (int j = 0; j < 3; j++) {
          sigma_dot_kernel_gradient_i[i] += pi->stress_tensor_sigma[i][j] * kernel_gradient_i[j];
          sigma_dot_kernel_gradient_j[i] += pj->stress_tensor_sigma[i][j] * kernel_gradient_j[j];
      }
      
  }  
    
  /* Use the force Luke! */  
  pi->a_hydro[0] += mj * ((sigma_dot_kernel_gradient_i[0]  + Q_term_i[0]) / rho_factor_i + (sigma_dot_kernel_gradient_j[0] + Q_term_j[0]) / rho_factor_j);
  pi->a_hydro[1] += mj * ((sigma_dot_kernel_gradient_i[1]  + Q_term_i[1]) / rho_factor_i + (sigma_dot_kernel_gradient_j[1] + Q_term_j[1]) / rho_factor_j);
  pi->a_hydro[2] += mj * ((sigma_dot_kernel_gradient_i[2]  + Q_term_i[2]) / rho_factor_i + (sigma_dot_kernel_gradient_j[2] + Q_term_j[2]) / rho_factor_j);    
    

  /* dx dot kernel gradient term needed for du/dt in e.g. eq 13 of Wadsley and
   * equiv*/
  const float dv_dot_sigma_dot_kernel_gradient_i = (pi->v[0] - pj->v[0]) * sigma_dot_kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * sigma_dot_kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * sigma_dot_kernel_gradient_i[2];
  const float dv_dot_Q_term_i = (pi->v[0] - pj->v[0]) * Q_term_i[0] +
                       (pi->v[1] - pj->v[1]) * Q_term_i[1] +
                       (pi->v[2] - pj->v[2]) * Q_term_i[2];
    

  /* Get the time derivative for u, including the viscosity */
  float du_dt_i = (dv_dot_sigma_dot_kernel_gradient_i + dv_dot_Q_term_i) / rho_factor_i;

#ifdef PLANETARY_FIXED_ENTROPY
  du_dt_i = dv_dot_sigma_dot_kernel_gradient_i / rho_factor_i;
#endif

  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
    
  const float dvdG_i = (pi->v[0] - pj->v[0]) * kernel_gradient_i[0] +
                       (pi->v[1] - pj->v[1]) * kernel_gradient_i[1] +
                       (pi->v[2] - pj->v[2]) * kernel_gradient_i[2];


  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdG_i / rhoj;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  pi->n_force += wi + wj;
  pi->N_force++;
#endif
}

#endif /* SWIFT_PLANETARY_HYDRO_IACT_H */
