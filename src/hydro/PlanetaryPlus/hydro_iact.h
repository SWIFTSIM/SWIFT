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

#ifdef PLANETARY_IMBALANCE
  /* Add contribution to kernel averages */
  pi->sum_wij += wi * mj;
  pj->sum_wij += wj * mi;

  /* Add contribution r_ij * mj * sqrt(Wij) */
  if (pi->mat_id == pj->mat_id) {
    pi->sum_rij[0] += -dx[0] * wi * mj;
    pi->sum_rij[1] += -dx[1] * wi * mj;
    pi->sum_rij[2] += -dx[2] * wi * mj;

    pj->sum_rij[0] += dx[0] * wj * mi;
    pj->sum_rij[1] += dx[1] * wj * mi;
    pj->sum_rij[2] += dx[2] * wj * mi;
  }

  if (pi->mat_id != pj->mat_id) {
    pi->sum_rij[0] += dx[0] * wi * mj;
    pi->sum_rij[1] += dx[1] * wi * mj;
    pi->sum_rij[2] += dx[2] * wi * mj;

    pj->sum_rij[0] += -dx[0] * wj * mi;
    pj->sum_rij[1] += -dx[1] * wj * mi;
    pj->sum_rij[2] += -dx[2] * wj * mi;
  }
#endif
    
#ifdef PLANETARY_SMOOTHING_CORRECTION
  pi->drho_dh -= pj->mass * (hydro_dimension * wi + ui * wi_dx);
  pj->drho_dh -= pi->mass * (hydro_dimension * wj + uj * wj_dx);
    
  pi->grad_rho[0] += dx[0]*wi_dx*r_inv*mj;
  pi->grad_rho[1] += dx[1]*wi_dx*r_inv*mj;
  pi->grad_rho[2] += dx[2]*wi_dx*r_inv*mj;

  pj->grad_rho[0] += -dx[0]*wj_dx*r_inv*mi;
  pj->grad_rho[1] += -dx[1]*wj_dx*r_inv*mi;
  pj->grad_rho[2] += -dx[2]*wj_dx*r_inv*mi;      
#endif

#ifdef PLANETARY_QUAD_VISC
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;

  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float wj_dr = hjd_inv * wj_dx;

  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Inverse of D matrix (eq 20 in Rosswog 2020) */
      pi->Dinv[i][j] += pj->mass * dx[i] * dx[j] * wi_dr * r_inv;
      pj->Dinv[i][j] += pi->mass * dx[i] * dx[j] * wj_dr * r_inv;

      /* E matrix (second part of eq 19 in Rosswog 2020) */
      pi->E[i][j] += pj->mass * (pi->v[i] - pj->v[i]) * dx[j] * wi_dr * r_inv;
      pj->E[i][j] += pi->mass * (pi->v[i] - pj->v[i]) * dx[j] * wj_dr * r_inv;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Dinv[2][2] = 1.f;
  pj->Dinv[2][2] = 1.f;

#endif
#endif
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

#ifdef PLANETARY_IMBALANCE

  /* Add contribution to kernel averages */
  pi->sum_wij += wi * mj;

  /* Add contribution r_ij*mj*sqrt(Wij) */
  if (pi->mat_id == pj->mat_id) {
    pi->sum_rij[0] += -dx[0] * wi * mj;
    pi->sum_rij[1] += -dx[1] * wi * mj;
    pi->sum_rij[2] += -dx[2] * wi * mj;
  }

  if (pi->mat_id != pj->mat_id) {
    pi->sum_rij[0] += dx[0] * wi * mj;
    pi->sum_rij[1] += dx[1] * wi * mj;
    pi->sum_rij[2] += dx[2] * wi * mj;
  }
#endif
    
#ifdef PLANETARY_SMOOTHING_CORRECTION
  pi->drho_dh -= pj->mass * (hydro_dimension * wi + ui * wi_dx);
    
  pi->grad_rho[0] += dx[0]*wi_dx*r_inv*mj;
  pi->grad_rho[1] += dx[1]*wi_dx*r_inv*mj;
  pi->grad_rho[2] += dx[2]*wi_dx*r_inv*mj;    
#endif

#ifdef PLANETARY_QUAD_VISC
  const float hid_inv = pow_dimension_plus_one(h_inv); /* 1/h^(d+1) */
  const float wi_dr = hid_inv * wi_dx;

  int i, j;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      /* Inverse of D matrix (eq 20 in Rosswog 2020) */
      pi->Dinv[i][j] += pj->mass * dx[i] * dx[j] * wi_dr * r_inv;

      /* E matrix (second part of eq 19 in Rosswog 2020) */
      pi->E[i][j] += pj->mass * (pi->v[i] - pj->v[i]) * dx[j] * wi_dr * r_inv;
    }
  }

#if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Dinv[2][2] = 1.f;

#endif
#endif
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
    

   
#ifdef PLANETARY_IMBALANCE
  /* Compute kernel averages */
  pi->sum_wij_exp += wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_P += pj->P * wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_T += pj->T * wi * expf(-pj->I * pj->I);

  pj->sum_wij_exp += wj * expf(-pi->I * pi->I);
  pj->sum_wij_exp_P += pi->P * wj * expf(-pi->I * pi->I);
  pj->sum_wij_exp_T += pi->T * wj * expf(-pi->I * pi->I);
#endif
    
#ifdef PLANETARY_SMOOTHING_CORRECTION
  float gi = (pi->h / pi->rho) * sqrtf(pi->grad_rho[0] * pi->grad_rho[0] + pi->grad_rho[1] * pi->grad_rho[1] + pi->grad_rho[2] * pi->grad_rho[2]);   
  float gj = (pj->h / pj->rho) * sqrtf(pj->grad_rho[0] * pj->grad_rho[0] + pj->grad_rho[1] * pj->grad_rho[1] + pj->grad_rho[2] * pj->grad_rho[2]);
   
  float f_gi = 0.5f * (1.f + tanhf(3.f - 3.f * gi / (0.5f * gj)));
  float f_gj = 0.5f * (1.f + tanhf(3.f - 3.f * gj / (0.5f * gi)));
    
  pi->P_tilde_numerator += wi * pj->P * f_gj;
  pj->P_tilde_numerator += wj * pi->P * f_gi;   
  pi->P_tilde_denominator += wi * f_gj;
  pj->P_tilde_denominator += wj * f_gi;
    
  pi->max_ngb_sph_rho = max(pi->max_ngb_sph_rho, pj->rho);
  pi->min_ngb_sph_rho = min(pi->min_ngb_sph_rho, pj->rho);
  pj->max_ngb_sph_rho = max(pj->max_ngb_sph_rho, pi->rho);
  pj->min_ngb_sph_rho = min(pj->min_ngb_sph_rho, pi->rho);
    
  pi->grad_drho_dh[0] += (pj->drho_dh - pi->drho_dh) * (dx[0]*wi_dx*r_inv) * (pj->mass / pj->rho);
  pi->grad_drho_dh[1] += (pj->drho_dh - pi->drho_dh) * (dx[1]*wi_dx*r_inv) * (pj->mass / pj->rho);
  pi->grad_drho_dh[2] += (pj->drho_dh - pi->drho_dh) * (dx[2]*wi_dx*r_inv) * (pj->mass / pj->rho);

  pj->grad_drho_dh[0] += (pi->drho_dh - pj->drho_dh) * (-dx[0]*wj_dx*r_inv) * (pi->mass / pi->rho);
  pj->grad_drho_dh[1] += (pi->drho_dh - pj->drho_dh) * (-dx[1]*wj_dx*r_inv) * (pi->mass / pi->rho);
  pj->grad_drho_dh[2] += (pi->drho_dh - pj->drho_dh) * (-dx[2]*wj_dx*r_inv) * (pi->mass / pi->rho);    
#endif


#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC 

  float volume_i = pi->mass * rho_inv_i;
  float volume_j = pj->mass * rho_inv_j;  

    
  #ifdef PLANETARY_SMOOTHING_CORRECTION

      if(pi->last_corrected_rho){
        volume_i = pi->mass / (pi->last_f_S * pi->rho + (1.f - pi->last_f_S) * pi->last_corrected_rho);
      }
      if(pj->last_corrected_rho){
        volume_j = pj->mass / (pj->last_f_S * pj->rho + (1.f - pj->last_f_S) * pj->last_corrected_rho);
      }   
  #endif  
    
    
  pi->sum_w_V += wi * volume_j;// * pj->mass / pj->rho;// 
  pj->sum_w_V += wj * volume_i;// * pi->mass / pi->rho;//
    
  pi->sum_r_w_V[0] += -dx[0] * wi * volume_j;// * pj->mass / pj->rho;//
  pi->sum_r_w_V[1] += -dx[1] * wi * volume_j;// * pj->mass / pj->rho;// * volume_j;//
  pi->sum_r_w_V[2] += -dx[2] * wi * volume_j;// * pj->mass / pj->rho;// * volume_j;//
    
  pj->sum_r_w_V[0] += dx[0] * wj * volume_i;// * pi->mass / pi->rho;// * volume_i;//
  pj->sum_r_w_V[1] += dx[1] * wj * volume_i;// * pi->mass / pi->rho;// * volume_i;//
  pj->sum_r_w_V[2] += dx[2] * wj * volume_i;// * pi->mass / pi->rho;// * volume_i;//    
    
    
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {

      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * volume_j;
      pj->Cinv[i][j] += dx[i] * dx[j] * wj * volume_i;
    }
  }
  
  #if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Cinv[2][2] = 1.f;
  pj->Cinv[2][2] = 1.f;

#endif

#endif

#ifdef PLANETARY_QUAD_VISC 
  
  int k;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {

      /* Gradients from eq 18 in Rosswog 2020 (without C multiplied)*/
      pi->dv[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wi * volume_j;
      pj->dv[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wj * volume_i;

      for (k = 0; k < 3; ++k) {

        /* Gradients from eq 18 in Rosswog 2020 (without C multiplied). Note
         * that we now use dv_aux to get second derivative*/
        pi->ddv[i][j][k] += (pi->dv_aux[i][j] - pj->dv_aux[i][j]) *
                            dx[k] * wi * volume_j;
        pj->ddv[i][j][k] += (pi->dv_aux[i][j] - pj->dv_aux[i][j]) *
                            dx[k] * wj * volume_i;
      }
    }
  }


  /* Number of neighbours. Needed for eta_crit factor in slope limiter */
  pi->N_grad += 1.f;
  pj->N_grad += 1.f;
#endif
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
    
 
#ifdef PLANETARY_IMBALANCE
  /* Compute kernel averages */
  pi->sum_wij_exp += wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_P += pj->P * wi * expf(-pj->I * pj->I);
  pi->sum_wij_exp_T += pj->T * wi * expf(-pj->I * pj->I);
#endif
    
#ifdef PLANETARY_SMOOTHING_CORRECTION
  float gi = (pi->h / pi->rho) * sqrtf(pi->grad_rho[0] * pi->grad_rho[0] + pi->grad_rho[1] * pi->grad_rho[1] + pi->grad_rho[2] * pi->grad_rho[2]);   
  float gj = (pj->h / pj->rho) * sqrtf(pj->grad_rho[0] * pj->grad_rho[0] + pj->grad_rho[1] * pj->grad_rho[1] + pj->grad_rho[2] * pj->grad_rho[2]);
   
  float f_gj = 0.5f * (1.f + tanhf(3.f - 3.f * gj / (0.5f * gi)));
    
  pi->P_tilde_numerator += wi * pj->P * f_gj; 
  pi->P_tilde_denominator += wi * f_gj;
   
  pi->max_ngb_sph_rho = max(pi->max_ngb_sph_rho, pj->rho);
  pi->min_ngb_sph_rho = min(pi->min_ngb_sph_rho, pj->rho);
    
  pi->grad_drho_dh[0] += (pj->drho_dh - pi->drho_dh) * (dx[0]*wi_dx*r_inv) * (pj->mass / pj->rho);
  pi->grad_drho_dh[1] += (pj->drho_dh - pi->drho_dh) * (dx[1]*wi_dx*r_inv) * (pj->mass / pj->rho);
  pi->grad_drho_dh[2] += (pj->drho_dh - pi->drho_dh) * (dx[2]*wi_dx*r_inv) * (pj->mass / pj->rho);  
#endif


#if defined PLANETARY_MATRIX_INVERSION || defined PLANETARY_QUAD_VISC

  float volume_j = pj->mass * rho_inv_j;  

   
  #ifdef PLANETARY_SMOOTHING_CORRECTION
      if(pj->last_corrected_rho){
        volume_j = pj->mass / (pj->last_f_S * pj->rho + (1.f - pj->last_f_S) * pj->last_corrected_rho);
      }   
  #endif    
    
  pi->sum_w_V += wi * volume_j;// * pj->mass / pj->rho;// 
    
  pi->sum_r_w_V[0] += -dx[0] * wi * volume_j;// * pj->mass / pj->rho;//
  pi->sum_r_w_V[1] += -dx[1] * wi * volume_j;// * pj->mass / pj->rho;// * volume_j;//
  pi->sum_r_w_V[2] += -dx[2] * wi * volume_j;// * pj->mass / pj->rho;// * volume_j;//
 
    
  int i, j;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {

      /* Inverse of C matrix (eq 6 in Rosswog 2020) */
      pi->Cinv[i][j] += dx[i] * dx[j] * wi * volume_j;
    }
  }
  
  #if defined(HYDRO_DIMENSION_2D)
  /* This is so we can do 3x3 matrix inverse even when 2D */
  pi->Cinv[2][2] = 1.f;
#endif

#endif

#ifdef PLANETARY_QUAD_VISC 
  
  int k;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {

      /* Gradients from eq 18 in Rosswog 2020 (without C multiplied)*/
      pi->dv[i][j] += (pi->v[i] - pj->v[i]) * dx[j] * wi * volume_j;
      for (k = 0; k < 3; ++k) {

        /* Gradients from eq 18 in Rosswog 2020 (without C multiplied). Note
         * that we now use dv_aux to get second derivative*/
        pi->ddv[i][j][k] += (pi->dv_aux[i][j] - pj->dv_aux[i][j]) *
                            dx[k] * wi * volume_j;
      }
    }
  }


  /* Number of neighbours. Needed for eta_crit factor in slope limiter */
  pi->N_grad += 1.f;
#endif
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

  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley+2017
    or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  /* For loops */
  int i;

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
#ifdef PLANETARY_MATRIX_INVERSION
  if (!pi->is_h_max && !pj->is_h_max && !pi->is_vacuum_boundary && !pj->is_vacuum_boundary) {
    for (i = 0; i < 3; ++i) {
       /* eq 4 and 5 in Rosswog 2020. These replace the gradient of the kernel */
      Gi[i] =
          -(pi->C[i][0] * dx[0] + pi->C[i][1] * dx[1] + pi->C[i][2] * dx[2]) *
          wi;
      Gj[i] =
          -(pj->C[i][0] * dx[0] + pj->C[i][1] * dx[1] + pj->C[i][2] * dx[2]) *
          wj;
    }
  }else{
    for (i = 0; i < 3; ++i) {
      /* If h=h_max use the standard kernel gradients */
      Gi[i] = wi_dr * dx[i] * r_inv;
      Gj[i] = wj_dr * dx[i] * r_inv;
    }
  }
#else

#ifdef PLANETARY_GDF
  /* Standard GDF kernel gradients, Wadsley+2017 Eqn. 7, in Rosswog2020
   * framework */
  /* Include the dx and r_inv here instead of later */
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
    Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  }

#else
  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
    
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * f_ij;
    Gj[i] = wj_dr * dx[i] * r_inv * f_ji;
  }
#endif
#endif
   
#ifdef PLANETARY_QUAD_VISC
  /* Quadratically reconstructed velocities at the halfway point between
   * particles */
  float vtilde_i[3], vtilde_j[3];

  /* Some parameters for artificial visc. Taken from Rosswog 2020 */
  float alpha = 1.f;
  float beta = 2.f;
  float epsilon = 0.1;

  /* Square of eta (eq 16 in Rosswog 2020) */
  float eta_i_2 = r2 * hi_inv * hi_inv;
  float eta_j_2 = r2 * hj_inv * hj_inv;

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!pi->is_h_max && !pj->is_h_max) {

    /* For loops */
    int j, k;

    /* eq 23 in Rosswog 2020 */
    float eta_ab = min(r * hi_inv, r * hj_inv);

    /* A numerators and denominators (eq 22 in Rosswog 2020) */
    float A_i_v = 0.f;
    float A_j_v = 0.f;

    /* Terms in square brackets in Rosswog 2020 eq 17 */
    float v_quad_i[3] = {0};
    float v_quad_j[3] = {0};

    /* eq 23 in Rosswog 2020 for 3D (rearranged to get something like 4/3 pi
     * eta^3) */
    float eta_crit = 2.172975 / cbrt((pi->N_grad + pj->N_grad) * 0.5);

#if defined(HYDRO_DIMENSION_2D)

    /* eq 23 in Rosswog 2020 for 2D (rearranged to get something like pi eta^2)
     */
    eta_crit = 2.754572 / sqrtf((pi->N_grad + pj->N_grad) * 0.5);

#endif

    for (i = 0; i < 3; ++i) {
      
      for (j = 0; j < 3; ++j) {

        /* Get the A numerators and denominators (eq 22 in Rosswog 2020). C_dv
         * is dv from eq 18 */
        A_i_v += pi->C_dv[i][j] * dx[i] * dx[j];
        A_j_v += pj->C_dv[i][j] * dx[i] * dx[j];

        /* Terms in square brackets in Rosswog 2020 eq 17. Add in FIRST
         * derivative terms */
        v_quad_i[j] -= 0.5 * pi->C_dv[i][j] * dx[i];
        v_quad_j[j] += 0.5 * pj->C_dv[i][j] * dx[i];

        for (k = 0; k < 3; ++k) {
          /* Terms in square brackets in Rosswog 2020 eq 17. Add in SECOND
           * derivative terms */
          v_quad_i[j] += 0.125 * pi->C_ddv[i][j][k] * dx[i] * dx[k];
          v_quad_j[j] += 0.125 * pj->C_ddv[i][j][k] * dx[i] * dx[k];
        }
      }
    }

    /* Slope limiter (eq 21 in Rosswog 2020) */
    float phi_i_v =
        min(1.f, 4 * A_i_v / A_j_v / (1 + A_i_v / A_j_v) / (1 + A_i_v / A_j_v));
    phi_i_v = max(0.f, phi_i_v);

    float phi_j_v =
        min(1.f, 4 * A_j_v / A_i_v / (1 + A_j_v / A_i_v) / (1 + A_j_v / A_i_v));
    phi_j_v = max(0.f, phi_j_v);

    if (eta_ab < eta_crit) {
      phi_i_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
      phi_j_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
    }
    /* These are here to catch division by 0. In this case phi tends to 0 anyway
     */
    if (isnan(phi_i_v) || isinf(phi_i_v)) {
      phi_i_v = 0.f;
    }
    if (isnan(phi_j_v) || isinf(phi_j_v)) {
      phi_j_v = 0.f;
    }

    for (i = 0; i < 3; ++i) {
      /* Assemble the reconstructed velocity (eq 17 in Rosswog 2020) */
      vtilde_i[i] = pi->v[i] + phi_i_v * v_quad_i[i];
      vtilde_j[i] = pj->v[i] + phi_j_v * v_quad_j[i];
    }

  } else {

    for (i = 0; i < 3; ++i) {

      /* If h=h_max don't reconstruct velocity */
      vtilde_i[i] = pi->v[i];
      vtilde_j[i] = pj->v[i];
    }
  }

  /* Finally assemble eq 15 in Rosswog 2020 */
  float mu_i = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hi_inv / (eta_i_2 + epsilon * epsilon));
  float mu_j = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hj_inv / (eta_j_2 + epsilon * epsilon));

  /* Get viscous pressure terms (eq 14 in Rosswog 2020) */
  float Q_i = rhoi * (-alpha * ci * mu_i + beta * mu_i * mu_i);
  float Q_j = rhoj * (-alpha * cj * mu_j + beta * mu_j * mu_j);

#else
    
  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara; 
    
#ifdef PLANETARY_GDF

  /* Artificial viscosity terms, as pressure in Rosswog2020 framework, S2.2.1 */
  float Q_i = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoj /
              (rhoi + rhoj);
  float Q_j = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoj /
              (rhoi + rhoj);
#else
  /* Artificial viscosity terms, as pressure in Rosswog2020 framework, S2.2.1 */
  float Q_i = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoi /
              (rhoi + rhoj);
  float Q_j = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoj * rhoj /
              (rhoi + rhoj);   
#endif
    
#endif
    
    float kernel_gradient_i[3], kernel_gradient_j[3];
    float P_i_term, P_j_term;
    float Q_i_term, Q_j_term;

    
#ifdef PLANETARY_GDF
  /* In GDF we use average of Gi and Gj. */
  kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);
    
  kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);

  P_i_term = pressurei / (rhoi * rhoj); 
  P_j_term = pressurej / (rhoi * rhoj);

  Q_i_term = Q_i / (rhoi * rhoj);
  Q_j_term = Q_j / (rhoi * rhoj);
#else 
    
  kernel_gradient_i[0] = Gi[0];
  kernel_gradient_i[1] = Gi[1];
  kernel_gradient_i[2] = Gi[2];
    
  kernel_gradient_j[0] = Gj[0];
  kernel_gradient_j[1] = Gj[1];
  kernel_gradient_j[2] = Gj[2];
    
      
  P_i_term = pressurei / (rhoi * rhoi);
  P_j_term = pressurej / (rhoj * rhoj);
    
  Q_i_term = Q_i / (rhoi * rhoi);
  Q_j_term = Q_j / (rhoj * rhoj);
#endif
    
float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];
#ifdef PLANETARY_QUAD_VISC
  Q_kernel_gradient_i[0] = kernel_gradient_i[0];
  Q_kernel_gradient_i[1] = kernel_gradient_i[1];
  Q_kernel_gradient_i[2] = kernel_gradient_i[2];
    
  Q_kernel_gradient_j[0] = kernel_gradient_j[0];
  Q_kernel_gradient_j[1] = kernel_gradient_j[1];
  Q_kernel_gradient_j[2] = kernel_gradient_j[2];
#else
  Q_kernel_gradient_i[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_i[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_i[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
    
  Q_kernel_gradient_j[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_j[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_j[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
#endif
    
    
    
  /* Use the force Luke! */
  pi->a_hydro[0] -= mj * (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] + Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pi->a_hydro[1] -= mj * (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] + Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pi->a_hydro[2] -= mj * (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] + Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);

  pj->a_hydro[0] += mi * (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] + Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pj->a_hydro[1] += mi * (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] + Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pj->a_hydro[2] += mi * (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] + Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);
    
    
// Use standard SPH P, rho and gradients for adiabatic heating    
#ifdef PLANETARY_SMOOTHING_CORRECTION   
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
    Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  }
    
    #ifdef PLANETARY_GDF
      /* In GDF we use average of Gi and Gj. */
      kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
      kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
      kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);

      kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
      kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
      kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);

      P_i_term = pi->P_sph / (pi->rho_sph * pj->rho_sph);
      P_j_term = pj->P_sph / (pi->rho_sph * pj->rho_sph);  

    #else 

      kernel_gradient_i[0] = Gi[0];
      kernel_gradient_i[1] = Gi[1];
      kernel_gradient_i[2] = Gi[2];

      kernel_gradient_j[0] = Gj[0];
      kernel_gradient_j[1] = Gj[1];
      kernel_gradient_j[2] = Gj[2];


      P_i_term = pi->P_sph / (pi->rho_sph * pi->rho_sph);
      P_j_term = pj->P_sph / (pj->rho_sph * pj->rho_sph);

    #endif    

#endif     
    
  
        
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

    
  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

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

  /* G[3] is the kernel gradient term. Takes the place of eq 7 in Wadsley 2017
   * or the average of eq 4 and 5 in Rosswog 2020 (as described below eq 11) */
  float Gj[3], Gi[3];
  /* For loops */
  int i;

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
#ifdef PLANETARY_MATRIX_INVERSION
  if (!pi->is_h_max && !pj->is_h_max && !pi->is_vacuum_boundary && !pj->is_vacuum_boundary) {
    for (i = 0; i < 3; ++i) {
       /* eq 4 and 5 in Rosswog 2020. These replace the gradient of the kernel */
      Gi[i] =
          -(pi->C[i][0] * dx[0] + pi->C[i][1] * dx[1] + pi->C[i][2] * dx[2]) *
          wi;
      Gj[i] =
          -(pj->C[i][0] * dx[0] + pj->C[i][1] * dx[1] + pj->C[i][2] * dx[2]) *
          wj;
    }
  }else{
    for (i = 0; i < 3; ++i) {
      /* If h=h_max use the standard kernel gradients */
      Gi[i] = wi_dr * dx[i] * r_inv;
      Gj[i] = wj_dr * dx[i] * r_inv;
    }
  }
#else

#ifdef PLANETARY_GDF
  /* Standard GDF kernel gradients, Wadsley+2017 Eqn. 7, in Rosswog2020
   * framework */
  /* Include the dx and r_inv here instead of later */
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
    Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  }

#else
  /* Variable smoothing length term */
  const float f_ij = 1.f - pi->force.f / mj;
  const float f_ji = 1.f - pj->force.f / mi;
    
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * f_ij;
    Gj[i] = wj_dr * dx[i] * r_inv * f_ji;
  }
#endif
#endif
   
#ifdef PLANETARY_QUAD_VISC
  /* Quadratically reconstructed velocities at the halfway point between
   * particles */
  float vtilde_i[3], vtilde_j[3];

  /* Some parameters for artificial visc. Taken from Rosswog 2020 */
  float alpha = 1.f;
  float beta = 2.f;
  float epsilon = 0.1;

  /* Square of eta (eq 16 in Rosswog 2020) */
  float eta_i_2 = r2 * hi_inv * hi_inv;
  float eta_j_2 = r2 * hj_inv * hj_inv;

  /* If h=h_max don't do anything fancy. Things like using m/rho to calculate
   * the volume stops working */
  if (!pi->is_h_max && !pj->is_h_max) {

    /* For loops */
    int j, k;

    /* eq 23 in Rosswog 2020 */
    float eta_ab = min(r * hi_inv, r * hj_inv);

    /* A numerators and denominators (eq 22 in Rosswog 2020) */
    float A_i_v = 0.f;
    float A_j_v = 0.f;

    /* Terms in square brackets in Rosswog 2020 eq 17 */
    float v_quad_i[3] = {0};
    float v_quad_j[3] = {0};

    /* eq 23 in Rosswog 2020 for 3D (rearranged to get something like 4/3 pi
     * eta^3) */
    float eta_crit = 2.172975 / cbrt((pi->N_grad + pj->N_grad) * 0.5);

#if defined(HYDRO_DIMENSION_2D)

    /* eq 23 in Rosswog 2020 for 2D (rearranged to get something like pi eta^2)
     */
    eta_crit = 2.754572 / sqrtf((pi->N_grad + pj->N_grad) * 0.5);

#endif

    for (i = 0; i < 3; ++i) {
      
      for (j = 0; j < 3; ++j) {

        /* Get the A numerators and denominators (eq 22 in Rosswog 2020). C_dv
         * is dv from eq 18 */
        A_i_v += pi->C_dv[i][j] * dx[i] * dx[j];
        A_j_v += pj->C_dv[i][j] * dx[i] * dx[j];

        /* Terms in square brackets in Rosswog 2020 eq 17. Add in FIRST
         * derivative terms */
        v_quad_i[j] -= 0.5 * pi->C_dv[i][j] * dx[i];
        v_quad_j[j] += 0.5 * pj->C_dv[i][j] * dx[i];

        for (k = 0; k < 3; ++k) {
          /* Terms in square brackets in Rosswog 2020 eq 17. Add in SECOND
           * derivative terms */
          v_quad_i[j] += 0.125 * pi->C_ddv[i][j][k] * dx[i] * dx[k];
          v_quad_j[j] += 0.125 * pj->C_ddv[i][j][k] * dx[i] * dx[k];
        }
      }
    }

    /* Slope limiter (eq 21 in Rosswog 2020) */
    float phi_i_v =
        min(1.f, 4 * A_i_v / A_j_v / (1 + A_i_v / A_j_v) / (1 + A_i_v / A_j_v));
    phi_i_v = max(0.f, phi_i_v);

    float phi_j_v =
        min(1.f, 4 * A_j_v / A_i_v / (1 + A_j_v / A_i_v) / (1 + A_j_v / A_i_v));
    phi_j_v = max(0.f, phi_j_v);

    if (eta_ab < eta_crit) {
      phi_i_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
      phi_j_v *= exp(-(eta_ab - eta_crit) * (eta_ab - eta_crit) * 25);
    }
    /* These are here to catch division by 0. In this case phi tends to 0 anyway
     */
    if (isnan(phi_i_v) || isinf(phi_i_v)) {
      phi_i_v = 0.f;
    }
    if (isnan(phi_j_v) || isinf(phi_j_v)) {
      phi_j_v = 0.f;
    }

    for (i = 0; i < 3; ++i) {
      /* Assemble the reconstructed velocity (eq 17 in Rosswog 2020) */
      vtilde_i[i] = pi->v[i] + phi_i_v * v_quad_i[i];
      vtilde_j[i] = pj->v[i] + phi_j_v * v_quad_j[i];
    }

  } else {

    for (i = 0; i < 3; ++i) {

      /* If h=h_max don't reconstruct velocity */
      vtilde_i[i] = pi->v[i];
      vtilde_j[i] = pj->v[i];
    }
  }

  /* Finally assemble eq 15 in Rosswog 2020 */
  float mu_i = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hi_inv / (eta_i_2 + epsilon * epsilon));
  float mu_j = min(0.f, ((vtilde_i[0] - vtilde_j[0]) * dx[0] +
                         (vtilde_i[1] - vtilde_j[1]) * dx[1] +
                         (vtilde_i[2] - vtilde_j[2]) * dx[2]) *
                            hj_inv / (eta_j_2 + epsilon * epsilon));

  /* Get viscous pressure terms (eq 14 in Rosswog 2020) */
  float Q_i = rhoi * (-alpha * ci * mu_i + beta * mu_i * mu_i);
  float Q_j = rhoj * (-alpha * cj * mu_j + beta * mu_j * mu_j);

#else
    
  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara; 
    
#ifdef PLANETARY_GDF

  /* Artificial viscosity terms, as pressure in Rosswog2020 framework, S2.2.1 */
  float Q_i = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoj /
              (rhoi + rhoj);
  float Q_j = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoj /
              (rhoi + rhoj);
#else
  /* Artificial viscosity terms, as pressure in Rosswog2020 framework, S2.2.1 */
  float Q_i = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoi * rhoi /
              (rhoi + rhoj);
  float Q_j = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) * rhoj * rhoj /
              (rhoi + rhoj);   
#endif
    
#endif
    
    float kernel_gradient_i[3], kernel_gradient_j[3];
    float P_i_term, P_j_term;
    float Q_i_term, Q_j_term;

    
#ifdef PLANETARY_GDF
  /* In GDF we use average of Gi and Gj. */
  kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);
    
  kernel_gradient_j[0] = 0.5f * (Gi[0] + Gj[0]);
  kernel_gradient_j[1] = 0.5f * (Gi[1] + Gj[1]);
  kernel_gradient_j[2] = 0.5f * (Gi[2] + Gj[2]);

  P_i_term = pressurei / (rhoi * rhoj); 
  P_j_term = pressurej / (rhoi * rhoj);

  Q_i_term = Q_i / (rhoi * rhoj);
  Q_j_term = Q_j / (rhoi * rhoj);
   
#else 
    
  kernel_gradient_i[0] = Gi[0];
  kernel_gradient_i[1] = Gi[1];
  kernel_gradient_i[2] = Gi[2];
    
  kernel_gradient_j[0] = Gj[0];
  kernel_gradient_j[1] = Gj[1];
  kernel_gradient_j[2] = Gj[2];
    
      
  P_i_term = pressurei / (rhoi * rhoi);
  P_j_term = pressurej / (rhoj * rhoj);
    
  Q_i_term = Q_i / (rhoi * rhoi);
  Q_j_term = Q_j / (rhoj * rhoj);
      
#endif    
    
float Q_kernel_gradient_i[3], Q_kernel_gradient_j[3];
#ifdef PLANETARY_QUAD_VISC
  Q_kernel_gradient_i[0] = kernel_gradient_i[0];
  Q_kernel_gradient_i[1] = kernel_gradient_i[1];
  Q_kernel_gradient_i[2] = kernel_gradient_i[2];
    
  Q_kernel_gradient_j[0] = kernel_gradient_j[0];
  Q_kernel_gradient_j[1] = kernel_gradient_j[1];
  Q_kernel_gradient_j[2] = kernel_gradient_j[2];
#else
  Q_kernel_gradient_i[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_i[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_i[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
    
  Q_kernel_gradient_j[0] = 0.5f * (kernel_gradient_i[0] + kernel_gradient_j[0]);
  Q_kernel_gradient_j[1] = 0.5f * (kernel_gradient_i[1] + kernel_gradient_j[1]);
  Q_kernel_gradient_j[2] = 0.5f * (kernel_gradient_i[2] + kernel_gradient_j[2]);
#endif
    
    
    
  /* Use the force Luke! */
  pi->a_hydro[0] -= mj * (P_i_term * kernel_gradient_i[0] + P_j_term * kernel_gradient_j[0] + Q_i_term * Q_kernel_gradient_i[0] + Q_j_term * Q_kernel_gradient_j[0]);
  pi->a_hydro[1] -= mj * (P_i_term * kernel_gradient_i[1] + P_j_term * kernel_gradient_j[1] + Q_i_term * Q_kernel_gradient_i[1] + Q_j_term * Q_kernel_gradient_j[1]);
  pi->a_hydro[2] -= mj * (P_i_term * kernel_gradient_i[2] + P_j_term * kernel_gradient_j[2] + Q_i_term * Q_kernel_gradient_i[2] + Q_j_term * Q_kernel_gradient_j[2]);
    


// Use standard SPH P, rho and gradients for adiabatic heating    
#ifdef PLANETARY_SMOOTHING_CORRECTION   
  for (i = 0; i < 3; i++) {
    Gi[i] = wi_dr * dx[i] * r_inv * pi->f_gdf;
    Gj[i] = wj_dr * dx[i] * r_inv * pj->f_gdf;
  }
    
    #ifdef PLANETARY_GDF
      /* In GDF we use average of Gi and Gj. */
      kernel_gradient_i[0] = 0.5f * (Gi[0] + Gj[0]);
      kernel_gradient_i[1] = 0.5f * (Gi[1] + Gj[1]);
      kernel_gradient_i[2] = 0.5f * (Gi[2] + Gj[2]);

      P_i_term = pi->P_sph / (pi->rho_sph * pj->rho_sph);

    #else 

      kernel_gradient_i[0] = Gi[0];
      kernel_gradient_i[1] = Gi[1];
      kernel_gradient_i[2] = Gi[2];

      P_i_term = pi->P_sph / (pi->rho_sph * pi->rho_sph);

    #endif    

#endif    
    
    
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
    
  /* Internal energy time derivative */
  pi->u_dt += du_dt_i * mj;

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
