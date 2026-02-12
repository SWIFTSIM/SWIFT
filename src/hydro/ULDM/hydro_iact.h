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
  pi->density.rot_v[1] += faci * curlvr[1];      // Check the sign
  pi->density.rot_v[2] += faci * curlvr[2];

  pj->density.rot_v[0] -= facj * curlvr[0];
  pj->density.rot_v[1] -= facj * curlvr[1];
  pj->density.rot_v[2] -= facj * curlvr[2];
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
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  float wi, wj, wi_dx, wj_dx, wi_d2x, wj_d2x;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get the masses and densities */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float sqrtrhoij_inv = 1/sqrtf(rhoj*rhoi);
  const float sqrtrhoij_inv_r_inv = sqrtrhoij_inv*r_inv;

  /* Compute quantities for pi. */
  const float hi_inv = 1.f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */  
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  kernel_d2eval(ui, &wi_d2x);
  const float wi_dr = hid_inv * wi_dx;  

  /* Density gradient */  
  pi->grad_rho[0] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[0];
  pi->grad_rho[1] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[1];
  pi->grad_rho[2] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[2];

  /* Lapacien of the gradient */
  pi->laplacian_rho += mj * (wi_d2x + 2/ui * wi_dx) * (rhoj-rhoi) *sqrtrhoij_inv;
    
  /* Compute quantities for pj. */
  const float hj_inv = 1.f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */  
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);
  kernel_d2eval(uj, &wj_d2x);
  const float wj_dr = hjd_inv * wj_dx;
  
  /* Density gradient */
  pj->grad_rho[0] -= mi * (rhoi-rhoj)*sqrtrhoij_inv_r_inv * wj_dr * dx[0];
  pj->grad_rho[1] -= mi * (rhoi-rhoj)*sqrtrhoij_inv_r_inv * wj_dr * dx[1];
  pj->grad_rho[2] -= mi * (rhoi-rhoj)*sqrtrhoij_inv_r_inv * wj_dr * dx[2];

  /* Lapacien of the gradient */
  pj->laplacian_rho += mi * (wj_d2x + 2/uj * wj_dx) * (rhoi-rhoj) *sqrtrhoij_inv;
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
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  float wi, wi_dx, wi_d2x;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Get the masses and densities */
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float sqrtrhoij_inv = 1/sqrtf(rhoj*rhoi);
  const float sqrtrhoij_inv_r_inv = sqrtrhoij_inv*r_inv;

  /* Compute quantities for pi. */
  const float hi_inv = 1.f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */  
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  kernel_d2eval(ui, &wi_d2x);
  const float wi_dr = hid_inv * wi_dx;  

  /* Density gradient */  
  pi->grad_rho[0] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[0];
  pi->grad_rho[1] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[1];
  pi->grad_rho[2] += mj * (rhoj-rhoi)*sqrtrhoij_inv_r_inv *wi_dr * dx[2];

  /* Lapacien of the gradient */
  pi->laplacian_rho += mj * (wi_d2x + 2/ui * wi_dx) * (rhoj-rhoi) *sqrtrhoij_inv;
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
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float rhoi2 = rhoi*rhoi;
  const float rhoj2 = rhoi*rhoj;
  const float lapacian_rhoi = pi->laplacian_rho;
  const float lapacian_rhoj = pj->laplacian_rho;
  const float norm_grad_rhoi2 = pi->norm_grad_rho2;
  const float norm_grad_rhoj2 = pj->norm_grad_rho2;

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

  /* Assemble the acceleration */
  const float QPi = lapacian_rhoj/2.0f/rhoj - norm_grad_rhoj2/4.0f/rhoj2; 
  const float QPj = lapacian_rhoi/2.0f/rhoi - norm_grad_rhoi2/4.0f/rhoi2; 
  
  const float fi = wi_dx/mj/rhoj*QPi*r_inv;
  const float fj = wj_dx/mi/rhoi*QPj*r_inv;

  /* Use the force Luke ! Last part of the acceleration is not included */
  pi->a_hydro[0] -= fi * dx[0];
  pi->a_hydro[1] -= fi * dx[1];   // check sign !!!
  pi->a_hydro[2] -= fi * dx[2];
                    
  pj->a_hydro[0] += fj * dx[0];
  pj->a_hydro[1] += fj * dx[1];
  pj->a_hydro[2] += fj * dx[2];
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
  const float r_inv = r ? 1.0f / r : 0.0f;

  /* Recover some data */
  const float mj = pj->mass;
  const float rhoj = pj->rho;
  const float rhoj2 = rhoj*rhoj;
  const float lapacian_rhoj = pj->laplacian_rho;
  const float norm_grad_rhoj2 = pj->norm_grad_rho2;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);

  /* Assemble the acceleration */
  const float QPi = lapacian_rhoj/2.0f/rhoj - norm_grad_rhoj2/4.0f/rhoj2; 
  
  const float fi = wi_dx/mj/rhoj*QPi*r_inv;

  /* Use the force Luke ! Last part of the acceleration is not included */
  pi->a_hydro[0] -= fi * dx[0];
  pi->a_hydro[1] -= fi * dx[1];   // check sign !!!
  pi->a_hydro[2] -= fi * dx[2];
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
