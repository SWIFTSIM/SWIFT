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
#ifndef SWIFT_RUNNER_IACT_MINIMAL_H
#define SWIFT_RUNNER_IACT_MINIMAL_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "part.h"
#include "vector.h"

/**
 * @brief Minimal conservative implementation of SPH
 *
 * The thermal variable is the internal energy (u). No viscosity nor
 * thermal conduction terms are implemented.
 */

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float r = sqrtf(r2);

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Compute density of pi. */
  const float hi_inv = 1.f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (3.f * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  /* Compute density of pj. */
  const float hj_inv = 1.f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->rho_dh -= mi * (3.f * wj + xj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= xj * wj_dx;
}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);

  const float h_inv = 1.f / hi;
  const float xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (3.f * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;
}

/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute dv dot r. */
  float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];
  dvdr *= r_inv;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  const float omega_ij = fminf(dvdr, 0.f);

  /* Compute sound speeds */
  const float ci = sqrtf(const_hydro_gamma * pressurei / rhoi);
  const float cj = sqrtf(const_hydro_gamma * pressurej / rhoj);
  const float v_sig = ci + cj + 3.f * omega_ij;

  /* SPH acceleration term */
  const float sph_term = (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_term * dx[0];
  pi->a_hydro[1] -= mj * sph_term * dx[1];
  pi->a_hydro[2] -= mj * sph_term * dx[2];

  pj->a_hydro[0] += mi * sph_term * dx[0];
  pj->a_hydro[1] += mi * sph_term * dx[1];
  pj->a_hydro[2] += mi * sph_term * dx[2];

  /* Get the time derivative for u. */
  pi->u_dt += P_over_rho_i * mj * dvdr * wi_dr;
  pj->u_dt += P_over_rho_j * mi * dvdr * wj_dr;

  /* Get the time derivative for h. */
  pi->h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);
}

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Recover some data */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float xi = r * hi_inv;
  float wi, wi_dx;
  kernel_deval(xi, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  float wj, wj_dx;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute dv dot r. */
  float dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
               (pi->v[2] - pj->v[2]) * dx[2];
  dvdr *= r_inv;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  const float omega_ij = fminf(dvdr, 0.f);

  /* Compute sound speeds */
  const float ci = sqrtf(const_hydro_gamma * pressurei / rhoi);
  const float cj = sqrtf(const_hydro_gamma * pressurej / rhoj);
  const float v_sig = ci + cj + 3.f * omega_ij;

  /* SPH acceleration term */
  const float sph_term = (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * sph_term * dx[0];
  pi->a_hydro[1] -= mj * sph_term * dx[1];
  pi->a_hydro[2] -= mj * sph_term * dx[2];

  /* Get the time derivative for u. */
  pi->u_dt += P_over_rho_i * mj * dvdr * wi_dr;

  /* Get the time derivative for h. */
  pi->h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
}

#endif /* SWIFT_RUNNER_IACT_MINIMAL_H */
