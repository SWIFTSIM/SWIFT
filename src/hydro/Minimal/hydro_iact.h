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
 * viscosity term without switches is implemented. No thermal conduction
 * term is implemented.
 *
 * This corresponds to equations (43), (44), (45), (101), (103)  and (104) with
 * \f$\beta=3\f$ and \f$\alpha_u=0\f$ of Price, D., Journal of Computational
 * Physics, 2012, Volume 231, Issue 3, pp. 759-794.
 */

#include "adiabatic_index.h"
#include "minmax.h"

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
  const float ui = r * h_inv;
  kernel_deval(ui, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);
}

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * pi->force.f;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * pj->force.f;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

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
  const float visc_du_term = 0.5f * visc_acc_term * dvdr;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;
  const float du_dt_j = sph_du_term_j + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;
  pj->u_dt += du_dt_j * mi;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  const float fac_mu = 1.f; /* Will change with cosmological integration */

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

  /* Compute gradient terms */
  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * pi->force.f;
  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * pj->force.f;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Compute sound speeds and signal velocity */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.5f * const_viscosity_alpha * v_sig * mu_ij / rho_ij;

  /* Convolve with the kernel */
  const float visc_acc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;

  /* SPH acceleration term */
  const float sph_acc_term =
      (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

  /* Assemble the acceleration */
  const float acc = sph_acc_term + visc_acc_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for u. */
  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

  /* Viscosity term */
  const float visc_du_term = 0.5f * visc_acc_term * dvdr;

  /* Assemble the energy equation term */
  const float du_dt_i = sph_du_term_i + visc_du_term;

  /* Internal energy time derivatibe */
  pi->u_dt += du_dt_i * mj;

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
