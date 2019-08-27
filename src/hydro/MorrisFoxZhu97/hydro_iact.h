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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

}

//#ifdef NOPE
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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  float r = sqrtf(r2);
  float xi, xj;
  float h_inv;
  float wi, wj, wi_dx, wj_dx;
  float mi, mj;
  float inv_hidim_pow_plus_one;
  float inv_hjdim_pow_plus_one;
  float dv[3];
  float rhoi, rhoj, rhoi_inv, rhoj_inv;
  float Pi_over_rhosq, Pj_over_rhosq;

/*  if( (pi->is_boundary && !pj->is_boundary) || (pj->is_boundary && !pi->is_boundary)){
    printf("We did a thing!\n");
  }*/

  mi = pi->mass;
  mj = pj->mass;
  rhoi = (float)pi->rho;
  rhoj = (float)pj->rho;
  rhoi_inv = 1.0f / rhoi;
  rhoj_inv = 1.0f / rhoj;

  /* Compute dv */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
 
  inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
  inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);


  Pi_over_rhosq = pi->pressure * (rhoi_inv * rhoi_inv);
  Pj_over_rhosq = pj->pressure * (rhoj_inv * rhoj_inv);

  /* Compute density of pi. */
  h_inv = 1.0f / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dx = wi_dx * inv_hidim_pow_plus_one;

  /* Compute density of pj. */
  h_inv = 1.0 / hj;
  xj = r * h_inv;
  kernel_deval(xj, &wj, &wj_dx);
  wj_dx = wj_dx * inv_hjdim_pow_plus_one;

  /* Compute the artificial viscosity term */
//  printf("pi->visc = %f\n", pi->dynamic_viscosity);
  const float mu_i = rhoi * pi->dynamic_viscosity;
  const float mu_j = rhoj * pj->dynamic_viscosity;
  const float mu_sum = mu_i + mu_j;
  const float dens_ij_inv = rhoi_inv * rhoj_inv;

  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
  const float acc = (Pi_over_rhosq + Pj_over_rhosq);

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];
//  if((pi->id == 86400 && pj->is_boundary) || (pi->is_boundary && pj->id == 86400)) printf("pi_ab = %f, v_dot_r = %f\n", pi_ab, v_dot_r);
/*  if( (!pj->is_boundary) || (!pi->is_boundary)){
    printf("r2 = %f, pressure = %f, xi=%f, rho=%f, acc=%f, wi_dx=%f, dx[1]=%f drho_dt=%f\n",r2,pi->pressure,xi, pi->rho, acc, wi_dx, dx[1], mj*dens);
    printf("hydro_i = %f, hydro_j = %f\n", -mj*acc*wi_dx*dx[1], mi*acc*wj_dx*dx[1]);
}*/
//  if(pi->id == 6144 || pj->id == 6144) printf("v_dot_r = %f, pi_ab = %f\n", v_dot_r, pi_ab); 
//  printf("dens=%f, dv[0]=%f, dx[0]=%f\n", dens,dv[0],dx[0]);
  pi->drho_dt += mj * dens;
  pi->a_hydro[0] -= mj * acc * wi_dx * dx[0];
  pi->a_hydro[0] += (mj * mu_sum * dv[0] * dens_ij_inv) * wi_dx /r;
  pi->a_hydro[1] -= mj * acc * wi_dx * dx[1];
  pi->a_hydro[1] += (mj * mu_sum * dv[1] * dens_ij_inv) * wi_dx /r;
  pi->a_hydro[2] -= mj * acc * wi_dx * dx[2];
  pi->a_hydro[2] += (mj * mu_sum * dv[2] * dens_ij_inv) * wi_dx /r;

  /* Compute density of pj. */
  pj->drho_dt +=  mi * dens;
  pj->a_hydro[0] += mi * acc * wj_dx * dx[0];
  pj->a_hydro[0] -= (mi * mu_sum * dv[0] * dens_ij_inv) * wj_dx /r;
  pj->a_hydro[1] += mi * acc * wj_dx * dx[1];
  pj->a_hydro[1] -= (mi * mu_sum * dv[1] * dens_ij_inv) * wj_dx /r;
  pj->a_hydro[2] += mi * acc * wj_dx * dx[2];
  pj->a_hydro[2] -= (mi * mu_sum * dv[2] * dens_ij_inv) * wj_dx /r;
  
//if(pi->id == 65744) printf("hydro = %e, visco = %e, sum = %e\n", -mj * acc * wi_dx * dx[0], (mj * mu_sum * dv[0] * dens_ij_inv) * wi_dx /r, pi->a_hydro[0]);
//if(pj->id == 65744) printf("hydro = %e, visco = %e, sum = %e\n", mi * acc * wj_dx * dx[0], -(mi * mu_sum * dv[0] * dens_ij_inv) * wj_dx /r, pj->a_hydro[0]);

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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    const struct part *restrict pj, float a, float H) {

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  float r = sqrtf(r2);
  float xi;
  float h_inv;
  float wi, wi_dx;
  float mj;
  float inv_hidim_pow_plus_one;
  //float inv_hjdim_pow_plus_one;
  float dv[3];
  float rhoi, rhoj, rhoi_inv, rhoj_inv;
  float Pi_over_rhosq, Pj_over_rhosq;

  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  rhoi_inv = 1.0f / rhoi;
  rhoj_inv = 1.0f / rhoj;

  /* Compute dv */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
 // inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);


  Pi_over_rhosq = pi->pressure * (rhoi_inv * rhoi_inv);
  Pj_over_rhosq = pj->pressure * (rhoj_inv * rhoj_inv);

  /* Compute density of pi. */
  h_inv = 1.0f / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dx = wi_dx * inv_hidim_pow_plus_one;

  /* Compute the artificial viscosity term */
  const float mu_i = rhoi * pi->dynamic_viscosity;
  const float mu_j = rhoj * pj->dynamic_viscosity;
  const float mu_sum = mu_i + mu_j;
  const float dens_ij_inv = rhoi_inv * rhoj_inv;


//  if((pi->id == 86400 && pj->is_boundary)) printf("pi_ab = %f, v_dot_r = %f\n", pi_ab, v_dot_r);
  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
  const float acc = (Pi_over_rhosq + Pj_over_rhosq);

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];
//  printf("dens=%f dv[0]=%f dx[0]=%f \n", dens, dv[0], dx[0]);
  pi->drho_dt += mj * dens;
  pi->a_hydro[0] -= mj * acc * wi_dx * dx[0];
  pi->a_hydro[0] += (mj * mu_sum * dv[0] * dens_ij_inv) * wi_dx /r;
  pi->a_hydro[1] -= mj * acc * wi_dx * dx[1];
  pi->a_hydro[1] += (mj * mu_sum * dv[1] * dens_ij_inv) * wi_dx /r;
  pi->a_hydro[2] -= mj * acc * wi_dx * dx[2];
  pi->a_hydro[2] += (mj * mu_sum * dv[2] * dens_ij_inv) * wi_dx /r;

}
//#endif

/**
 * @brief Timestep limiter loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

  /* Nothing to do here if both particles are active */
}

/**
 * @brief Timestep limiter loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_limiter(
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj, float a, float H) {

}

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
