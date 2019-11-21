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
#include "boundary.h"

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
  const int r_inv = 1.0 / r;

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

/* Now we can do the calculation of the equations of motion */

  /* Do the 'regular' hydro forces... */
  const float acc_sph_i = (Pi_over_rhosq + Pj_over_rhosq) * r_inv * wi_dx;
  const float acc_sph_j = (Pi_over_rhosq + Pj_over_rhosq) * r_inv * wj_dx;
  
  /* Viscosity implementation */
  const float mu_i = rhoi * pi->dynamic_viscosity;
  const float mu_j = rhoj * pj->dynamic_viscosity;
  const float mu_sum = mu_i + mu_j;
  const float dens_ij_inv = rhoi_inv * rhoj_inv;
  const float acc_visc_i = mu_sum * wi_dx * dens_ij_inv * r_inv;
  const float acc_visc_j = mu_sum * wj_dx * dens_ij_inv * r_inv;
  const float dens_i = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];
  const float dens_j = dv[0] * wj_dx*dx[0] + dv[1] * wj_dx*dx[1] + dv[2] * wj_dx*dx[2];
  
  /* Apply forces */
  pi->drho_dt += mj * dens_i;
  pi->a_hydro[0] += mj * (acc_visc_i * dv[0] - acc_sph_i * dx[0]);
  pi->a_hydro[1] += mj * (acc_visc_i * dv[1] - acc_sph_i * dx[1]);
  pi->a_hydro[2] += mj * (acc_visc_i * dv[2] - acc_sph_i * dx[2]);

  pi->a_visc[0] += mj * acc_visc_i * dv[0];
  pi->a_visc[1] += mj * acc_visc_i * dv[1];
  pi->a_visc[2] += mj * acc_visc_i * dv[2];

  pj->drho_dt += mi * dens_j;
  pj->a_hydro[0] -= mi * (acc_visc_j * dv[0] - acc_sph_j * dx[0]);
  pj->a_hydro[1] -= mi * (acc_visc_j * dv[1] - acc_sph_j * dx[1]);
  pj->a_hydro[2] -= mi * (acc_visc_j * dv[2] - acc_sph_j * dx[2]);

  pj->a_visc[0] -= mi * acc_visc_j * dv[0];
  pj->a_visc[1] -= mi * acc_visc_j * dv[1];
  pj->a_visc[2] -= mi * acc_visc_j * dv[2];


  boundary_fluid_interaction(pi, pj, r, r2, dx);

 // if(pi->id == 3374 && pj->is_boundary) printf("i: hydro = [%e %e %e, dv[0] = %e, dx[0] = %e, dx[1]=%e, r=%e P=%e, dens=%e, P/rho=%e, w_dx=%e\n", -mj * acc_sph_i * dx[0], -mj * acc_sph_i * dx[1], -mj * acc_sph_i * dx[2], dv[0], dx[0], dx[1], r, pj->pressure, mj*dens_i,(Pi_over_rhosq + Pj_over_rhosq),wi_dx );
//  if(pj->id == 3374 && pi->is_boundary) printf("j: hydro = [%e %e %e], dv[0] = %e, dx[0] = %e, dx[1]=%e, r=%e P=%e, dens=%e,  P/rho=%e, w_dx=%e\n", mi * acc_sph_j * dx[0], mi * acc_sph_j * dx[1], mi * acc_sph_j * dx[2], -dv[0], -dx[0], dx[1], r, pi->pressure, mi*dens_j, (Pi_over_rhosq + Pj_over_rhosq),wj_dx );  
//  if(pi->id == 3374) printf("viscosity a[0] = [ %e %e %e], hydro = [%e %e %e]\n", mj * acc_visc_i * dv[0],mj * acc_visc_i * dv[1], mj * acc_visc_i * dv[2], -mj * acc_sph_i * dx[0], -mj * acc_sph_i * dx[1], -mj * acc_sph_i * dx[2] );  
//  if(pj->id == 3374) printf("viscosity a[0] = [ %e %e %e], hydro = [%e %e %e]\n", -mi * acc_visc_j * dv[0], -mi * acc_visc_j * dv[1], -mi * acc_visc_j * dv[2], mi * acc_sph_j * dx[0], mi * acc_sph_j * dx[1], mi * acc_sph_j * dx[2] );  

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
  const int r_inv = 1.0 / r;

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

  /* Do the 'regular' hydro forces... */
  const float acc_sph_i = (Pi_over_rhosq + Pj_over_rhosq) * r_inv * wi_dx;

  /* Viscosity implementation */
  const float mu_i = rhoi * pi->dynamic_viscosity;
  const float mu_j = rhoj * pj->dynamic_viscosity;
  const float mu_sum = mu_i + mu_j;
  const float dens_ij_inv = rhoi_inv * rhoj_inv;
  const float acc_visc_i = mu_sum * wi_dx * dens_ij_inv * r_inv;
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];

  /* Apply forces */
  pi->drho_dt += mj * dens;
  pi->a_hydro[0] += mj * (acc_visc_i * dv[0] - acc_sph_i * dx[0]); 
  pi->a_hydro[1] += mj * (acc_visc_i * dv[1] - acc_sph_i * dx[1]);
  pi->a_hydro[2] += mj * (acc_visc_i * dv[2] - acc_sph_i * dx[2]);

  pi->a_visc[0] += mj * acc_visc_i * dv[0];
  pi->a_visc[1] += mj * acc_visc_i * dv[1];
  pi->a_visc[2] += mj * acc_visc_i * dv[2];

  boundary_fluid_interaction_nonsym(pi, pj, r, r2, dx);

//  if(pi->id == 3374 && pj->is_boundary) printf("i: hydro = [%e %e %e, dv[0] = %e, dx[0] = %e, dx[1]=%e, r=%e P=%e, dens=%e, P/rho=%e, w_dx=%e\n", -mj * acc_sph_i * dx[0], -mj * acc_sph_i * dx[1], -mj * acc_sph_i * dx[2], dv[0], dx[0], dx[1], r, pj->pressure, mj*dens,(Pi_over_rhosq + Pj_over_rhosq),wi_dx );
  ///*if(pi->id == 3374)*/ printf("viscosity a[0] = [ %e %e %e], hydro = [%e %e %e]\n", mj * acc_visc_i * dv[0],mj * acc_visc_i * dv[1], mj * acc_visc_i * dv[2], -mj * acc_sph_i * dx[0], -mj * acc_sph_i * dx[1], -mj * acc_sph_i * dx[2] );  
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
