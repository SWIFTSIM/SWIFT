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
#include "density_diffusion.h"
#include "boundary.h"
#include "particle_shifting.h"

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
  //float Pi_over_rhosq, Pj_over_rhosq;

/*  if( (pi->is_boundary && !pj->is_boundary) || (pj->is_boundary && !pi->is_boundary)){
    printf("We did a thing!\n");
  }*/

  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  rhoi_inv = 1.0f / rhoi;
  if(rhoi_inv != rhoi_inv){
    error("NaN rho = %e\n", rhoi);
  }
  rhoj_inv = 1.0f / rhoj;
  if(rhoj_inv != rhoj_inv){
    error("NaN rho = %e\n", rhoj);
  }
  const float r_inv = 1.0 / r;
  if(r_inv != r_inv){
    error("Zero distance: %e %lli %lli\n", r, pi->id, pj->id);
  }

  /* Compute dv */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  if(dv[0] != dv[0]){
    error("dv0 %e %e\n", pi->v[0], pj->v[0]);
  }
  if(dv[1] != dv[1]){
    error("dv1 %e %e\n", pi->v[1], pj->v[1]);
  }
  if(dv[2] != dv[2]){
    error("dv2 %e %e\n", pi->v[2], pj->v[2]);
  }
 
  inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
  inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);
  if(inv_hidim_pow_plus_one != inv_hidim_pow_plus_one){
    error("inv_hidim_pow_plus_one=%e\n", inv_hidim_pow_plus_one);
  }
  if(inv_hjdim_pow_plus_one != inv_hjdim_pow_plus_one){
    error("inv_hjdim_pow_plus_one=%e\n", inv_hjdim_pow_plus_one);
  }


//  Pi_over_rhosq = pi->pressure * (rhoi_inv * rhoi_inv);
//  Pj_over_rhosq = pj->pressure * (rhoj_inv * rhoj_inv);

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

  
  if(pi->a_hydro[0] != pi->a_hydro[0] || pi->a_hydro[1] != pi->a_hydro[1] || pi->a_hydro[2] != pi->a_hydro[2]){
    error("Unknown reason for NaN atm");
  }

  /* Compute the laminar viscosity term*/
  const float visc_i = 4.0 * pi->viscosity;
  const float visc_j = 4.0 * pj->viscosity;
  const float inv_dens_sum = rhoi + rhoj;
  const float r2_eta2 = r2 + 0.01*0.01;
  const float inv_r2eta2 = 1.0 / r2_eta2;
  const float temp_i = visc_i * inv_r2eta2 * inv_dens_sum;
  const float temp_j = visc_j * inv_r2eta2 * inv_dens_sum;
  const float multiplier_i = dx[0]*dx[0]*wi_dx + dx[1]*dx[1]*wi_dx + dx[2]*dx[2]*wi_dx;
  const float multiplier_j = dx[0]*dx[0]*wj_dx + dx[1]*dx[1]*wj_dx + dx[2]*dx[2]*wj_dx;
  //const float multiplier_i = dx[0]*wi_dx + dx[1]*wi_dx + dx[2]*wi_dx;
  //const float multiplier_j = dx[0]*wj_dx + dx[1]*wj_dx + dx[2]*wj_dx;
  const float prea0 = pi->a_hydro[0];
  const float prea1 = pi->a_hydro[1];
  const float prea2 = pi->a_hydro[2];
  pi->a_hydro[0] += mj*temp_i*multiplier_i*dv[0];
  pi->a_hydro[1] += mj*temp_i*multiplier_i*dv[1];
  pi->a_hydro[2] += mj*temp_i*multiplier_i*dv[2];

  if(pi->a_hydro[0] != pi->a_hydro[0] || pi->a_hydro[1] != pi->a_hydro[1] || pi->a_hydro[2] != pi->a_hydro[2]){
    error("NaN with %e %e %e %e %e %e %e %e", mj*temp_i*multiplier_i, dv[0], dv[1], dv[2], r2, prea0, prea1, prea2);
  }

  pj->a_hydro[0] -= mi*temp_j*multiplier_j*dv[0];
  pj->a_hydro[1] -= mi*temp_j*multiplier_j*dv[1];
  pj->a_hydro[2] -= mi*temp_j*multiplier_j*dv[2];
  if(pj->a_hydro[0] != pj->a_hydro[0] || pj->a_hydro[1] != pj->a_hydro[1] || pj->a_hydro[2] != pj->a_hydro[2]){
    error("NaN with %e %e %e %e", mi*temp_j*multiplier_j, dv[0], dv[1], dv[2]);
  }

  if(0 && /*!pi->is_boundary*/pi->id == 2531/*4129*/ && pj->is_boundary){
    printf("visc_i = %e laminar %e\n", visc_i, mj*temp_i*multiplier_i*dv[0]);
//    printf("%lli %e\n", pi->id, pj->v[0]);
  }
  
  /* SPS turbulence term */
  const float tau_xx = pi->tau_xx + pj->tau_xx;
  const float tau_xy = pi->tau_xy + pj->tau_xy;
  const float tau_xz = pi->tau_xz + pj->tau_xz;
  const float tau_yy = pi->tau_yy + pj->tau_yy;
  const float tau_yz = pi->tau_yz + pj->tau_yz;
  const float tau_zz = pi->tau_zz + pj->tau_zz;
  pi->a_hydro[0] += mi*mj*(tau_xx*dx[0] + tau_xy * dx[1] + tau_xz * dx[2]);
  pi->a_hydro[1] += mi*mj*(tau_xy*dx[0] + tau_yy * dx[1] + tau_yz * dx[2]);
  pi->a_hydro[2] += mi*mj*(tau_xz*dx[0] + tau_yz * dx[1] + tau_zz * dx[2]);
  if(0 &&/*!pi->is_boundary*/pi->id == 2531/*4129*/ && pj->is_boundary){
    printf("SPS %e\n", mi*mj*(tau_xx*dx[0] + tau_xy * dx[1] + tau_xz * dx[2]));
//    printf("%lli %e\n", pi->id, pj->v[0]);
  }
  if(pi->a_hydro[0] != pi->a_hydro[0] || pi->a_hydro[1] != pi->a_hydro[1] || pi->a_hydro[2] != pi->a_hydro[2]){
    error("NaN with %e %e %e %e %e %e", tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz );
  }

  pj->a_hydro[0] -= mi*mj*(tau_xx*dx[0] + tau_xy * dx[1] + tau_xz * dx[2]);
  pj->a_hydro[1] -= mi*mj*(tau_xy*dx[0] + tau_yy * dx[1] + tau_yz * dx[2]);
  pj->a_hydro[2] -= mi*mj*(tau_xz*dx[0] + tau_yz * dx[1] + tau_zz * dx[2]);


  if(pj->a_hydro[0] != pj->a_hydro[0] || pj->a_hydro[1] != pj->a_hydro[1] || pj->a_hydro[2] != pj->a_hydro[2]){
    error("NaN with %e %e %e %e %e %e", tau_xx, tau_xy, tau_xz, tau_yy, tau_yz, tau_zz );
  }


  /* Compute Velocity gradients */
  const float mj_over_rhoj = mj * rhoj_inv;
  const float mi_over_rhoi = mi * rhoi_inv;
  const float dv_over_rho_x_i = dv[0]*mj_over_rhoj;
  const float dv_over_rho_y_i = dv[1]*mj_over_rhoj;
  const float dv_over_rho_z_i = dv[2]*mj_over_rhoj;
  const float dv_over_rho_x_j = dv[0]*mi_over_rhoi;
  const float dv_over_rho_y_j = dv[1]*mi_over_rhoi;
  const float dv_over_rho_z_j = dv[2]*mi_over_rhoi;

  pi->grad_v_xx += dv_over_rho_x_i * dx[0];
  pi->grad_v_xy += dv_over_rho_x_i * dx[1];
  pi->grad_v_xz += dv_over_rho_x_i * dx[2];
  pi->grad_v_xy += dv_over_rho_y_i * dx[0];
  pi->grad_v_yy += dv_over_rho_y_i * dx[1];
  pi->grad_v_yz += dv_over_rho_y_i * dx[2];
  pi->grad_v_xz += dv_over_rho_z_i * dx[0];
  pi->grad_v_yz += dv_over_rho_z_i * dx[1];
  pi->grad_v_zz += dv_over_rho_z_i * dx[2];


  pj->grad_v_xx += dv_over_rho_x_j * dx[0];
  pj->grad_v_xy += dv_over_rho_x_j * dx[1];
  pj->grad_v_xz += dv_over_rho_x_j * dx[2];
  pj->grad_v_xy += dv_over_rho_y_j * dx[0];
  pj->grad_v_yy += dv_over_rho_y_j * dx[1];
  pj->grad_v_yz += dv_over_rho_y_j * dx[2];
  pj->grad_v_xz += dv_over_rho_z_j * dx[0];
  pj->grad_v_yz += dv_over_rho_z_j * dx[1];
  pj->grad_v_zz += dv_over_rho_z_j * dx[2];

  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
//  const float acc = (Pi_over_rhosq + Pj_over_rhosq) * r_inv;
  const float acc = (pi->pressure + pj->pressure) * rhoi_inv * rhoj_inv * r_inv;

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];

  pi->drho_dt += mj * dens;
  pi->a_hydro[0] -= mj * acc * wi_dx * dx[0];
  pi->a_hydro[1] -= mj * acc * wi_dx * dx[1];
  pi->a_hydro[2] -= mj * acc * wi_dx * dx[2];
   
  if(0 &&/*!pi->is_boundary*/pi->id == 2531/*4129*/ && pj->is_boundary){
    printf("hydro %e\n",  -mj * acc * wi_dx * dx[0]);
//    printf("%lli %e\n", pi->id, pj->v[0]);
  } 
/*  if(pi->id == 462){
    printf("hydro %e\n",  -mj * acc * wi_dx * dx[1]);
  } 
  if(pj->id == 462){
    printf("hydro %e\n",  mi * acc * wj_dx * dx[1]);
  } */
  if(pi->a_hydro[0] != pi->a_hydro[0] || pi->a_hydro[1] != pi->a_hydro[1] || pi->a_hydro[2] != pi->a_hydro[2]){
    error("NaN with %e %e %e %e %e %e %e %e %i %i %e %e %e %llu %llu", mj, acc, wi_dx, rhoi, rhoj, r_inv, r2, r, pi->is_boundary, pj->is_boundary, pi->x[0], pi->x[1], pi->x[2], pi->id, pj->id );
  }

  /* Compute density of pj. */
  pj->drho_dt +=  mi * dens;
  pj->a_hydro[0] += mi * acc * wj_dx * dx[0];
  pj->a_hydro[1] += mi * acc * wj_dx * dx[1];
  pj->a_hydro[2] += mi * acc * wj_dx * dx[2];

/*  if(pi->id == 4281){
    printf("a_hydro= %e %e %e %e\n", -mj * acc * wi_dx * dx[0], pj->x[0], pj->x[1], pj->x[2]);
  }
  if(pj->id == 4281){
    printf("a_hydro= %e %e %e %e\n", mi * acc * wj_dx * dx[0], pi->x[0], pi->x[1], pi->x[2]);
  }*/

  if(!pi->is_boundary/*pi->id == 4090*/ && pj->is_boundary){
    //printf("%e %e %e\n", pj->v[0], pj->v[1], pj->v[2]);
    //printf("%e %e %e\n", dv_over_rho_x_i, mj_over_rhoj, dv[0]);
//    printf("%lli %e\n", pi->id, pj->v[0]);
   // printf("%lli %e\n", pi->id, r2);
  }  

  /* Store viscous effect towards CFL condition */
  const float dvdr = dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2];
  const float dvdr_rr2 = dvdr * inv_r2eta2;
  pi->max_visc = fmaxf(pi->max_visc, dvdr_rr2);
  pj->max_visc = fmaxf(pj->max_visc, dvdr_rr2);

  boundary_fluid_interaction(pi, pj, r, r2, dx);
  compute_density_diffusive_term(pi, pj, r2, r, wi_dx, wj_dx, dx);
  compute_shifting_term(pi, pj, r2, r, wi_dx, wj_dx, dx);
  if(pi->a_hydro[0] != pi->a_hydro[0] || pi->a_hydro[1] != pi->a_hydro[1] || pi->a_hydro[2] != pi->a_hydro[2]){
    error("NaN with ?!");
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
  float mi,mj;
  float inv_hidim_pow_plus_one;
  //float inv_hjdim_pow_plus_one;
  float dv[3];
  float rhoi, rhoj, rhoi_inv, rhoj_inv;
//  float Pi_over_rhosq, Pj_over_rhosq;

  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  rhoi_inv = 1.0f / rhoi;
  rhoj_inv = 1.0f / rhoj;
  const float r_inv = 1.0 / r;

  /* Compute dv */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];

  inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
 // inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);


//  Pi_over_rhosq = pi->pressure * (rhoi_inv * rhoi_inv);
//  Pj_over_rhosq = pj->pressure * (rhoj_inv * rhoj_inv);

  /* Compute density of pi. */
  h_inv = 1.0f / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dx = wi_dx * inv_hidim_pow_plus_one;

  /* Compute the laminar viscosity term*/
  const float visc_i = 4.0 * pi->viscosity;
  const float inv_dens_sum = rhoi + rhoj;
  const float r2_eta2 = r2 + 0.01*0.01;
  const float inv_r2eta2 = 1.0 / r2_eta2;
  const float temp_i = visc_i * inv_r2eta2 * inv_dens_sum;
  const float multiplier_i = dx[0]*dx[0]*wi_dx + dx[1]*dx[1]*wi_dx + dx[2]*dx[2]*wi_dx;
  pi->a_hydro[0] += mj*temp_i*multiplier_i*dv[0];
  pi->a_hydro[1] += mj*temp_i*multiplier_i*dv[1];
  pi->a_hydro[2] += mj*temp_i*multiplier_i*dv[2];

  /* SPS turbulence term */
  const float tau_xx = pi->tau_xx + pj->tau_xx;
  const float tau_xy = pi->tau_xy + pj->tau_xy;
  const float tau_xz = pi->tau_xz + pj->tau_xz;
  const float tau_yy = pi->tau_yy + pj->tau_yy;
  const float tau_yz = pi->tau_yz + pj->tau_yz;
  const float tau_zz = pi->tau_zz + pj->tau_zz;
  pi->a_hydro[0] += mi*mj*(tau_xx*dx[0] + tau_xy * dx[1] + tau_xz * dx[2]);
  pi->a_hydro[1] += mi*mj*(tau_xy*dx[0] + tau_yy * dx[1] + tau_yz * dx[2]);
  pi->a_hydro[2] += mi*mj*(tau_xz*dx[0] + tau_yz * dx[1] + tau_zz * dx[2]);

  /* Compute Velocity gradients */
  const float mj_over_rhoj = mj * rhoj_inv;
  const float dv_over_rho_x_i = dv[0]*mj_over_rhoj;
  const float dv_over_rho_y_i = dv[1]*mj_over_rhoj;
  const float dv_over_rho_z_i = dv[2]*mj_over_rhoj;

  pi->grad_v_xx += dv_over_rho_x_i * dx[0];
  pi->grad_v_xy += dv_over_rho_x_i * dx[1];
  pi->grad_v_xz += dv_over_rho_x_i * dx[2];
  pi->grad_v_xy += dv_over_rho_y_i * dx[0];
  pi->grad_v_yy += dv_over_rho_y_i * dx[1];
  pi->grad_v_yz += dv_over_rho_y_i * dx[2];
  pi->grad_v_xz += dv_over_rho_z_i * dx[0];
  pi->grad_v_yz += dv_over_rho_z_i * dx[1];
  pi->grad_v_zz += dv_over_rho_z_i * dx[2];


  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
//  const float acc = (Pi_over_rhosq + Pj_over_rhosq) * r_inv;
  const float acc = (pi->pressure + pj->pressure) * rhoi_inv * rhoj_inv * r_inv;

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];
//  printf("dens=%f dv[0]=%f dx[0]=%f \n", dens, dv[0], dx[0]);

  pi->drho_dt += mj * dens;
  pi->a_hydro[0] -= mj * acc * wi_dx * dx[0];
  pi->a_hydro[1] -= mj * acc * wi_dx * dx[1];
  pi->a_hydro[2] -= mj * acc * wi_dx * dx[2];

  if(pi->id == 4281){
    printf("a_hydro= %e\n", -mj * acc * wi_dx * dx[0]);
  }
  /* Store viscous effect towards CFL condition */
  const float dvdr = dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2];
  const float dvdr_rr2 = dvdr * inv_r2eta2;
  pi->max_visc = fmaxf(pi->max_visc, dvdr_rr2);

  boundary_fluid_interaction_nonsym(pi, pj, r, r2, dx);
  compute_density_diffusive_term_asym(pi, pj, r2, r, wi_dx, dx);
  compute_shifting_term_nonsym(pi, pj, r2, r, wi_dx, dx);

}
//#endif

#endif /* SWIFT_MINIMAL_HYDRO_IACT_H */
