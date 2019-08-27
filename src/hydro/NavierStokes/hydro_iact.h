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
#ifndef SWIFT_NAVIER_STOK_HYDRO_IACT_H
#define SWIFT_NAVIER_STOK_HYDRO_IACT_H

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

    float r = sqrtf(r2);
    float wi, wi_dx, wj, wj_dx;
    float hi_inv, hj_inv;
    float xi, xj;
    float rhoi, rhoj, rhoi_inv, rhoj_inv;
    float inv_hidim_pow_plus_one;
    float inv_hjdim_pow_plus_one;
    float dv[3];
    float mi, mj;


    pi->neighbours++;
    pj->neighbours++;
    mi = pi->mass;
    mj = pj->mass;
    rhoi = pi->rho;
    rhoj = pj->rho;
    rhoi_inv = 1.0f / rhoi;
    rhoj_inv = 1.0f / rhoj;

    dv[0] = pi->v[0] - pj->v[0];
    dv[1] = pi->v[1] - pj->v[1];
    dv[2] = pi->v[2] - pj->v[2];

    inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
    inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);

    hi_inv = 1.0f / hi;
    xi = r * hi_inv;
    kernel_deval(xi, &wi, &wi_dx);
    wi_dx = wi_dx * inv_hidim_pow_plus_one;

    hj_inv = 1.0 / hj;
    xj = r * hj_inv;
    kernel_deval(xj, &wj, &wj_dx);
    wj_dx = wj_dx * inv_hjdim_pow_plus_one;

   const float dvdr = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];


    pi->div_v += rhoi_inv * wi_dx * mj * ( dvdr );
    pj->div_v += rhoj_inv * wj_dx * mi * ( dvdr );

    const float faci = rhoi_inv * wi_dx * mj;
    const float facj = rhoj_inv * wj_dx * mi;
    pi->dvx_xx += faci * ( dv[0] * dx[0] ); //*dx[0] i think.
    pj->dvx_xx += facj * ( dv[0] * dx[0] );
//    float test;
 //   kernel_deval(xi, &wi, &test);
    pi->dvx_xy += faci * ( dv[0] * dx[1]);
    pj->dvx_xy += facj * ( dv[0] * dx[1]);

  //  if((pi->id == 6144 || pi->id == 7008) && dv[0] != 0.0) printf("%llu,%e,%e\n", pi->id, faci * (dv[0] * dx[0]),pi->dvx_xx);
 //   if((pj->id == 6144 || pj->id == 7008) && dv[0] != 0.0) printf("%llu,%e,%e\n", pj->id, facj * (dv[0] * dx[0]),pj->dvx_xx);
    pi->dvx_xz += faci * ( dv[0] * dx[2]);
    pj->dvx_xz += facj * ( dv[0] * dx[2]);

    pi->dvy_xx += faci * ( dv[1] * dx[0]);
    pj->dvy_xx += facj * ( dv[1] * dx[0]);

    pi->dvy_xy += faci * ( dv[1] * dx[1]);
    pj->dvy_xy += facj * ( dv[1] * dx[1]);
   
    pi->dvy_xz += faci * ( dv[1] * dx[2]);
    pj->dvy_xz += facj * ( dv[1] * dx[2]);
 
    pi->dvz_xx += faci * ( dv[2] * dx[0]);
    pj->dvz_xx += facj * ( dv[2] * dx[0]);

    pi->dvz_xy += faci * ( dv[2] * dx[1]);
    pj->dvz_xy += facj * ( dv[2] * dx[1]);
 
    pi->dvz_xz += faci * ( dv[2] * dx[2]);
    pj->dvz_xz += facj * ( dv[2] * dx[2]);

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

    float r = sqrtf(r2);
    float wi, wi_dx;
    float hi_inv;
    float xi;
    float rhoi, rhoi_inv;
    float inv_hidim_pow_plus_one;
    float dv[3];
    const float mj = pj->mass;

    pi->neighbours++;
    rhoi = pi->rho;
    rhoi_inv = 1.0f / rhoi;

    dv[0] = pi->v[0] - pj->v[0];
    dv[1] = pi->v[1] - pj->v[1];
    dv[2] = pi->v[2] - pj->v[2];

    inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);

    hi_inv = 1.0f / hi;
    xi = r * hi_inv;
    kernel_deval(xi, &wi, &wi_dx);
    wi_dx = wi_dx * inv_hidim_pow_plus_one;

    //const float dvdr = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];
    const float dvdr = dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

    pi->div_v += rhoi_inv * wi_dx * mj * ( dvdr );

//    if(pi->id == 57600 && pi->is_boundary && !pj->is_boundary) printf("neighbour = %llu\n", pj->id);
    const float faci = rhoi_inv * wi_dx * mj;
    float test;
    kernel_deval(xi, &wi, &test);
//    if(pi->id == 115200 || pi->id == 153120 ) printf("%llu %llu %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",pi->id, pj->id, faci * ( dv[0] * dx[0]), faci * dx[1] * dv[1], faci * dx[2] * dv[2], dx[0], dx[1], dx[2],dv[0], r, pj->x[0], pj->x[1], pj->x[2], test, wi, dvdr);
    pi->dvx_xx += faci * ( dv[0] * dx[0] );
//    if((pi->id == 6144 || pi->id == 7008) && dv[0] != 0.0) printf("%llu,%e,%e\n", pi->id, faci * (dv[0] * dx[0]),pi->dvx_xx);
    pi->dvx_xy += faci * ( dv[0] * dx[1] );
    pi->dvx_xz += faci * ( dv[0] * dx[2] );
    pi->dvy_xx += faci * ( dv[1] * dx[0] );
    pi->dvy_xy += faci * ( dv[1] * dx[1] );
    pi->dvy_xz += faci * ( dv[1] * dx[2] );
    pi->dvz_xx += faci * ( dv[2] * dx[0] );
    pi->dvz_xy += faci * ( dv[2] * dx[1] );
    pi->dvz_xz += faci * ( dv[2] * dx[2] );

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
  float rhoi, rhoj, rhoi_inv, rhoj_inv, rhoi_inv_sq,rhoj_inv_sq;
  float Pi_over_rhosq, Pj_over_rhosq;

/*  if( (pi->is_boundary && !pj->is_boundary) || (pj->is_boundary && !pi->is_boundary)){
    printf("We did a thing!\n");
  }*/

  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  rhoi_inv = 1.0f / rhoi;
  rhoj_inv = 1.0f / rhoj;
  rhoi_inv_sq = rhoi_inv * rhoi_inv;
  rhoj_inv_sq = rhoj_inv * rhoj_inv;

  /* Compute dv */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
 
  inv_hidim_pow_plus_one = 1.0f / pow_dimension_plus_one(hi);
  inv_hjdim_pow_plus_one = 1.0f / pow_dimension_plus_one(hj);


  Pi_over_rhosq = pi->pressure * (rhoi_inv_sq);
  Pj_over_rhosq = pj->pressure * (rhoj_inv_sq);

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

  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
  const float acc = (Pi_over_rhosq + Pj_over_rhosq);

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];

/*  if( (!pj->is_boundary) || (!pi->is_boundary)){
    printf("r2 = %f, pressure = %f, xi=%f, rho=%f, acc=%f, wi_dx=%f, dx[1]=%f drho_dt=%f\n",r2,pi->pressure,xi, pi->rho, acc, wi_dx, dx[1], mj*dens);
    printf("hydro_i = %f, hydro_j = %f\n", -mj*acc*wi_dx*dx[1], mi*acc*wj_dx*dx[1]);
}*/
  const double taui_over_rhosq_xx = pi->tau_xx * rhoi_inv_sq;
  const double tauj_over_rhosq_xx = pj->tau_xx * rhoj_inv_sq;
  const double taui_over_rhosq_xy = pi->tau_xy * rhoi_inv_sq;
  const double tauj_over_rhosq_xy = pj->tau_xy * rhoj_inv_sq;
  const double taui_over_rhosq_xz = pi->tau_xz * rhoi_inv_sq;
  const double tauj_over_rhosq_xz = pj->tau_xz * rhoj_inv_sq;
  const double taui_over_rhosq_yy = pi->tau_yy * rhoi_inv_sq;
  const double tauj_over_rhosq_yy = pj->tau_yy * rhoj_inv_sq;
  const double taui_over_rhosq_yz = pi->tau_yz * rhoi_inv_sq;
  const double tauj_over_rhosq_yz = pj->tau_yz * rhoj_inv_sq;
  const double taui_over_rhosq_zz = pi->tau_zz * rhoi_inv_sq;
  const double tauj_over_rhosq_zz = pj->tau_zz * rhoj_inv_sq;

//  printf("dens=%f, dv[0]=%f, dx[0]=%f\n", dens,dv[0],dx[0]);
//    if(pi->id == 115440 && pj->is_boundary) printf("xx = %e, xy = %e, xz = %e\n", taui_over_rhosq_xx, taui_over_rhosq_xy, taui_over_rhosq_xz);
//    if(pj->id == 115440 && pi->is_boundary) printf("xx = %e, xy = %e, xz = %e\n", tauj_over_rhosq_xx, tauj_over_rhosq_xy, tauj_over_rhosq_xz);

  const float i_x_visc_acc = mj *( (wi_dx * dx[0] * (taui_over_rhosq_xx + tauj_over_rhosq_xx)) + (wi_dx * dx[1] * (taui_over_rhosq_xy + tauj_over_rhosq_xy) ) + (wi_dx * dx[2] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)));
  const double i_y_visc_acc = mj *( ( wi_dx * dx[0] * (taui_over_rhosq_xy + tauj_over_rhosq_xy)) + (wi_dx * dx[1] * (taui_over_rhosq_yy + tauj_over_rhosq_yy)) + (wi_dx * dx[2] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) );
  const float i_z_visc_acc = mj *( ( wi_dx * dx[0] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)) + (wi_dx * dx[1] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) + (wi_dx * dx[2]  * (taui_over_rhosq_zz + tauj_over_rhosq_zz)) );
  const float i_x_hydro_acc = mj * acc * wi_dx * dx[0];
  const float i_y_hydro_acc = mj * acc * wi_dx * dx[1];
  const float i_z_hydro_acc = mj * acc * wi_dx * dx[2];

  pi->drho_dt += mj * dens;
  pi->a_viscosity[0] -= i_x_visc_acc;
  pi->a_viscosity[1] -= i_y_visc_acc;
  if(pi->id == 6144) printf("%e,%e,%e,%e,%e\n", i_y_visc_acc, pi->a_viscosity[1],dx[0],dx[1],dx[2]);
  pi->a_viscosity[2] -= i_z_visc_acc;
//  pi->a_hydro[0] += i_x_visc_acc - i_x_hydro_acc;
 // pi->a_hydro[1] += i_y_visc_acc - i_y_hydro_acc;
  //pi->a_hydro[2] += i_z_visc_acc - i_z_hydro_acc;
  pi->a_hydro[0] -= i_x_hydro_acc;
  pi->a_hydro[1] -= i_y_hydro_acc;
  pi->a_hydro[2] -= i_z_hydro_acc;

  const float j_x_visc_acc = mi *( (wj_dx * -dx[0] * (taui_over_rhosq_xx + tauj_over_rhosq_xx)) + (wi_dx * -dx[1] * (taui_over_rhosq_xy + tauj_over_rhosq_xy) ) + (wi_dx * -dx[2] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)));
  const double j_y_visc_acc = mi *( ( wj_dx * -dx[0] * (taui_over_rhosq_xy + tauj_over_rhosq_xy)) + (wi_dx * -dx[1] * (taui_over_rhosq_yy + tauj_over_rhosq_yy)) + (wi_dx * -dx[2] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) );
  const float j_z_visc_acc = mi *( ( wj_dx * -dx[0] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)) + (wi_dx * -dx[1] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) + (wi_dx * -dx[2]  * (taui_over_rhosq_zz + tauj_over_rhosq_zz)) );
  const float j_x_hydro_acc = mi * acc * wj_dx * -dx[0];
  const float j_y_hydro_acc = mi * acc * wj_dx * -dx[1];
  const float j_z_hydro_acc = mi * acc * wj_dx * -dx[2];
  /* Compute density of pj. */
  pj->drho_dt +=  mi * dens;
  pj->a_viscosity[0] -= j_x_visc_acc;
  pj->a_viscosity[1] -= j_y_visc_acc;
  if(pj->id == 6144) printf("%e,%e,%e,%e,%e\n", j_y_visc_acc, pj->a_viscosity[1],dx[0],dx[1],dx[2]);
  pj->a_viscosity[2] -= j_z_visc_acc;
//  pj->a_hydro[0] += j_x_visc_acc - j_x_hydro_acc;
//  pj->a_hydro[1] += j_y_visc_acc - j_y_hydro_acc;
//  pj->a_hydro[2] += j_z_visc_acc - j_z_hydro_acc;
  pj->a_hydro[0] -= j_x_hydro_acc;
  pj->a_hydro[1] -= j_y_hydro_acc;
  pj->a_hydro[2] -= j_z_hydro_acc;
  
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
  float rhoi, rhoj, rhoi_inv, rhoj_inv, rhoi_inv_sq, rhoj_inv_sq;
  float Pi_over_rhosq, Pj_over_rhosq;

  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  rhoi_inv = 1.0f / rhoi;
  rhoj_inv = 1.0f / rhoj;
  rhoi_inv_sq = rhoi_inv * rhoi_inv;
  rhoj_inv_sq = rhoj_inv * rhoj_inv;

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

  /* Acceleration term (needs multiplying by the dx in relevant dimension to correctly
 *   calculate the grad W term. */
  const float acc = (Pi_over_rhosq + Pj_over_rhosq);

  /* drho_dt term exclusing particle masses */
  const float dens = dv[0] * wi_dx*dx[0] + dv[1] * wi_dx*dx[1] + dv[2] * wi_dx*dx[2];
//  printf("dens=%f dv[0]=%f dx[0]=%f \n", dens, dv[0], dx[0]);

  const double taui_over_rhosq_xx = pi->tau_xx * rhoi_inv_sq;
  const double tauj_over_rhosq_xx = pj->tau_xx * rhoj_inv_sq;
  const double taui_over_rhosq_xy = pi->tau_xy * rhoi_inv_sq;
  const double tauj_over_rhosq_xy = pj->tau_xy * rhoj_inv_sq;
  const double taui_over_rhosq_xz = pi->tau_xz * rhoi_inv_sq;
  const double tauj_over_rhosq_xz = pj->tau_xz * rhoj_inv_sq;
  const double taui_over_rhosq_yy = pi->tau_yy * rhoi_inv_sq;
  const double tauj_over_rhosq_yy = pj->tau_yy * rhoj_inv_sq;
  const double taui_over_rhosq_yz = pi->tau_yz * rhoi_inv_sq;
  const double tauj_over_rhosq_yz = pj->tau_yz * rhoj_inv_sq;
  const double taui_over_rhosq_zz = pi->tau_zz * rhoi_inv_sq;
  const double tauj_over_rhosq_zz = pj->tau_zz * rhoj_inv_sq;


  const float i_x_visc_acc = mj *( (wi_dx * dx[0] * (taui_over_rhosq_xx + tauj_over_rhosq_xx)) + (wi_dx * dx[1] * (taui_over_rhosq_xy + tauj_over_rhosq_xy) ) + (wi_dx * dx[2] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)));
  const double i_y_visc_acc = mj *( ( wi_dx * dx[0] * (taui_over_rhosq_xy + tauj_over_rhosq_xy)) + (wi_dx * dx[1] * (taui_over_rhosq_yy + tauj_over_rhosq_yy)) + (wi_dx * dx[2] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) );
  const float i_z_visc_acc = mj *( ( wi_dx * dx[0] * (taui_over_rhosq_xz + tauj_over_rhosq_xz)) + (wi_dx * dx[1] * (taui_over_rhosq_yz + tauj_over_rhosq_yz)) + (wi_dx * dx[2]  * (taui_over_rhosq_zz + tauj_over_rhosq_zz)) );
  const float i_x_hydro_acc = mj * acc * wi_dx * dx[0];
  const float i_y_hydro_acc = mj * acc * wi_dx * dx[1];
  const float i_z_hydro_acc = mj * acc * wi_dx * dx[2];
  pi->drho_dt += mj * dens;
  pi->a_viscosity[0] -= i_x_visc_acc;
  pi->a_viscosity[1] -= i_y_visc_acc;
  if(pi->id == 6144) printf("%e,%e,%e,%e,%e\n", i_y_visc_acc, pi->a_viscosity[1],dx[0],dx[1],dx[2]);
  pi->a_viscosity[2] -= i_z_visc_acc;
//  pi->a_hydro[0] += i_x_visc_acc - i_x_hydro_acc;
 // pi->a_hydro[1] += i_y_visc_acc - i_y_hydro_acc;
  //pi->a_hydro[2] += i_z_visc_acc - i_z_hydro_acc;
  pi->a_hydro[0] -= i_x_hydro_acc;
  pi->a_hydro[1] -= i_y_hydro_acc;
  pi->a_hydro[2] -= i_z_hydro_acc;
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
