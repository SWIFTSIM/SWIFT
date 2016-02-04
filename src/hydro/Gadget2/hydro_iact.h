/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_IACT_LEGACY_H
#define SWIFT_RUNNER_IACT_LEGACY_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "part.h"
#include "vector.h"


/**
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 *and use the same
 * numerical coefficients as the Gadget-2 code. When used with the Spline-3
 *kernel, the results
 * should be equivalent to the ones obtained with Gadget-2 up to the rounding
 *errors and interactions
 * missed by the Gadget-2 tree-code neighbours search.
 *
 * The code uses internal energy instead of entropy as a thermodynamical
 *variable.
 */

/**
 * @brief Density loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float wj, wj_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mi = pj->mass;
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Compute the kernel function for pi */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->rho_dh -= mj * kernel_igamma * (3.f * wi + ui * wi_dx);
  
  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= ui * wi_dx;


  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  pj->rho_dh -= mi * kernel_igamma * (3.f * wj + uj * wj_dx);
  
  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= uj * wj_dx;

  
  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;
  
  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->div_v -= faci * dvdr;
  pj->div_v -= facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->rot_v[0] += faci * curlvr[0];
  pi->rot_v[1] += faci * curlvr[1];
  pi->rot_v[2] += faci * curlvr[2];

  pj->rot_v[0] += facj * curlvr[0];
  pj->rot_v[1] += facj * curlvr[1];
  pj->rot_v[2] += facj * curlvr[2];

  /* if(pi->id == 1000) */
  /*   message("Interacting with %lld. r=%f\n", pj->id, r); */

  /* if(pj->id == 1000) */
  /*   message("Interacting with %lld. r=%f\n", pi->id, r); */

}


/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  const float ri = 1.0f / r;

  if(pi->id == 1000 && pj->id == 1103) hi = 0.2234976 / 2.f;
  
  /* Compute the kernel function */
  const float h_inv = 1.0f / hi;
  const float u = r * h_inv;
  kernel_deval(u, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->rho_dh -= mj * kernel_igamma * (3.f * wi + u * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= u * wi_dx;

  const float ih3 = h_inv * h_inv * h_inv;
  const float ih4 = h_inv * h_inv * h_inv * h_inv;
  
  if(pi->id == 1000 && pj->id == 1103)
    message("Interacting with %lld. r=%e hi=%e u=%e W=%e dW/dx=%e dh_drho1=%e dh_drho2=%e %e\n",
	    pj->id,
	    r,
	    hi,
	    u,
	    wi * ih3,
	    wi_dx * ih4,
	    -mj * (3.f * kernel_igamma * wi) * ih4,
	    -mj * u * wi_dx * kernel_igamma * ih4,
	    kernel_igamma
	    );

  const float fac = mj * wi_dx * ri;
  
  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  pi->div_v -= fac * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->rot_v[0] += fac * curlvr[0];
  pi->rot_v[1] += fac * curlvr[1];
  pi->rot_v[2] += fac * curlvr[2];  
}


/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* float r = sqrtf(r2), ri = 1.0f / r; */
  /* float xi, xj; */
  /* float hi_inv, hi2_inv; */
  /* float hj_inv, hj2_inv; */
  /* float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr; */
  /* float mi, mj, POrho2i, POrho2j, rhoi, rhoj; */
  /* float v_sig, omega_ij, Pi_ij; */
  /* float f; */
  /* int k; */

  /* /\* Get some values in local variables. *\/ */
  /* mi = pi->mass; */
  /* mj = pj->mass; */
  /* rhoi = pi->rho; */
  /* rhoj = pj->rho; */
  /* //POrho2i = pi->force.POrho2; */
  /* //POrho2j = pj->force.POrho2; */

  /* /\* Get the kernel for hi. *\/ */
  /* hi_inv = 1.0f / hi; */
  /* hi2_inv = hi_inv * hi_inv; */
  /* xi = r * hi_inv; */
  /* kernel_deval(xi, &wi, &wi_dx); */
  /* wi_dr = hi2_inv * hi2_inv * wi_dx; */

  /* /\* Get the kernel for hj. *\/ */
  /* hj_inv = 1.0f / hj; */
  /* hj2_inv = hj_inv * hj_inv; */
  /* xj = r * hj_inv; */
  /* kernel_deval(xj, &wj, &wj_dx); */
  /* wj_dr = hj2_inv * hj2_inv * wj_dx; */

  /* /\* Compute dv dot r. *\/ */
  /* dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] + */
  /*        (pi->v[2] - pj->v[2]) * dx[2]; */
  /* dvdr *= ri; */

  /* /\* Compute the relative velocity. (This is 0 if the particles move away from */
  /*  * each other and negative otherwise) *\/ */
  /* omega_ij = fminf(dvdr, 0.f); */

  /* /\* Compute signal velocity *\/ */
  /* v_sig = pi->force.c + pj->force.c - 3.0f * omega_ij; */

  /* /\* Compute viscosity tensor *\/ */
  /* Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (rhoi + rhoj); */

  /* /\* Apply balsara switch *\/ */
  /* Pi_ij *= (pi->force.balsara + pj->force.balsara); */

  /* /\* Get the common factor out. *\/ */
  /* w = ri * */
  /*     ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr)); */

  /* /\* Use the force, Luke! *\/ */
  /* for (k = 0; k < 3; k++) { */
  /*   f = dx[k] * w; */
  /*   pi->a[k] -= mj * f; */
  /*   pj->a[k] += mi * f; */
  /* } */

  /* /\* Get the time derivative for u. *\/ */
  /* pi->force.u_dt += */
  /*     mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr)); */
  /* pj->force.u_dt += */
  /*     mi * dvdr * (POrho2j * wj_dr + 0.125f * Pi_ij * (wi_dr + wj_dr)); */

  /* /\* Get the time derivative for h. *\/ */
  /* pi->force.h_dt -= mj * dvdr / rhoj * wi_dr; */
  /* pj->force.h_dt -= mi * dvdr / rhoi * wj_dr; */

  /* /\* Update the signal velocity. *\/ */
  /* pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig); */
  /* pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig); */
}


/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  /* float r = sqrtf(r2), ri = 1.0f / r; */
  /* float xi, xj; */
  /* float hi_inv, hi2_inv; */
  /* float hj_inv, hj2_inv; */
  /* float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr; */
  /* float /\*mi,*\/ mj, POrho2i, POrho2j, rhoi, rhoj; */
  /* float v_sig, omega_ij, Pi_ij; */
  /* float f; */
  /* int k; */

  /* /\* Get some values in local variables. *\/ */
  /* // mi = pi->mass; */
  /* mj = pj->mass; */
  /* rhoi = pi->rho; */
  /* rhoj = pj->rho; */
  /* POrho2i = pi->force.POrho2; */
  /* POrho2j = pj->force.POrho2; */

  /* /\* Get the kernel for hi. *\/ */
  /* hi_inv = 1.0f / hi; */
  /* hi2_inv = hi_inv * hi_inv; */
  /* xi = r * hi_inv; */
  /* kernel_deval(xi, &wi, &wi_dx); */
  /* wi_dr = hi2_inv * hi2_inv * wi_dx; */

  /* /\* Get the kernel for hj. *\/ */
  /* hj_inv = 1.0f / hj; */
  /* hj2_inv = hj_inv * hj_inv; */
  /* xj = r * hj_inv; */
  /* kernel_deval(xj, &wj, &wj_dx); */
  /* wj_dr = hj2_inv * hj2_inv * wj_dx; */

  /* /\* Compute dv dot r. *\/ */
  /* dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] + */
  /*        (pi->v[2] - pj->v[2]) * dx[2]; */
  /* dvdr *= ri; */

  /* /\* Compute the relative velocity. (This is 0 if the particles move away from */
  /*  * each other and negative otherwise) *\/ */
  /* omega_ij = fminf(dvdr, 0.f); */

  /* /\* Compute signal velocity *\/ */
  /* v_sig = pi->force.c + pj->force.c - 3.0f * omega_ij; */

  /* /\* Compute viscosity tensor *\/ */
  /* Pi_ij = -const_viscosity_alpha * v_sig * omega_ij / (rhoi + rhoj); */

  /* /\* Apply balsara switch *\/ */
  /* Pi_ij *= (pi->force.balsara + pj->force.balsara); */

  /* /\* Get the common factor out. *\/ */
  /* w = ri * */
  /*     ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr)); */

  /* /\* Use the force, Luke! *\/ */
  /* for (k = 0; k < 3; k++) { */
  /*   f = dx[k] * w; */
  /*   pi->a[k] -= mj * f; */
  /* } */

  /* /\* Get the time derivative for u. *\/ */
  /* pi->force.u_dt += */
  /*     mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr)); */

  /* /\* Get the time derivative for h. *\/ */
  /* pi->force.h_dt -= mj * dvdr / rhoj * wi_dr; */

  /* /\* Update the signal velocity. *\/ */
  /* pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig); */
  /* pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig); */
}


#endif /* SWIFT_RUNNER_IACT_LEGACY_H */
