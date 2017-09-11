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
#ifndef SWIFT_GADGET2_HYDRO_IACT_H
#define SWIFT_GADGET2_HYDRO_IACT_H

/**
 * @file Gadget2/hydro_iact.h
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 * Springel, V., MNRAS, Volume 364, Issue 4, pp. 1105-1134.
 * We use the same numerical coefficients as the Gadget-2 code. When used with
 * the Spline-3 kernel, the results should be equivalent to the ones obtained
 * with Gadget-2 up to the rounding errors and interactions missed by the
 * Gadget-2 tree-code neighbours search.
 */

#include "cache.h"
#include "minmax.h"

/**
 * @brief Density loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float wj, wj_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mi = pi->mass;
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
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  pj->density.rho_dh -= mi * (hydro_dimension * wj + uj * wj_dx);

  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= (hydro_dimension * wj + uj * wj_dx);

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
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
  const float r_inv = 1.0f / r;

  /* Compute the kernel function */
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->density.rho_dh -= mj * (hydro_dimension * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  const float fac = mj * wi_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  pi->density.div_v -= fac * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->density.rot_v[0] += fac * curlvr[0];
  pi->density.rot_v[1] += fac * curlvr[1];
  pi->density.rot_v[2] += fac * curlvr[2];
}

#ifdef WITH_VECTORIZATION

/**
 * @brief Density interaction computed using 1 vector
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_1_vec_density(vector *r2, vector *dx, vector *dy, vector *dz,
                                 vector hi_inv, vector vix, vector viy,
                                 vector viz, float *Vjx, float *Vjy, float *Vjz,
                                 float *Mj, vector *rhoSum, vector *rho_dhSum,
                                 vector *wcountSum, vector *wcount_dhSum,
                                 vector *div_vSum, vector *curlvxSum,
                                 vector *curlvySum, vector *curlvzSum,
                                 mask_t mask) {

  vector r, ri, ui, wi, wi_dx;
  vector mj;
  vector dvx, dvy, dvz;
  vector vjx, vjy, vjz;
  vector dvdr;
  vector curlvrx, curlvry, curlvrz;

  /* Fill the vectors. */
  mj.v = vec_load(Mj);
  vjx.v = vec_load(Vjx);
  vjy.v = vec_load(Vjy);
  vjz.v = vec_load(Vjz);

  /* Get the radius and inverse radius. */
  ri = vec_reciprocal_sqrt(*r2);
  r.v = vec_mul(r2->v, ri.v);

  ui.v = vec_mul(r.v, hi_inv.v);

  /* Calculate the kernel for two particles. */
  kernel_deval_1_vec(&ui, &wi, &wi_dx);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvz.v = vec_sub(viz.v, vjz.v);

  /* Compute dv dot r */
  dvdr.v = vec_fma(dvx.v, dx->v, vec_fma(dvy.v, dy->v, vec_mul(dvz.v, dz->v)));
  dvdr.v = vec_mul(dvdr.v, ri.v);

  /* Compute dv cross r */
  curlvrx.v =
      vec_fma(dvy.v, dz->v, vec_mul(vec_set1(-1.0f), vec_mul(dvz.v, dy->v)));
  curlvry.v =
      vec_fma(dvz.v, dx->v, vec_mul(vec_set1(-1.0f), vec_mul(dvx.v, dz->v)));
  curlvrz.v =
      vec_fma(dvx.v, dy->v, vec_mul(vec_set1(-1.0f), vec_mul(dvy.v, dx->v)));
  curlvrx.v = vec_mul(curlvrx.v, ri.v);
  curlvry.v = vec_mul(curlvry.v, ri.v);
  curlvrz.v = vec_mul(curlvrz.v, ri.v);

  vector wcount_dh_update;
  wcount_dh_update.v =
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v));

  /* Mask updates to intermediate vector sums for particle pi. */
  rhoSum->v = vec_mask_add(rhoSum->v, vec_mul(mj.v, wi.v), mask);
  rho_dhSum->v =
      vec_mask_sub(rho_dhSum->v, vec_mul(mj.v, wcount_dh_update.v), mask);
  wcountSum->v = vec_mask_add(wcountSum->v, wi.v, mask);
  wcount_dhSum->v = vec_mask_sub(wcount_dhSum->v, wcount_dh_update.v, mask);
  div_vSum->v =
      vec_mask_sub(div_vSum->v, vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)), mask);
  curlvxSum->v = vec_mask_add(curlvxSum->v,
                              vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)), mask);
  curlvySum->v = vec_mask_add(curlvySum->v,
                              vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)), mask);
  curlvzSum->v = vec_mask_add(curlvzSum->v,
                              vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)), mask);
}

/**
 * @brief Density interaction computed using 2 interleaved vectors
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_2_vec_density(float *R2, float *Dx, float *Dy, float *Dz,
                                 vector hi_inv, vector vix, vector viy,
                                 vector viz, float *Vjx, float *Vjy, float *Vjz,
                                 float *Mj, vector *rhoSum, vector *rho_dhSum,
                                 vector *wcountSum, vector *wcount_dhSum,
                                 vector *div_vSum, vector *curlvxSum,
                                 vector *curlvySum, vector *curlvzSum,
                                 mask_t mask, mask_t mask2, short mask_cond) {

  vector r, ri, r2, ui, wi, wi_dx;
  vector mj;
  vector dx, dy, dz, dvx, dvy, dvz;
  vector vjx, vjy, vjz;
  vector dvdr;
  vector curlvrx, curlvry, curlvrz;
  vector r_2, ri2, r2_2, ui2, wi2, wi_dx2;
  vector mj2;
  vector dx2, dy2, dz2, dvx2, dvy2, dvz2;
  vector vjx2, vjy2, vjz2;
  vector dvdr2;
  vector curlvrx2, curlvry2, curlvrz2;

  /* Fill the vectors. */
  mj.v = vec_load(Mj);
  mj2.v = vec_load(&Mj[VEC_SIZE]);
  vjx.v = vec_load(Vjx);
  vjx2.v = vec_load(&Vjx[VEC_SIZE]);
  vjy.v = vec_load(Vjy);
  vjy2.v = vec_load(&Vjy[VEC_SIZE]);
  vjz.v = vec_load(Vjz);
  vjz2.v = vec_load(&Vjz[VEC_SIZE]);
  dx.v = vec_load(Dx);
  dx2.v = vec_load(&Dx[VEC_SIZE]);
  dy.v = vec_load(Dy);
  dy2.v = vec_load(&Dy[VEC_SIZE]);
  dz.v = vec_load(Dz);
  dz2.v = vec_load(&Dz[VEC_SIZE]);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  r2_2.v = vec_load(&R2[VEC_SIZE]);
  ri = vec_reciprocal_sqrt(r2);
  ri2 = vec_reciprocal_sqrt(r2_2);
  r.v = vec_mul(r2.v, ri.v);
  r_2.v = vec_mul(r2_2.v, ri2.v);

  ui.v = vec_mul(r.v, hi_inv.v);
  ui2.v = vec_mul(r_2.v, hi_inv.v);

  /* Calculate the kernel for two particles. */
  kernel_deval_2_vec(&ui, &wi, &wi_dx, &ui2, &wi2, &wi_dx2);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvx2.v = vec_sub(vix.v, vjx2.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvy2.v = vec_sub(viy.v, vjy2.v);
  dvz.v = vec_sub(viz.v, vjz.v);
  dvz2.v = vec_sub(viz.v, vjz2.v);

  /* Compute dv dot r */
  dvdr.v = vec_fma(dvx.v, dx.v, vec_fma(dvy.v, dy.v, vec_mul(dvz.v, dz.v)));
  dvdr2.v =
      vec_fma(dvx2.v, dx2.v, vec_fma(dvy2.v, dy2.v, vec_mul(dvz2.v, dz2.v)));
  dvdr.v = vec_mul(dvdr.v, ri.v);
  dvdr2.v = vec_mul(dvdr2.v, ri2.v);

  /* Compute dv cross r */
  curlvrx.v =
      vec_fma(dvy.v, dz.v, vec_mul(vec_set1(-1.0f), vec_mul(dvz.v, dy.v)));
  curlvrx2.v =
      vec_fma(dvy2.v, dz2.v, vec_mul(vec_set1(-1.0f), vec_mul(dvz2.v, dy2.v)));
  curlvry.v =
      vec_fma(dvz.v, dx.v, vec_mul(vec_set1(-1.0f), vec_mul(dvx.v, dz.v)));
  curlvry2.v =
      vec_fma(dvz2.v, dx2.v, vec_mul(vec_set1(-1.0f), vec_mul(dvx2.v, dz2.v)));
  curlvrz.v =
      vec_fma(dvx.v, dy.v, vec_mul(vec_set1(-1.0f), vec_mul(dvy.v, dx.v)));
  curlvrz2.v =
      vec_fma(dvx2.v, dy2.v, vec_mul(vec_set1(-1.0f), vec_mul(dvy2.v, dx2.v)));
  curlvrx.v = vec_mul(curlvrx.v, ri.v);
  curlvrx2.v = vec_mul(curlvrx2.v, ri2.v);
  curlvry.v = vec_mul(curlvry.v, ri.v);
  curlvry2.v = vec_mul(curlvry2.v, ri2.v);
  curlvrz.v = vec_mul(curlvrz.v, ri.v);
  curlvrz2.v = vec_mul(curlvrz2.v, ri2.v);

  vector wcount_dh_update, wcount_dh_update2;
  wcount_dh_update.v =
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v));
  wcount_dh_update2.v =
      vec_fma(vec_set1(hydro_dimension), wi2.v, vec_mul(ui2.v, wi_dx2.v));

  /* Mask updates to intermediate vector sums for particle pi. */
  /* Mask only when needed. */
  if (mask_cond) {
    rhoSum->v = vec_mask_add(rhoSum->v, vec_mul(mj.v, wi.v), mask);
    rhoSum->v = vec_mask_add(rhoSum->v, vec_mul(mj2.v, wi2.v), mask2);
    rho_dhSum->v =
        vec_mask_sub(rho_dhSum->v, vec_mul(mj.v, wcount_dh_update.v), mask);
    rho_dhSum->v =
        vec_mask_sub(rho_dhSum->v, vec_mul(mj2.v, wcount_dh_update2.v), mask2);
    wcountSum->v = vec_mask_add(wcountSum->v, wi.v, mask);
    wcountSum->v = vec_mask_add(wcountSum->v, wi2.v, mask2);
    wcount_dhSum->v = vec_mask_sub(wcount_dhSum->v, wcount_dh_update.v, mask);
    wcount_dhSum->v = vec_mask_sub(wcount_dhSum->v, wcount_dh_update2.v, mask2);
    div_vSum->v = vec_mask_sub(div_vSum->v,
                               vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)), mask);
    div_vSum->v = vec_mask_sub(
        div_vSum->v, vec_mul(mj2.v, vec_mul(dvdr2.v, wi_dx2.v)), mask2);
    curlvxSum->v = vec_mask_add(
        curlvxSum->v, vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)), mask);
    curlvxSum->v = vec_mask_add(
        curlvxSum->v, vec_mul(mj2.v, vec_mul(curlvrx2.v, wi_dx2.v)), mask2);
    curlvySum->v = vec_mask_add(
        curlvySum->v, vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)), mask);
    curlvySum->v = vec_mask_add(
        curlvySum->v, vec_mul(mj2.v, vec_mul(curlvry2.v, wi_dx2.v)), mask2);
    curlvzSum->v = vec_mask_add(
        curlvzSum->v, vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)), mask);
    curlvzSum->v = vec_mask_add(
        curlvzSum->v, vec_mul(mj2.v, vec_mul(curlvrz2.v, wi_dx2.v)), mask2);
  } else {
    rhoSum->v = vec_add(rhoSum->v, vec_mul(mj.v, wi.v));
    rhoSum->v = vec_add(rhoSum->v, vec_mul(mj2.v, wi2.v));
    rho_dhSum->v = vec_sub(rho_dhSum->v, vec_mul(mj.v, wcount_dh_update.v));
    rho_dhSum->v = vec_sub(rho_dhSum->v, vec_mul(mj2.v, wcount_dh_update2.v));
    wcountSum->v = vec_add(wcountSum->v, wi.v);
    wcountSum->v = vec_add(wcountSum->v, wi2.v);
    wcount_dhSum->v = vec_sub(wcount_dhSum->v, wcount_dh_update.v);
    wcount_dhSum->v = vec_sub(wcount_dhSum->v, wcount_dh_update2.v);
    div_vSum->v = vec_sub(div_vSum->v, vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)));
    div_vSum->v =
        vec_sub(div_vSum->v, vec_mul(mj2.v, vec_mul(dvdr2.v, wi_dx2.v)));
    curlvxSum->v =
        vec_add(curlvxSum->v, vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)));
    curlvxSum->v =
        vec_add(curlvxSum->v, vec_mul(mj2.v, vec_mul(curlvrx2.v, wi_dx2.v)));
    curlvySum->v =
        vec_add(curlvySum->v, vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)));
    curlvySum->v =
        vec_add(curlvySum->v, vec_mul(mj2.v, vec_mul(curlvry2.v, wi_dx2.v)));
    curlvzSum->v =
        vec_add(curlvzSum->v, vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)));
    curlvzSum->v =
        vec_add(curlvzSum->v, vec_mul(mj2.v, vec_mul(curlvrz2.v, wi_dx2.v)));
  }
}
#endif

/**
 * @brief Force loop
 */
__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute h-gradient terms */
  const float f_i = pi->force.f;
  const float f_j = pj->force.f;

  /* Compute pressure terms */
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = (dvdr < 0.f) ? dvdr : 0.f;
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij *
                     (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      (f_i * P_over_rho2_i * wi_dr + f_j * P_over_rho2_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  pj->a_hydro[0] += mi * acc * dx[0];
  pj->a_hydro[1] += mi * acc * dx[1];
  pj->a_hydro[2] += mi * acc * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = (pi->force.v_sig > v_sig) ? pi->force.v_sig : v_sig;
  pj->force.v_sig = (pj->force.v_sig > v_sig) ? pj->force.v_sig : v_sig;

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr;
  pj->entropy_dt += mi * visc_term * dvdr;
}

/**
 * @brief Force loop (non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hjd_inv * wj_dx;

  /* Compute h-gradient terms */
  const float f_i = pi->force.f;
  const float f_j = pj->force.f;

  /* Compute pressure terms */
  const float P_over_rho2_i = pi->force.P_over_rho2;
  const float P_over_rho2_j = pj->force.P_over_rho2;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = (dvdr < 0.f) ? dvdr : 0.f;
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij *
                     (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      (f_i * P_over_rho2_i * wi_dr + f_j * P_over_rho2_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = (pi->force.v_sig > v_sig) ? pi->force.v_sig : v_sig;

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr;
}

#ifdef WITH_VECTORIZATION
static const vector const_viscosity_alpha_fac =
    FILL_VEC(-0.25f * const_viscosity_alpha);

/**
 * @brief Force interaction computed using 1 vector
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_1_vec_force(
    vector *r2, vector *dx, vector *dy, vector *dz, vector vix, vector viy,
    vector viz, vector pirho, vector grad_hi, vector piPOrho2, vector balsara_i,
    vector ci, float *Vjx, float *Vjy, float *Vjz, float *Pjrho, float *Grad_hj,
    float *PjPOrho2, float *Balsara_j, float *Cj, float *Mj, vector hi_inv,
    vector hj_inv, vector *a_hydro_xSum, vector *a_hydro_ySum,
    vector *a_hydro_zSum, vector *h_dtSum, vector *v_sigSum,
    vector *entropy_dtSum, mask_t mask) {

#ifdef WITH_VECTORIZATION

  vector r, ri;
  vector vjx, vjy, vjz, dvx, dvy, dvz;
  vector pjrho, grad_hj, pjPOrho2, balsara_j, cj, mj;
  vector xi, xj;
  vector hid_inv, hjd_inv;
  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector piax, piay, piaz;
  vector pih_dt;
  vector v_sig;
  vector omega_ij, mu_ij, fac_mu, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;

  /* Fill vectors. */
  vjx.v = vec_load(Vjx);
  vjy.v = vec_load(Vjy);
  vjz.v = vec_load(Vjz);
  mj.v = vec_load(Mj);

  pjrho.v = vec_load(Pjrho);
  grad_hj.v = vec_load(Grad_hj);
  pjPOrho2.v = vec_load(PjPOrho2);
  balsara_j.v = vec_load(Balsara_j);
  cj.v = vec_load(Cj);

  fac_mu.v = vec_set1(1.f); /* Will change with cosmological integration */

  /* Load stuff. */
  balsara.v = vec_add(balsara_i.v, balsara_j.v);

  /* Get the radius and inverse radius. */
  ri = vec_reciprocal_sqrt(*r2);
  r.v = vec_mul(r2->v, ri.v);

  /* Get the kernel for hi. */
  hid_inv = pow_dimension_plus_one_vec(hi_inv);
  xi.v = vec_mul(r.v, hi_inv.v);
  kernel_eval_dWdx_force_vec(&xi, &wi_dx);
  wi_dr.v = vec_mul(hid_inv.v, wi_dx.v);

  /* Get the kernel for hj. */
  hjd_inv = pow_dimension_plus_one_vec(hj_inv);
  xj.v = vec_mul(r.v, hj_inv.v);

  /* Calculate the kernel for two particles. */
  kernel_eval_dWdx_force_vec(&xj, &wj_dx);

  wj_dr.v = vec_mul(hjd_inv.v, wj_dx.v);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvz.v = vec_sub(viz.v, vjz.v);

  /* Compute dv dot r. */
  dvdr.v = vec_fma(dvx.v, dx->v, vec_fma(dvy.v, dy->v, vec_mul(dvz.v, dz->v)));

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_setzero());
  mu_ij.v =
      vec_mul(fac_mu.v, vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v = vec_fnma(vec_set1(3.f), mu_ij.v, vec_add(ci.v, cj.v));

  /* Now construct the full viscosity term */
  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho.v));
  visc.v = vec_div(vec_mul(const_viscosity_alpha_fac.v,
                           vec_mul(v_sig.v, vec_mul(mu_ij.v, balsara.v))),
                   rho_ij.v);

  /* Now, convolve with the kernel */
  visc_term.v =
      vec_mul(vec_set1(0.5f),
              vec_mul(visc.v, vec_mul(vec_add(wi_dr.v, wj_dr.v), ri.v)));

  sph_term.v =
      vec_mul(vec_fma(vec_mul(grad_hi.v, piPOrho2.v), wi_dr.v,
                      vec_mul(grad_hj.v, vec_mul(pjPOrho2.v, wj_dr.v))),
              ri.v);

  /* Eventually get the acceleration */
  acc.v = vec_add(visc_term.v, sph_term.v);

  /* Use the force, Luke! */
  piax.v = vec_mul(mj.v, vec_mul(dx->v, acc.v));
  piay.v = vec_mul(mj.v, vec_mul(dy->v, acc.v));
  piaz.v = vec_mul(mj.v, vec_mul(dz->v, acc.v));

  /* Get the time derivative for h. */
  pih_dt.v =
      vec_div(vec_mul(mj.v, vec_mul(dvdr.v, vec_mul(ri.v, wi_dr.v))), pjrho.v);

  /* Change in entropy */
  entropy_dt.v = vec_mul(mj.v, vec_mul(visc_term.v, dvdr.v));

  /* Store the forces back on the particles. */
  a_hydro_xSum->v = vec_mask_sub(a_hydro_xSum->v, piax.v, mask);
  a_hydro_ySum->v = vec_mask_sub(a_hydro_ySum->v, piay.v, mask);
  a_hydro_zSum->v = vec_mask_sub(a_hydro_zSum->v, piaz.v, mask);
  h_dtSum->v = vec_mask_sub(h_dtSum->v, pih_dt.v, mask);
  v_sigSum->v = vec_fmax(v_sigSum->v, vec_and_mask(v_sig.v, mask));
  entropy_dtSum->v = vec_mask_add(entropy_dtSum->v, entropy_dt.v, mask);

#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
      "the vectorised version should have been used.");

#endif
}

/**
 * @brief Force interaction computed using 2 interleaved vectors
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_2_vec_force(
    float *R2, float *Dx, float *Dy, float *Dz, vector vix, vector viy,
    vector viz, vector pirho, vector grad_hi, vector piPOrho2, vector balsara_i,
    vector ci, float *Vjx, float *Vjy, float *Vjz, float *Pjrho, float *Grad_hj,
    float *PjPOrho2, float *Balsara_j, float *Cj, float *Mj, vector hi_inv,
    float *Hj_inv, vector *a_hydro_xSum, vector *a_hydro_ySum,
    vector *a_hydro_zSum, vector *h_dtSum, vector *v_sigSum,
    vector *entropy_dtSum, mask_t mask, mask_t mask_2, short mask_cond) {

#ifdef WITH_VECTORIZATION

  vector r, r2, ri;
  vector dx, dy, dz, dvx, dvy, dvz;
  vector vjx, vjy, vjz;
  vector pjrho, grad_hj, pjPOrho2, balsara_j, cj, mj, hj_inv;
  vector ui, uj;
  vector hid_inv, hjd_inv;
  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector piax, piay, piaz;
  vector pih_dt;
  vector v_sig;
  vector omega_ij, mu_ij, fac_mu, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;

  vector r_2, r2_2, ri_2;
  vector dx_2, dy_2, dz_2, dvx_2, dvy_2, dvz_2;
  vector vjx_2, vjy_2, vjz_2;
  vector pjrho_2, grad_hj_2, pjPOrho2_2, balsara_j_2, cj_2, mj_2, hj_inv_2;
  vector ui_2, uj_2;
  vector hjd_inv_2;
  vector wi_dx_2, wj_dx_2, wi_dr_2, wj_dr_2, dvdr_2;
  vector piax_2, piay_2, piaz_2;
  vector pih_dt_2;
  vector v_sig_2;
  vector omega_ij_2, mu_ij_2, balsara_2;
  vector rho_ij_2, visc_2, visc_term_2, sph_term_2, acc_2, entropy_dt_2;

  /* Fill vectors. */
  mj.v = vec_load(Mj);
  mj_2.v = vec_load(&Mj[VEC_SIZE]);
  vjx.v = vec_load(Vjx);
  vjx_2.v = vec_load(&Vjx[VEC_SIZE]);
  vjy.v = vec_load(Vjy);
  vjy_2.v = vec_load(&Vjy[VEC_SIZE]);
  vjz.v = vec_load(Vjz);
  vjz_2.v = vec_load(&Vjz[VEC_SIZE]);
  dx.v = vec_load(Dx);
  dx_2.v = vec_load(&Dx[VEC_SIZE]);
  dy.v = vec_load(Dy);
  dy_2.v = vec_load(&Dy[VEC_SIZE]);
  dz.v = vec_load(Dz);
  dz_2.v = vec_load(&Dz[VEC_SIZE]);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  r2_2.v = vec_load(&R2[VEC_SIZE]);
  ri = vec_reciprocal_sqrt(r2);
  ri_2 = vec_reciprocal_sqrt(r2_2);
  r.v = vec_mul(r2.v, ri.v);
  r_2.v = vec_mul(r2_2.v, ri_2.v);

  /* Get remaining properties. */
  pjrho.v = vec_load(Pjrho);
  pjrho_2.v = vec_load(&Pjrho[VEC_SIZE]);
  grad_hj.v = vec_load(Grad_hj);
  grad_hj_2.v = vec_load(&Grad_hj[VEC_SIZE]);
  pjPOrho2.v = vec_load(PjPOrho2);
  pjPOrho2_2.v = vec_load(&PjPOrho2[VEC_SIZE]);
  balsara_j.v = vec_load(Balsara_j);
  balsara_j_2.v = vec_load(&Balsara_j[VEC_SIZE]);
  cj.v = vec_load(Cj);
  cj_2.v = vec_load(&Cj[VEC_SIZE]);
  hj_inv.v = vec_load(Hj_inv);
  hj_inv_2.v = vec_load(&Hj_inv[VEC_SIZE]);

  fac_mu.v = vec_set1(1.f); /* Will change with cosmological integration */

  /* Find the balsara switch. */
  balsara.v = vec_add(balsara_i.v, balsara_j.v);
  balsara_2.v = vec_add(balsara_i.v, balsara_j_2.v);

  /* Get the kernel for hi. */
  hid_inv = pow_dimension_plus_one_vec(hi_inv);
  ui.v = vec_mul(r.v, hi_inv.v);
  ui_2.v = vec_mul(r_2.v, hi_inv.v);
  kernel_eval_dWdx_force_vec(&ui, &wi_dx);
  kernel_eval_dWdx_force_vec(&ui_2, &wi_dx_2);
  wi_dr.v = vec_mul(hid_inv.v, wi_dx.v);
  wi_dr_2.v = vec_mul(hid_inv.v, wi_dx_2.v);

  /* Get the kernel for hj. */
  hjd_inv = pow_dimension_plus_one_vec(hj_inv);
  hjd_inv_2 = pow_dimension_plus_one_vec(hj_inv_2);
  uj.v = vec_mul(r.v, hj_inv.v);
  uj_2.v = vec_mul(r_2.v, hj_inv_2.v);

  /* Calculate the kernel for two particles. */
  kernel_eval_dWdx_force_vec(&uj, &wj_dx);
  kernel_eval_dWdx_force_vec(&uj_2, &wj_dx_2);

  wj_dr.v = vec_mul(hjd_inv.v, wj_dx.v);
  wj_dr_2.v = vec_mul(hjd_inv_2.v, wj_dx_2.v);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvx_2.v = vec_sub(vix.v, vjx_2.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvy_2.v = vec_sub(viy.v, vjy_2.v);
  dvz.v = vec_sub(viz.v, vjz.v);
  dvz_2.v = vec_sub(viz.v, vjz_2.v);

  /* Compute dv dot r. */
  dvdr.v = vec_fma(dvx.v, dx.v, vec_fma(dvy.v, dy.v, vec_mul(dvz.v, dz.v)));
  dvdr_2.v = vec_fma(dvx_2.v, dx_2.v,
                     vec_fma(dvy_2.v, dy_2.v, vec_mul(dvz_2.v, dz_2.v)));

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_setzero());
  omega_ij_2.v = vec_fmin(dvdr_2.v, vec_setzero());
  mu_ij.v =
      vec_mul(fac_mu.v, vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */
  mu_ij_2.v = vec_mul(
      fac_mu.v, vec_mul(ri_2.v, omega_ij_2.v)); /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v = vec_fnma(vec_set1(3.f), mu_ij.v, vec_add(ci.v, cj.v));
  v_sig_2.v = vec_fnma(vec_set1(3.f), mu_ij_2.v, vec_add(ci.v, cj_2.v));

  /* Now construct the full viscosity term */
  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho.v));
  rho_ij_2.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho_2.v));

  visc.v = vec_div(vec_mul(const_viscosity_alpha_fac.v,
                           vec_mul(v_sig.v, vec_mul(mu_ij.v, balsara.v))),
                   rho_ij.v);
  visc_2.v =
      vec_div(vec_mul(const_viscosity_alpha_fac.v,
                      vec_mul(v_sig_2.v, vec_mul(mu_ij_2.v, balsara_2.v))),
              rho_ij_2.v);

  /* Now, convolve with the kernel */
  visc_term.v =
      vec_mul(vec_set1(0.5f),
              vec_mul(visc.v, vec_mul(vec_add(wi_dr.v, wj_dr.v), ri.v)));
  visc_term_2.v = vec_mul(
      vec_set1(0.5f),
      vec_mul(visc_2.v, vec_mul(vec_add(wi_dr_2.v, wj_dr_2.v), ri_2.v)));

  vector grad_hi_mul_piPOrho2;
  grad_hi_mul_piPOrho2.v = vec_mul(grad_hi.v, piPOrho2.v);

  sph_term.v =
      vec_mul(vec_fma(grad_hi_mul_piPOrho2.v, wi_dr.v,
                      vec_mul(grad_hj.v, vec_mul(pjPOrho2.v, wj_dr.v))),
              ri.v);
  sph_term_2.v =
      vec_mul(vec_fma(grad_hi_mul_piPOrho2.v, wi_dr_2.v,
                      vec_mul(grad_hj_2.v, vec_mul(pjPOrho2_2.v, wj_dr_2.v))),
              ri_2.v);

  /* Eventually get the acceleration */
  acc.v = vec_add(visc_term.v, sph_term.v);
  acc_2.v = vec_add(visc_term_2.v, sph_term_2.v);

  /* Use the force, Luke! */
  piax.v = vec_mul(mj.v, vec_mul(dx.v, acc.v));
  piax_2.v = vec_mul(mj_2.v, vec_mul(dx_2.v, acc_2.v));
  piay.v = vec_mul(mj.v, vec_mul(dy.v, acc.v));
  piay_2.v = vec_mul(mj_2.v, vec_mul(dy_2.v, acc_2.v));
  piaz.v = vec_mul(mj.v, vec_mul(dz.v, acc.v));
  piaz_2.v = vec_mul(mj_2.v, vec_mul(dz_2.v, acc_2.v));

  /* Get the time derivative for h. */
  pih_dt.v =
      vec_div(vec_mul(mj.v, vec_mul(dvdr.v, vec_mul(ri.v, wi_dr.v))), pjrho.v);
  pih_dt_2.v =
      vec_div(vec_mul(mj_2.v, vec_mul(dvdr_2.v, vec_mul(ri_2.v, wi_dr_2.v))),
              pjrho_2.v);

  /* Change in entropy */
  entropy_dt.v = vec_mul(mj.v, vec_mul(visc_term.v, dvdr.v));
  entropy_dt_2.v = vec_mul(mj_2.v, vec_mul(visc_term_2.v, dvdr_2.v));

  /* Store the forces back on the particles. */
  if (mask_cond) {
    a_hydro_xSum->v = vec_mask_sub(a_hydro_xSum->v, piax.v, mask);
    a_hydro_xSum->v = vec_mask_sub(a_hydro_xSum->v, piax_2.v, mask_2);
    a_hydro_ySum->v = vec_mask_sub(a_hydro_ySum->v, piay.v, mask);
    a_hydro_ySum->v = vec_mask_sub(a_hydro_ySum->v, piay_2.v, mask_2);
    a_hydro_zSum->v = vec_mask_sub(a_hydro_zSum->v, piaz.v, mask);
    a_hydro_zSum->v = vec_mask_sub(a_hydro_zSum->v, piaz_2.v, mask_2);
    h_dtSum->v = vec_mask_sub(h_dtSum->v, pih_dt.v, mask);
    h_dtSum->v = vec_mask_sub(h_dtSum->v, pih_dt_2.v, mask_2);
    v_sigSum->v = vec_fmax(v_sigSum->v, vec_and_mask(v_sig.v, mask));
    v_sigSum->v = vec_fmax(v_sigSum->v, vec_and_mask(v_sig_2.v, mask_2));
    entropy_dtSum->v = vec_mask_add(entropy_dtSum->v, entropy_dt.v, mask);
    entropy_dtSum->v = vec_mask_add(entropy_dtSum->v, entropy_dt_2.v, mask_2);
  } else {
    a_hydro_xSum->v = vec_sub(a_hydro_xSum->v, piax.v);
    a_hydro_xSum->v = vec_sub(a_hydro_xSum->v, piax_2.v);
    a_hydro_ySum->v = vec_sub(a_hydro_ySum->v, piay.v);
    a_hydro_ySum->v = vec_sub(a_hydro_ySum->v, piay_2.v);
    a_hydro_zSum->v = vec_sub(a_hydro_zSum->v, piaz.v);
    a_hydro_zSum->v = vec_sub(a_hydro_zSum->v, piaz_2.v);
    h_dtSum->v = vec_sub(h_dtSum->v, pih_dt.v);
    h_dtSum->v = vec_sub(h_dtSum->v, pih_dt_2.v);
    v_sigSum->v = vec_fmax(v_sigSum->v, v_sig.v);
    v_sigSum->v = vec_fmax(v_sigSum->v, v_sig_2.v);
    entropy_dtSum->v = vec_add(entropy_dtSum->v, entropy_dt.v);
    entropy_dtSum->v = vec_add(entropy_dtSum->v, entropy_dt_2.v);
  }
#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
      "the vectorised version should have been used.");

#endif
}

#endif

#endif /* SWIFT_GADGET2_HYDRO_IACT_H */
