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
 * @brief Density loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_density(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef WITH_OLD_VECTORIZATION

  vector r, ri, r2, ui, uj, hi, hj, hi_inv, hj_inv, wi, wj, wi_dx, wj_dx;
  vector rhoi, rhoj, rhoi_dh, rhoj_dh, wcounti, wcountj, wcounti_dh, wcountj_dh;
  vector mi, mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr, div_vi, div_vj;
  vector curlvr[3], curl_vi[3], curl_vj[3];
  int k, j;

#if VEC_SIZE == 8
  /* Get the masses. */
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  /* Get each velocity component. */
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  /* Get each component of particle separation.
   * (Dx={dx1,dy1,dz1,dx2,dy2,dz2,...,dxn,dyn,dzn})*/
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
#else
  error("Unknown vector size.");
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri = vec_reciprocal_sqrt(r2);
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv = vec_reciprocal(hi);
  ui.v = r.v * hi_inv.v;

  hj.v = vec_load(Hj);
  hj_inv = vec_reciprocal(hj);
  uj.v = r.v * hj_inv.v;

  /* Compute the kernel function. */
  kernel_deval_vec(&ui, &wi, &wi_dx);
  kernel_deval_vec(&uj, &wj, &wj_dx);

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;
  dv[1].v = vi[1].v - vj[1].v;
  dv[2].v = vi[2].v - vj[2].v;

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;

  /* Compute density of pi. */
  rhoi.v = mj.v * wi.v;
  rhoi_dh.v = mj.v * (vec_set1(hydro_dimension) * wi.v + ui.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = (vec_set1(hydro_dimension) * wi.v + ui.v * wi_dx.v);
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  /* Compute density of pj. */
  rhoj.v = mi.v * wj.v;
  rhoj_dh.v = mi.v * (vec_set1(hydro_dimension) * wj.v + uj.v * wj_dx.v);
  wcountj.v = wj.v;
  wcountj_dh.v = (vec_set1(hydro_dimension) * wj.v + uj.v * wj_dx.v);
  div_vj.v = mi.v * dvdr.v * wj_dx.v;
  for (k = 0; k < 3; k++) curl_vj[k].v = mi.v * curlvr[k].v * wj_dx.v;

  /* Update particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->density.rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->density.div_v -= div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.rot_v[j] += curl_vi[j].f[k];
    pj[k]->rho += rhoj.f[k];
    pj[k]->density.rho_dh -= rhoj_dh.f[k];
    pj[k]->density.wcount += wcountj.f[k];
    pj[k]->density.wcount_dh -= wcountj_dh.f[k];
    pj[k]->density.div_v -= div_vj.f[k];
    for (j = 0; j < 3; j++) pj[k]->density.rot_v[j] += curl_vj[j].f[k];
  }

#else

  error(
      "The Gadget2 serial version of runner_iact_density was called when the "
      "vectorised version should have been used.");

#endif
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

/**
 * @brief Density loop (non-symmetric vectorized version)
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_vec_density(float *R2, float *Dx, float *Hi, float *Hj,
                               struct part **pi, struct part **pj) {

#ifdef WITH_OLD_VECTORIZATION

  vector r, ri, r2, ui, hi, hi_inv, wi, wi_dx;
  vector rhoi, rhoi_dh, wcounti, wcounti_dh, div_vi;
  vector mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr;
  vector curlvr[3], curl_vi[3];
  int k, j;

#if VEC_SIZE == 8
  /* Get the masses. */
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  /* Get each velocity component. */
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  /* Get each component of particle separation.
   * (Dx={dx1,dy1,dz1,dx2,dy2,dz2,...,dxn,dyn,dzn})*/
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
#else
  error("Unknown vector size.");
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri = vec_reciprocal_sqrt(r2);
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv = vec_reciprocal(hi);
  ui.v = r.v * hi_inv.v;

  kernel_deval_vec(&ui, &wi, &wi_dx);

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;
  dv[1].v = vi[1].v - vj[1].v;
  dv[2].v = vi[2].v - vj[2].v;

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);
  dvdr.v = dvdr.v * ri.v;

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;

  /* Compute density of pi. */
  rhoi.v = mj.v * wi.v;
  rhoi_dh.v = mj.v * (vec_set1(hydro_dimension) * wi.v + ui.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = (vec_set1(hydro_dimension) * wi.v + ui.v * wi_dx.v);
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  /* Update particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->density.rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->density.div_v -= div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.rot_v[j] += curl_vi[j].f[k];
  }

#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_density was called "
      "when the vectorised version should have been used.");

#endif
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
                                 vector mask, int knlMask) {

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

/* Mask updates to intermediate vector sums for particle pi. */
#ifdef HAVE_AVX512_F
  rhoSum->v =
      _mm512_mask_add_ps(rhoSum->v, knlMask, vec_mul(mj.v, wi.v), rhoSum->v);

  rho_dhSum->v =
      _mm512_mask_sub_ps(rho_dhSum->v, knlMask, rho_dhSum->v,
                         vec_mul(mj.v, vec_fma(vec_set1(hydro_dimension), wi.v,
                                               vec_mul(ui.v, wi_dx.v))));

  wcountSum->v = _mm512_mask_add_ps(wcountSum->v, knlMask, wi.v, wcountSum->v);

  wcount_dhSum->v = _mm512_mask_sub_ps(
      rho_dhSum->v, knlMask, rho_dhSum->v,
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v)));

  div_vSum->v = _mm512_mask_sub_ps(div_vSum->v, knlMask, div_vSum->v,
                                   vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)));

  curlvxSum->v = _mm512_mask_add_ps(curlvxSum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)),
                                    curlvxSum->v);

  curlvySum->v = _mm512_mask_add_ps(curlvySum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)),
                                    curlvySum->v);

  curlvzSum->v = _mm512_mask_add_ps(curlvzSum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)),
                                    curlvzSum->v);
#else
  rhoSum->v += vec_and(vec_mul(mj.v, wi.v), mask.v);
  rho_dhSum->v -= vec_and(vec_mul(mj.v, vec_fma(vec_set1(hydro_dimension), wi.v,
                                                vec_mul(ui.v, wi_dx.v))),
                          mask.v);
  wcountSum->v += vec_and(wi.v, mask.v);
  wcount_dhSum->v -= vec_and(
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v)), mask.v);
  div_vSum->v -= vec_and(vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)), mask.v);
  curlvxSum->v += vec_and(vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)), mask.v);
  curlvySum->v += vec_and(vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)), mask.v);
  curlvzSum->v += vec_and(vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)), mask.v);
#endif
}

/**
 * @brief Density interaction computed using 2 interleaved vectors
 * (non-symmetric vectorized version).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_2_vec_density(
    float *R2, float *Dx, float *Dy, float *Dz, vector hi_inv, vector vix,
    vector viy, vector viz, float *Vjx, float *Vjy, float *Vjz, float *Mj,
    vector *rhoSum, vector *rho_dhSum, vector *wcountSum, vector *wcount_dhSum,
    vector *div_vSum, vector *curlvxSum, vector *curlvySum, vector *curlvzSum,
    vector mask, vector mask2, int knlMask, int knlMask2) {

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

/* Mask updates to intermediate vector sums for particle pi. */
#ifdef HAVE_AVX512_F
  rhoSum->v =
      _mm512_mask_add_ps(rhoSum->v, knlMask, vec_mul(mj.v, wi.v), rhoSum->v);
  rhoSum->v =
      _mm512_mask_add_ps(rhoSum->v, knlMask2, vec_mul(mj2.v, wi2.v), rhoSum->v);

  rho_dhSum->v =
      _mm512_mask_sub_ps(rho_dhSum->v, knlMask, rho_dhSum->v,
                         vec_mul(mj.v, vec_fma(vec_set1(hydro_dimension), wi.v,
                                               vec_mul(ui.v, wi_dx.v))));
  rho_dhSum->v =
      _mm512_mask_sub_ps(rho_dhSum->v, knlMask2, rho_dhSum->v,
                         vec_mul(mj.v, vec_fma(vec_set1(hydro_dimension), wi2.v,
                                               vec_mul(ui2.v, wi_dx2.v))));

  wcountSum->v = _mm512_mask_add_ps(wcountSum->v, knlMask, wi.v, wcountSum->v);
  wcountSum->v =
      _mm512_mask_add_ps(wcountSum->v, knlMask2, wi2.v, wcountSum->v);

  wcount_dhSum->v = _mm512_mask_sub_ps(
      rho_dhSum->v, knlMask, rho_dhSum->v,
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v)));

  wcount_dhSum->v = _mm512_mask_sub_ps(
      rho_dhSum->v, knlMask2, rho_dhSum->v,
      vec_fma(vec_set1(hydro_dimension), wi2.v, vec_mul(ui2.v, wi_dx2.v)));

  div_vSum->v = _mm512_mask_sub_ps(div_vSum->v, knlMask, div_vSum->v,
                                   vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)));
  div_vSum->v = _mm512_mask_sub_ps(div_vSum->v, knlMask2, div_vSum->v,
                                   vec_mul(mj2.v, vec_mul(dvdr2.v, wi_dx2.v)));

  curlvxSum->v = _mm512_mask_add_ps(curlvxSum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)),
                                    curlvxSum->v);
  curlvxSum->v = _mm512_mask_add_ps(
      curlvxSum->v, knlMask2, vec_mul(mj2.v, vec_mul(curlvrx2.v, wi_dx2.v)),
      curlvxSum->v);

  curlvySum->v = _mm512_mask_add_ps(curlvySum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)),
                                    curlvySum->v);
  curlvySum->v = _mm512_mask_add_ps(
      curlvySum->v, knlMask2, vec_mul(mj2.v, vec_mul(curlvry2.v, wi_dx2.v)),
      curlvySum->v);

  curlvzSum->v = _mm512_mask_add_ps(curlvzSum->v, knlMask,
                                    vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)),
                                    curlvzSum->v);
  curlvzSum->v = _mm512_mask_add_ps(
      curlvzSum->v, knlMask2, vec_mul(mj2.v, vec_mul(curlvrz2.v, wi_dx2.v)),
      curlvzSum->v);
#else
  rhoSum->v += vec_and(vec_mul(mj.v, wi.v), mask.v);
  rhoSum->v += vec_and(vec_mul(mj2.v, wi2.v), mask2.v);
  rho_dhSum->v -= vec_and(vec_mul(mj.v, vec_fma(vec_set1(hydro_dimension), wi.v,
                                                vec_mul(ui.v, wi_dx.v))),
                          mask.v);
  rho_dhSum->v -=
      vec_and(vec_mul(mj2.v, vec_fma(vec_set1(hydro_dimension), wi2.v,
                                     vec_mul(ui2.v, wi_dx2.v))),
              mask2.v);
  wcountSum->v += vec_and(wi.v, mask.v);
  wcountSum->v += vec_and(wi2.v, mask2.v);
  wcount_dhSum->v -= vec_and(
      vec_fma(vec_set1(hydro_dimension), wi.v, vec_mul(ui.v, wi_dx.v)), mask.v);
  wcount_dhSum->v -= vec_and(
      vec_fma(vec_set1(hydro_dimension), wi2.v, vec_mul(ui2.v, wi_dx2.v)),
      mask2.v);
  div_vSum->v -= vec_and(vec_mul(mj.v, vec_mul(dvdr.v, wi_dx.v)), mask.v);
  div_vSum->v -= vec_and(vec_mul(mj2.v, vec_mul(dvdr2.v, wi_dx2.v)), mask2.v);
  curlvxSum->v += vec_and(vec_mul(mj.v, vec_mul(curlvrx.v, wi_dx.v)), mask.v);
  curlvxSum->v +=
      vec_and(vec_mul(mj2.v, vec_mul(curlvrx2.v, wi_dx2.v)), mask2.v);
  curlvySum->v += vec_and(vec_mul(mj.v, vec_mul(curlvry.v, wi_dx.v)), mask.v);
  curlvySum->v +=
      vec_and(vec_mul(mj2.v, vec_mul(curlvry2.v, wi_dx2.v)), mask2.v);
  curlvzSum->v += vec_and(vec_mul(mj.v, vec_mul(curlvrz.v, wi_dx.v)), mask.v);
  curlvzSum->v +=
      vec_and(vec_mul(mj2.v, vec_mul(curlvrz2.v, wi_dx2.v)), mask2.v);
#endif
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
 * @brief Force loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef WITH_VECTORIZATION

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hid_inv, hjd_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector piPOrho2, pjPOrho2, pirho, pjrho;
  vector mi, mj;
  vector f;
  vector grad_hi, grad_hj;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3], pja[3];
  vector pih_dt, pjh_dt;
  vector ci, cj, v_sig;
  vector omega_ij, mu_ij, fac_mu, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;
  int j, k;

  fac_mu.v = vec_set1(1.f); /* Will change with cosmological integration */

/* Load stuff. */
#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v = vec_set(pi[0]->force.P_over_rho2, pi[1]->force.P_over_rho2,
                       pi[2]->force.P_over_rho2, pi[3]->force.P_over_rho2,
                       pi[4]->force.P_over_rho2, pi[5]->force.P_over_rho2,
                       pi[6]->force.P_over_rho2, pi[7]->force.P_over_rho2);
  pjPOrho2.v = vec_set(pj[0]->force.P_over_rho2, pj[1]->force.P_over_rho2,
                       pj[2]->force.P_over_rho2, pj[3]->force.P_over_rho2,
                       pj[4]->force.P_over_rho2, pj[5]->force.P_over_rho2,
                       pj[6]->force.P_over_rho2, pj[7]->force.P_over_rho2);
  grad_hi.v =
      vec_set(pi[0]->force.f, pi[1]->force.f, pi[2]->force.f, pi[3]->force.f,
              pi[4]->force.f, pi[5]->force.f, pi[6]->force.f, pi[7]->force.f);
  grad_hj.v =
      vec_set(pj[0]->force.f, pj[1]->force.f, pj[2]->force.f, pj[3]->force.f,
              pj[4]->force.f, pj[5]->force.f, pj[6]->force.f, pj[7]->force.f);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  ci.v = vec_set(pi[0]->force.soundspeed, pi[1]->force.soundspeed,
                 pi[2]->force.soundspeed, pi[3]->force.soundspeed,
                 pi[4]->force.soundspeed, pi[5]->force.soundspeed,
                 pi[6]->force.soundspeed, pi[7]->force.soundspeed);
  cj.v = vec_set(pj[0]->force.soundspeed, pj[1]->force.soundspeed,
                 pj[2]->force.soundspeed, pj[3]->force.soundspeed,
                 pj[4]->force.soundspeed, pj[5]->force.soundspeed,
                 pj[6]->force.soundspeed, pj[7]->force.soundspeed);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.P_over_rho2, pi[1]->force.P_over_rho2,
                       pi[2]->force.P_over_rho2, pi[3]->force.P_over_rho2);
  pjPOrho2.v = vec_set(pj[0]->force.P_over_rho2, pj[1]->force.P_over_rho2,
                       pj[2]->force.P_over_rho2, pj[3]->force.P_over_rho2);
  grad_hi.v =
      vec_set(pi[0]->force.f, pi[1]->force.f, pi[2]->force.f, pi[3]->force.f);
  grad_hj.v =
      vec_set(pj[0]->force.f, pj[1]->force.f, pj[2]->force.f, pj[3]->force.f);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  ci.v = vec_set(pi[0]->force.soundspeed, pi[1]->force.soundspeed,
                 pi[2]->force.soundspeed, pi[3]->force.soundspeed);
  cj.v = vec_set(pj[0]->force.soundspeed, pj[1]->force.soundspeed,
                 pj[2]->force.soundspeed, pj[3]->force.soundspeed);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
#else
  error("Unknown vector size.");
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri = vec_reciprocal_sqrt(r2);
  r.v = r2.v * ri.v;

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv = vec_reciprocal(hi);
  hid_inv = pow_dimension_plus_one_vec(hi_inv); /* 1/h^(d+1) */
  xi.v = r.v * hi_inv.v;
  kernel_deval_vec(&xi, &wi, &wi_dx);
  wi_dr.v = hid_inv.v * wi_dx.v;

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv = vec_reciprocal(hj);
  hjd_inv = pow_dimension_plus_one_vec(hj_inv); /* 1/h^(d+1) */
  xj.v = r.v * hj_inv.v;
  kernel_deval_vec(&xj, &wj, &wj_dx);
  wj_dr.v = hjd_inv.v * wj_dx.v;

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +
           ((vi[2].v - vj[2].v) * dx[2].v);
  // dvdr.v = dvdr.v * ri.v;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));
  mu_ij.v = fac_mu.v * ri.v * omega_ij.v; /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(3.0f) * mu_ij.v;

  /* Now construct the full viscosity term */
  rho_ij.v = vec_set1(0.5f) * (pirho.v + pjrho.v);
  visc.v = vec_set1(-0.25f) * vec_set1(const_viscosity_alpha) * v_sig.v *
           mu_ij.v * balsara.v / rho_ij.v;

  /* Now, convolve with the kernel */
  visc_term.v = vec_set1(0.5f) * visc.v * (wi_dr.v + wj_dr.v) * ri.v;
  sph_term.v =
      (grad_hi.v * piPOrho2.v * wi_dr.v + grad_hj.v * pjPOrho2.v * wj_dr.v) *
      ri.v;

  /* Eventually get the acceleration */
  acc.v = visc_term.v + sph_term.v;

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f.v = dx[k].v * acc.v;
    pia[k].v = mj.v * f.v;
    pja[k].v = mi.v * f.v;
  }

  /* Get the time derivative for h. */
  pih_dt.v = mj.v * dvdr.v * ri.v / pjrho.v * wi_dr.v;
  pjh_dt.v = mi.v * dvdr.v * ri.v / pirho.v * wj_dr.v;

  /* Change in entropy */
  entropy_dt.v = visc_term.v * dvdr.v;

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    for (j = 0; j < 3; j++) {
      pi[k]->a_hydro[j] -= pia[j].f[k];
      pj[k]->a_hydro[j] += pja[j].f[k];
    }
    pi[k]->force.h_dt -= pih_dt.f[k];
    pj[k]->force.h_dt -= pjh_dt.f[k];
    pi[k]->force.v_sig = max(pi[k]->force.v_sig, v_sig.f[k]);
    pj[k]->force.v_sig = max(pj[k]->force.v_sig, v_sig.f[k]);
    pi[k]->entropy_dt += entropy_dt.f[k] * mj.f[k];
    pj[k]->entropy_dt += entropy_dt.f[k] * mi.f[k];
  }

#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
      "the vectorised version should have been used.");

#endif
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

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef WITH_VECTORIZATION

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hid_inv, hjd_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector piPOrho2, pjPOrho2, pirho, pjrho;
  vector mj;
  vector f;
  vector grad_hi, grad_hj;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3];
  vector pih_dt;
  vector ci, cj, v_sig;
  vector omega_ij, mu_ij, fac_mu, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;
  int j, k;

  fac_mu.v = vec_set1(1.f); /* Will change with cosmological integration */

/* Load stuff. */
#if VEC_SIZE == 8
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v = vec_set(pi[0]->force.P_over_rho2, pi[1]->force.P_over_rho2,
                       pi[2]->force.P_over_rho2, pi[3]->force.P_over_rho2,
                       pi[4]->force.P_over_rho2, pi[5]->force.P_over_rho2,
                       pi[6]->force.P_over_rho2, pi[7]->force.P_over_rho2);
  pjPOrho2.v = vec_set(pj[0]->force.P_over_rho2, pj[1]->force.P_over_rho2,
                       pj[2]->force.P_over_rho2, pj[3]->force.P_over_rho2,
                       pj[4]->force.P_over_rho2, pj[5]->force.P_over_rho2,
                       pj[6]->force.P_over_rho2, pj[7]->force.P_over_rho2);
  grad_hi.v =
      vec_set(pi[0]->force.f, pi[1]->force.f, pi[2]->force.f, pi[3]->force.f,
              pi[4]->force.f, pi[5]->force.f, pi[6]->force.f, pi[7]->force.f);
  grad_hj.v =
      vec_set(pj[0]->force.f, pj[1]->force.f, pj[2]->force.f, pj[3]->force.f,
              pj[4]->force.f, pj[5]->force.f, pj[6]->force.f, pj[7]->force.f);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  ci.v = vec_set(pi[0]->force.soundspeed, pi[1]->force.soundspeed,
                 pi[2]->force.soundspeed, pi[3]->force.soundspeed,
                 pi[4]->force.soundspeed, pi[5]->force.soundspeed,
                 pi[6]->force.soundspeed, pi[7]->force.soundspeed);
  cj.v = vec_set(pj[0]->force.soundspeed, pj[1]->force.soundspeed,
                 pj[2]->force.soundspeed, pj[3]->force.soundspeed,
                 pj[4]->force.soundspeed, pj[5]->force.soundspeed,
                 pj[6]->force.soundspeed, pj[7]->force.soundspeed);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k], Dx[12 + k],
                      Dx[15 + k], Dx[18 + k], Dx[21 + k]);
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.P_over_rho2, pi[1]->force.P_over_rho2,
                       pi[2]->force.P_over_rho2, pi[3]->force.P_over_rho2);
  pjPOrho2.v = vec_set(pj[0]->force.P_over_rho2, pj[1]->force.P_over_rho2,
                       pj[2]->force.P_over_rho2, pj[3]->force.P_over_rho2);
  grad_hi.v =
      vec_set(pi[0]->force.f, pi[1]->force.f, pi[2]->force.f, pi[3]->force.f);
  grad_hj.v =
      vec_set(pj[0]->force.f, pj[1]->force.f, pj[2]->force.f, pj[3]->force.f);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  ci.v = vec_set(pi[0]->force.soundspeed, pi[1]->force.soundspeed,
                 pi[2]->force.soundspeed, pi[3]->force.soundspeed);
  cj.v = vec_set(pj[0]->force.soundspeed, pj[1]->force.soundspeed,
                 pj[2]->force.soundspeed, pj[3]->force.soundspeed);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  for (k = 0; k < 3; k++)
    dx[k].v = vec_set(Dx[0 + k], Dx[3 + k], Dx[6 + k], Dx[9 + k]);
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
#else
  error("Unknown vector size.");
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri = vec_reciprocal_sqrt(r2);
  r.v = r2.v * ri.v;

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv = vec_reciprocal(hi);
  hid_inv = pow_dimension_plus_one_vec(hi_inv);
  xi.v = r.v * hi_inv.v;
  kernel_deval_vec(&xi, &wi, &wi_dx);
  wi_dr.v = hid_inv.v * wi_dx.v;

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv = vec_reciprocal(hj);
  hjd_inv = pow_dimension_plus_one_vec(hj_inv);
  xj.v = r.v * hj_inv.v;
  kernel_deval_vec(&xj, &wj, &wj_dx);
  wj_dr.v = hjd_inv.v * wj_dx.v;

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +
           ((vi[2].v - vj[2].v) * dx[2].v);
  // dvdr.v = dvdr.v * ri.v;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));
  mu_ij.v = fac_mu.v * ri.v * omega_ij.v; /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(3.0f) * mu_ij.v;

  /* Now construct the full viscosity term */
  rho_ij.v = vec_set1(0.5f) * (pirho.v + pjrho.v);
  visc.v = vec_set1(-0.25f) * vec_set1(const_viscosity_alpha) * v_sig.v *
           mu_ij.v * balsara.v / rho_ij.v;

  /* Now, convolve with the kernel */
  visc_term.v = vec_set1(0.5f) * visc.v * (wi_dr.v + wj_dr.v) * ri.v;
  sph_term.v =
      (grad_hi.v * piPOrho2.v * wi_dr.v + grad_hj.v * pjPOrho2.v * wj_dr.v) *
      ri.v;

  /* Eventually get the acceleration */
  acc.v = visc_term.v + sph_term.v;

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f.v = dx[k].v * acc.v;
    pia[k].v = mj.v * f.v;
  }

  /* Get the time derivative for h. */
  pih_dt.v = mj.v * dvdr.v * ri.v / pjrho.v * wi_dr.v;

  /* Change in entropy */
  entropy_dt.v = mj.v * visc_term.v * dvdr.v;

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    for (j = 0; j < 3; j++) pi[k]->a_hydro[j] -= pia[j].f[k];
    pi[k]->force.h_dt -= pih_dt.f[k];
    pi[k]->force.v_sig = max(pi[k]->force.v_sig, v_sig.f[k]);
    pi[k]->entropy_dt += entropy_dt.f[k];
  }

#else

  error(
      "The Gadget2 serial version of runner_iact_nonsym_force was called when "
      "the vectorised version should have been used.");

#endif
}

#endif /* SWIFT_GADGET2_HYDRO_IACT_H */
