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

  float wi, wi_dx;
  float wj, wj_dx;
  float dv[3], curlvr[3];

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mi = pi->mass;
  const float mj = pj->mass;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

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

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;

  if (pj->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pj->ids_ngbs_density[pj->num_ngb_density] = pi->id;
  ++pj->num_ngb_density;
#endif
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

  float wi, wi_dx;
  float dv[3], curlvr[3];

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

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

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_density < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_density[pi->num_ngb_density] = pj->id;
  ++pi->num_ngb_density;
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
                                 mask_t mask) {

  vector r, ri, ui, wi, wi_dx;
  vector dvx, dvy, dvz;
  vector dvdr;
  vector curlvrx, curlvry, curlvrz;

  /* Fill the vectors. */
  const vector mj = vector_load(Mj);
  const vector vjx = vector_load(Vjx);
  const vector vjy = vector_load(Vjy);
  const vector vjz = vector_load(Vjz);

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
                                 mask_t mask, mask_t mask2, int mask_cond) {

  vector r, ri, ui, wi, wi_dx;
  vector dvx, dvy, dvz;
  vector dvdr;
  vector curlvrx, curlvry, curlvrz;
  vector r_2, ri2, ui2, wi2, wi_dx2;
  vector dvx2, dvy2, dvz2;
  vector dvdr2;
  vector curlvrx2, curlvry2, curlvrz2;

  /* Fill the vectors. */
  const vector mj = vector_load(Mj);
  const vector mj2 = vector_load(&Mj[VEC_SIZE]);
  const vector vjx = vector_load(Vjx);
  const vector vjx2 = vector_load(&Vjx[VEC_SIZE]);
  const vector vjy = vector_load(Vjy);
  const vector vjy2 = vector_load(&Vjy[VEC_SIZE]);
  const vector vjz = vector_load(Vjz);
  const vector vjz2 = vector_load(&Vjz[VEC_SIZE]);
  const vector dx = vector_load(Dx);
  const vector dx2 = vector_load(&Dx[VEC_SIZE]);
  const vector dy = vector_load(Dy);
  const vector dy2 = vector_load(&Dy[VEC_SIZE]);
  const vector dz = vector_load(Dz);
  const vector dz2 = vector_load(&Dz[VEC_SIZE]);

  /* Get the radius and inverse radius. */
  const vector r2 = vector_load(R2);
  const vector r2_2 = vector_load(&R2[VEC_SIZE]);
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

  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

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

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

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
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr_Hubble;
  pj->entropy_dt += mi * visc_term * dvdr_Hubble;

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;

  if (pj->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pj->ids_ngbs_force[pj->num_ngb_force] = pi->id;
  ++pj->num_ngb_force;
#endif
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

  float wi, wj, wi_dx, wj_dx;

#ifdef SWIFT_DEBUG_CHECKS
  if (pi->time_bin >= time_bin_inhibited)
    error("Inhibited pi in interaction function!");
  if (pj->time_bin >= time_bin_inhibited)
    error("Inhibited pj in interaction function!");
#endif

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get r and 1/r. */
  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

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

  /* Add Hubble flow */
  const float dvdr_Hubble = dvdr + a2_Hubble * r2;

  /* Balsara term */
  const float balsara_i = pi->force.balsara;
  const float balsara_j = pj->force.balsara;

  /* Are the particles moving towards each others ? */
  const float omega_ij = min(dvdr_Hubble, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * v_sig * mu_ij * (balsara_i + balsara_j) / rho_ij;

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
  pi->force.v_sig = max(pi->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += mj * visc_term * dvdr_Hubble;

#ifdef DEBUG_INTERACTIONS_SPH
  /* Update ngb counters */
  if (pi->num_ngb_force < MAX_NUM_OF_NEIGHBOURS)
    pi->ids_ngbs_force[pi->num_ngb_force] = pj->id;
  ++pi->num_ngb_force;
#endif
}

#ifdef WITH_VECTORIZATION

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
    vector hj_inv, const float a, const float H, vector *a_hydro_xSum,
    vector *a_hydro_ySum, vector *a_hydro_zSum, vector *h_dtSum,
    vector *v_sigSum, vector *entropy_dtSum, mask_t mask) {

#ifdef WITH_VECTORIZATION

  vector r, ri;
  vector dvx, dvy, dvz;
  vector xi, xj;
  vector hid_inv, hjd_inv;
  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr, dvdr_Hubble;
  vector piax, piay, piaz;
  vector pih_dt;
  vector v_sig;
  vector omega_ij, mu_ij, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;

  /* Fill vectors. */
  const vector vjx = vector_load(Vjx);
  const vector vjy = vector_load(Vjy);
  const vector vjz = vector_load(Vjz);
  const vector mj = vector_load(Mj);
  const vector pjrho = vector_load(Pjrho);
  const vector grad_hj = vector_load(Grad_hj);
  const vector pjPOrho2 = vector_load(PjPOrho2);
  const vector balsara_j = vector_load(Balsara_j);
  const vector cj = vector_load(Cj);

  /* Cosmological terms */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;
  const vector v_fac_mu = vector_set1(fac_mu);
  const vector v_a2_Hubble = vector_set1(a2_Hubble);

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

  /* Calculate the kernel. */
  kernel_eval_dWdx_force_vec(&xj, &wj_dx);

  wj_dr.v = vec_mul(hjd_inv.v, wj_dx.v);

  /* Compute dv. */
  dvx.v = vec_sub(vix.v, vjx.v);
  dvy.v = vec_sub(viy.v, vjy.v);
  dvz.v = vec_sub(viz.v, vjz.v);

  /* Compute dv dot r. */
  dvdr.v = vec_fma(dvx.v, dx->v, vec_fma(dvy.v, dy->v, vec_mul(dvz.v, dz->v)));

  /* Add Hubble flow */
  dvdr_Hubble.v = vec_add(dvdr.v, vec_mul(v_a2_Hubble.v, r2->v));

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr_Hubble.v, vec_setzero());
  mu_ij.v = vec_mul(v_fac_mu.v,
                    vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v =
      vec_fnma(vec_set1(const_viscosity_beta), mu_ij.v, vec_add(ci.v, cj.v));

  /* Now construct the full viscosity term */
  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho.v));
  visc.v = vec_div(
      vec_mul(vec_set1(-0.25f), vec_mul(v_sig.v, vec_mul(mu_ij.v, balsara.v))),
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
  entropy_dt.v = vec_mul(mj.v, vec_mul(visc_term.v, dvdr_Hubble.v));

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
    float *Hj_inv, const float a, const float H, vector *a_hydro_xSum,
    vector *a_hydro_ySum, vector *a_hydro_zSum, vector *h_dtSum,
    vector *v_sigSum, vector *entropy_dtSum, mask_t mask, mask_t mask_2,
    short mask_cond) {

#ifdef WITH_VECTORIZATION

  vector r, ri;
  vector dvx, dvy, dvz;
  vector ui, uj;
  vector hid_inv, hjd_inv;
  vector wi_dx, wj_dx, wi_dr, wj_dr, dvdr, dvdr_Hubble;
  vector piax, piay, piaz;
  vector pih_dt;
  vector v_sig;
  vector omega_ij, mu_ij, balsara;
  vector rho_ij, visc, visc_term, sph_term, acc, entropy_dt;

  vector r_2, ri_2;
  vector dvx_2, dvy_2, dvz_2;
  vector ui_2, uj_2;
  vector hjd_inv_2;
  vector wi_dx_2, wj_dx_2, wi_dr_2, wj_dr_2, dvdr_2, dvdr_Hubble_2;
  vector piax_2, piay_2, piaz_2;
  vector pih_dt_2;
  vector v_sig_2;
  vector omega_ij_2, mu_ij_2, balsara_2;
  vector rho_ij_2, visc_2, visc_term_2, sph_term_2, acc_2, entropy_dt_2;

  /* Fill vectors. */
  const vector mj = vector_load(Mj);
  const vector mj_2 = vector_load(&Mj[VEC_SIZE]);
  const vector vjx = vector_load(Vjx);
  const vector vjx_2 = vector_load(&Vjx[VEC_SIZE]);
  const vector vjy = vector_load(Vjy);
  const vector vjy_2 = vector_load(&Vjy[VEC_SIZE]);
  const vector vjz = vector_load(Vjz);
  const vector vjz_2 = vector_load(&Vjz[VEC_SIZE]);
  const vector dx = vector_load(Dx);
  const vector dx_2 = vector_load(&Dx[VEC_SIZE]);
  const vector dy = vector_load(Dy);
  const vector dy_2 = vector_load(&Dy[VEC_SIZE]);
  const vector dz = vector_load(Dz);
  const vector dz_2 = vector_load(&Dz[VEC_SIZE]);

  /* Get the radius and inverse radius. */
  const vector r2 = vector_load(R2);
  const vector r2_2 = vector_load(&R2[VEC_SIZE]);
  ri = vec_reciprocal_sqrt(r2);
  ri_2 = vec_reciprocal_sqrt(r2_2);
  r.v = vec_mul(r2.v, ri.v);
  r_2.v = vec_mul(r2_2.v, ri_2.v);

  /* Get remaining properties. */
  const vector pjrho = vector_load(Pjrho);
  const vector pjrho_2 = vector_load(&Pjrho[VEC_SIZE]);
  const vector grad_hj = vector_load(Grad_hj);
  const vector grad_hj_2 = vector_load(&Grad_hj[VEC_SIZE]);
  const vector pjPOrho2 = vector_load(PjPOrho2);
  const vector pjPOrho2_2 = vector_load(&PjPOrho2[VEC_SIZE]);
  const vector balsara_j = vector_load(Balsara_j);
  const vector balsara_j_2 = vector_load(&Balsara_j[VEC_SIZE]);
  const vector cj = vector_load(Cj);
  const vector cj_2 = vector_load(&Cj[VEC_SIZE]);
  const vector hj_inv = vector_load(Hj_inv);
  const vector hj_inv_2 = vector_load(&Hj_inv[VEC_SIZE]);

  /* Cosmological terms */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;
  const vector v_fac_mu = vector_set1(fac_mu);
  const vector v_a2_Hubble = vector_set1(a2_Hubble);

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

  /* Add the Hubble flow */
  dvdr_Hubble.v = vec_add(dvdr.v, vec_mul(v_a2_Hubble.v, r2.v));
  dvdr_Hubble_2.v = vec_add(dvdr_2.v, vec_mul(v_a2_Hubble.v, r2_2.v));

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr_Hubble.v, vec_setzero());
  omega_ij_2.v = vec_fmin(dvdr_Hubble_2.v, vec_setzero());
  mu_ij.v = vec_mul(v_fac_mu.v,
                    vec_mul(ri.v, omega_ij.v)); /* This is 0 or negative */
  mu_ij_2.v = vec_mul(
      v_fac_mu.v, vec_mul(ri_2.v, omega_ij_2.v)); /* This is 0 or negative */

  /* Compute signal velocity */
  v_sig.v =
      vec_fnma(vec_set1(const_viscosity_beta), mu_ij.v, vec_add(ci.v, cj.v));
  v_sig_2.v = vec_fnma(vec_set1(const_viscosity_beta), mu_ij_2.v,
                       vec_add(ci.v, cj_2.v));

  /* Now construct the full viscosity term */
  rho_ij.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho.v));
  rho_ij_2.v = vec_mul(vec_set1(0.5f), vec_add(pirho.v, pjrho_2.v));

  visc.v = vec_div(
      vec_mul(vec_set1(-0.25f), vec_mul(v_sig.v, vec_mul(mu_ij.v, balsara.v))),
      rho_ij.v);
  visc_2.v =
      vec_div(vec_mul(vec_set1(-0.25f),
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
  entropy_dt.v = vec_mul(mj.v, vec_mul(visc_term.v, dvdr_Hubble.v));
  entropy_dt_2.v = vec_mul(mj_2.v, vec_mul(visc_term_2.v, dvdr_Hubble_2.v));

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

  /* Wake up the neighbour? */
  if (pi->force.v_sig > const_limiter_max_v_sig_ratio * pj->force.v_sig) {

    pj->wakeup = time_bin_awake;

    // MATTHIEU
    // if (pj->wakeup == time_bin_not_awake)
    // pj->wakeup = time_bin_awake;
    // else if (pj->wakeup > 0)
    // pj->wakeup = -pj->wakeup;
  }
}

#endif /* SWIFT_GADGET2_HYDRO_IACT_H */
