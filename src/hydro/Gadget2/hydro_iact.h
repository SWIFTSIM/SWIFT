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
  pi->rho_dh -= mj * (3.f * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= ui * wi_dx;

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  pj->rho_dh -= mi * (3.f * wj + uj * wj_dx);

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

#ifdef WITH_VECTORIZATION

  vector r, ri, r2, xi, xj, hi, hj, hi_inv, hj_inv, wi, wj, wi_dx, wj_dx;
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
  /* Get each component of particle separation. (Dx={dx1,dy1,dz1,dx2,dy2,dz2,...,dxn,dyn,dzn})*/
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
  error("Unknown vector size.")
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  /*vec_rsqrt does not have the level of accuracy we need, so an extra term is added below.*/
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));
  xi.v = r.v * hi_inv.v;

  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);
  hj_inv.v = hj_inv.v - hj_inv.v * (hj_inv.v * hj.v - vec_set1(1.0f));
  xj.v = r.v * hj_inv.v;

  /* Compute the kernel function. */
  kernel_deval_vec(&xi, &wi, &wi_dx);
  kernel_deval_vec(&xj, &wj, &wj_dx);

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
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  /* Compute density of pj. */
  rhoj.v = mi.v * wj.v;
  rhoj_dh.v = mi.v * (vec_set1(3.0f) * wj.v + xj.v * wj_dx.v);
  wcountj.v = wj.v;
  wcountj_dh.v = xj.v * wj_dx.v;
  div_vj.v = mi.v * dvdr.v * wj_dx.v;
  for (k = 0; k < 3; k++) curl_vj[k].v = mi.v * curlvr[k].v * wj_dx.v;

  /* Update particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->div_v -= div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.rot_v[j] += curl_vi[j].f[k];
    pj[k]->rho += rhoj.f[k];
    pj[k]->rho_dh -= rhoj_dh.f[k];
    pj[k]->density.wcount += wcountj.f[k];
    pj[k]->density.wcount_dh -= wcountj_dh.f[k];
    pj[k]->div_v -= div_vj.f[k];
    for (j = 0; j < 3; j++) pj[k]->density.rot_v[j] += curl_vj[j].f[k];
  }

#else

    error("The Gadget2 serial version of runner_iact_density was called when the vectorised version should have been used.")

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
  const float ri = 1.0f / r;

  /* Compute the kernel function */
  const float h_inv = 1.0f / hi;
  const float u = r * h_inv;
  kernel_deval(u, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->rho_dh -= mj * (3.f * wi + u * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= u * wi_dx;

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

#ifdef WITH_VECTORIZATION

  vector r, ri, r2, xi, hi, hi_inv, wi, wi_dx;
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
  /* Get each component of particle separation. (Dx={dx1,dy1,dz1,dx2,dy2,dz2,...,dxn,dyn,dzn})*/
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
  error("Unknown vector size.")
#endif

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);
  /*vec_rsqrt does not have the level of accuracy we need, so an extra term is added below.*/
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));
  r.v = r2.v * ri.v;

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));
  xi.v = r.v * hi_inv.v;

  kernel_deval_vec(&xi, &wi, &wi_dx);

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
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;
  div_vi.v = mj.v * dvdr.v * wi_dx.v;
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;

  /* Update particles. */
  for (k = 0; k < VEC_SIZE; k++) {
    pi[k]->rho += rhoi.f[k];
    pi[k]->rho_dh -= rhoi_dh.f[k];
    pi[k]->density.wcount += wcounti.f[k];
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];
    pi[k]->div_v -= div_vi.f[k];
    for (j = 0; j < 3; j++) pi[k]->density.rot_v[j] += curl_vi[j].f[k];
  }

#else

    error("The Gadget2 serial version of runner_iact_nonsym_density was called when the vectorised version should have been used.")

#endif

}

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
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  const float balsara_i =
      fabsf(pi->div_v) /
      (fabsf(pi->div_v) + pi->force.curl_v + 0.0001f * ci / fac_mu / hi);
  const float balsara_j =
      fabsf(pj->div_v) /
      (fabsf(pj->div_v) + pj->force.curl_v + 0.0001f * cj / fac_mu / hj);

  /* Are the particles moving towards each others ? */
  const float omega_ij = fminf(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij *
                     (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term = (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

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
  pi->h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;
  pj->h_dt -= mi * dvdr * r_inv / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += 0.5f * mj * visc_term * dvdr;
  pj->entropy_dt -= 0.5f * mi * visc_term * dvdr;
}

/**
 * @brief Force loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {
  error(
      "A vectorised version of the Gadget2 force interaction function does not "
      "exist yet!");
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
  const float pressurei = pi->force.pressure;
  const float pressurej = pj->force.pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute sound speeds */
  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Balsara term */
  const float balsara_i =
      fabsf(pi->div_v) /
      (fabsf(pi->div_v) + pi->force.curl_v + 0.0001f * ci / fac_mu / hi);
  const float balsara_j =
      fabsf(pj->div_v) /
      (fabsf(pj->div_v) + pj->force.curl_v + 0.0001f * cj / fac_mu / hj);

  /* Are the particles moving towards each others ? */
  const float omega_ij = fminf(dvdr, 0.f);
  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

  /* Signal velocity */
  const float v_sig = ci + cj - 3.f * mu_ij;

  /* Now construct the full viscosity term */
  const float rho_ij = 0.5f * (rhoi + rhoj);
  const float visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij *
                     (balsara_i + balsara_j) / rho_ij;

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term = (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a_hydro[0] -= mj * acc * dx[0];
  pi->a_hydro[1] -= mj * acc * dx[1];
  pi->a_hydro[2] -= mj * acc * dx[2];

  /* Get the time derivative for h. */
  pi->h_dt -= mj * dvdr * r_inv / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += 0.5f * mj * visc_term * dvdr;
}

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {
  error(
      "A vectorised version of the Gadget2 non-symmetric force interaction "
      "function does not exist yet!");
}

#endif /* SWIFT_RUNNER_IACT_LEGACY_H */
