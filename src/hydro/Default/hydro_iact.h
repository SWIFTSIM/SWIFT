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
#ifndef SWIFT_DEFAULT_HYDRO_IACT_H
#define SWIFT_DEFAULT_HYDRO_IACT_H

#include "adiabatic_index.h"

/**
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 * and use the same
 * numerical coefficients as the Gadget-2 code. When used with the Spline-3
 * kernel, the results
 * should be equivalent to the ones obtained with Gadget-2 up to the rounding
 * errors and interactions
 * missed by the Gadget-2 tree-code neighbours search.
 *
 * The code uses internal energy instead of entropy as a thermodynamical
 * variable.
 */

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

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float h_inv;
  float wi, wj, wi_dx, wj_dx;
  float mi, mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mi = pi->mass;
  mj = pj->mass;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  dvdr *= ri;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];
  for (k = 0; k < 3; k++) curlvr[k] *= ri;

  /* Compute density of pi. */
  h_inv = 1.0f / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (hydro_dimension * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  pi->density.div_v -= mj * dvdr * wi_dx;
  for (k = 0; k < 3; k++) pi->density.rot_v[k] += mj * curlvr[k] * wi_dx;

  /* Compute density of pj. */
  h_inv = 1.0f / hj;
  xj = r * h_inv;
  kernel_deval(xj, &wj, &wj_dx);

  pj->rho += mi * wj;
  pj->rho_dh -= mi * (hydro_dimension * wj + xj * wj_dx);
  pj->density.wcount += wj;
  pj->density.wcount_dh -= xj * wj_dx;

  pj->density.div_v -= mi * dvdr * wj_dx;
  for (k = 0; k < 3; k++) pj->density.rot_v[k] += mi * curlvr[k] * wj_dx;
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

  float r, ri;
  float xi;
  float h_inv;
  float wi, wi_dx;
  float mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mj = pj->mass;

  /* Get r and r inverse. */
  r = sqrtf(r2);
  ri = 1.0f / r;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  dvdr *= ri;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];
  for (k = 0; k < 3; k++) curlvr[k] *= ri;

  h_inv = 1.0f / hi;
  xi = r * h_inv;
  kernel_deval(xi, &wi, &wi_dx);

  pi->rho += mj * wi;
  pi->rho_dh -= mj * (hydro_dimension * wi + xi * wi_dx);
  pi->density.wcount += wi;
  pi->density.wcount_dh -= xi * wi_dx;

  pi->density.div_v -= mj * dvdr * wi_dx;
  for (k = 0; k < 3; k++) pi->density.rot_v[k] += mj * curlvr[k] * wi_dx;
}

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

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hid_inv;
  float hj_inv, hjd_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float mi, mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  float f;
  int k;

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get some values in local variables. */
  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.P_over_rho2;
  POrho2j = pj->force.P_over_rho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;
  hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;
  hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
         (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;
  dvdr *= ri;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = min(fac_mu * dvdr, 0.f);

  /* Compute signal velocity */
  v_sig = pi->force.soundspeed + pj->force.soundspeed -
          const_viscosity_beta * omega_ij;

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * hydro_gamma_minus_one *
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);
  tc *= (wi_dr + wj_dr);

  /* Get the common factor out. */
  w = ri *
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f = dx[k] * w;
    pi->a_hydro[k] -= mj * f;
    pj->a_hydro[k] += mi * f;
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));
  pj->force.u_dt +=
      mi * dvdr * (POrho2j * wj_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);
  pj->force.u_dt += mi * tc * (pj->u - pi->u);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;
  pj->force.h_dt -= mi * dvdr / rhoi * wj_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
  pj->force.v_sig = max(pj->force.v_sig, v_sig);
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

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hid_inv;
  float hj_inv, hjd_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float /*mi,*/ mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  float f;
  int k;

  /* Cosmological factors entering the EoMs */
  const float fac_mu = pow_three_gamma_minus_five_over_two(a);
  const float a2_Hubble = a * a * H;

  /* Get some values in local variables. */
  // mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.P_over_rho2;
  POrho2j = pj->force.P_over_rho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;
  hid_inv = pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);
  wi_dr = hid_inv * wi_dx;

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;
  hjd_inv = pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  wj_dr = hjd_inv * wj_dx;

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +
         (pi->v[2] - pj->v[2]) * dx[2] + a2_Hubble * r2;
  dvdr *= ri;

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = min(fac_mu * dvdr, 0.f);

  /* Compute signal velocity */
  v_sig = pi->force.soundspeed + pj->force.soundspeed -
          const_viscosity_beta * omega_ij;

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * hydro_gamma_minus_one *
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);
  tc *= (wi_dr + wj_dr);

  /* Get the common factor out. */
  w = ri *
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {
    f = dx[k] * w;
    pi->a_hydro[k] -= mj * f;
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);

  /* Get the time derivative for h. */
  pi->force.h_dt -= mj * dvdr / rhoj * wi_dr;

  /* Update the signal velocity. */
  pi->force.v_sig = max(pi->force.v_sig, v_sig);
}

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
  }
}

#endif /* SWIFT_DEFAULT_HYDRO_IACT_H */
