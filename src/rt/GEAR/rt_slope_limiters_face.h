/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H
#define SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H

#include "sign.h"

/**
 * @file src/rt/GEAR/rt_slope_limiters_face.h
 * @brief File containing routines concerning the face slope
 * limiter for the GEAR RT scheme, i.e. limiting the actual
 * interparticle flux. */

/**
 * @brief Slope limit a single quantity at the interface. CURRENTLY NOT IN USE.
 *
 * @param phi_i Value of the quantity at the particle position.
 * @param phi_j Value of the quantity at the neighbouring particle position.
 * @param phi_mid0 Extrapolated value of the quantity at the interface position.
 * @param xij_norm Distance between the particle position and the interface
 * position.
 * @param r_inv Inverse distance between the particle and its neighbour.
 * @return The slope limited difference between the quantity at the particle
 * position and the quantity at the interface position.
 */
__attribute__((always_inline)) INLINE static float rt_slope_limit_face_quantity(
    float phi_i, float phi_j, float phi_mid0, float xij_norm, float r_inv) {

  const float dphi = phi_i - phi_j;
  const float gamma1 = 0.5f;
  const float gamma2 = 0.25f;

  const float delta1 = gamma1 * fabsf(dphi);
  const float delta2 = gamma2 * fabsf(dphi);

  const float phimin = min(phi_i, phi_j);
  const float phimax = max(phi_i, phi_j);

  const float phibar = phi_i + xij_norm * r_inv * dphi;

  float phiplus, phiminus, phi_mid;

  if (same_signf(phimax + delta1, phimax)) {
    phiplus = phimax + delta1;
  } else {
    phiplus = phimax / (1.0f + delta1 / (fabsf(phimax) + FLT_MIN));
  }

  if (same_signf(phimin - delta1, phimin)) {
    phiminus = phimin - delta1;
  } else {
    phiminus = phimin / (1.0f + delta1 / (fabsf(phimin) + FLT_MIN));
  }

  if (phi_i < phi_j) {
    const float temp = min(phibar + delta2, phi_mid0);
    phi_mid = max(phiminus, temp);
  } else {
    const float temp = max(phibar - delta2, phi_mid0);
    phi_mid = min(phiplus, temp);
  }

  return phi_mid - phi_i;
}

/**
 * the minmod slope limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
 */
__attribute__((always_inline)) INLINE static void rt_limiter_minmod(
    float* dQi, float* dQj) {

  if (*dQi * *dQj > 0.f) {
    if (fabsf(*dQi) < fabsf(*dQj)) {
      *dQj = *dQi;
    } else {
      *dQi = *dQj;
    }
  } else {
    *dQi = 0.f;
    *dQj = 0.f;
  }
}

/**
 * the monotonized cenetral limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
 * @return factor to slope limit the slope dQi
 */
__attribute__((always_inline)) INLINE static float rt_limiter_mc(
    const float dQi, const float dQj) {

  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float minterm = min3(0.5 * (1. + r), 2.f, 2.f * r);
  return max(0.f, minterm);
}

/**
 * the van Leer limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
 * @return factor to slope limit the slope dQi
 */
__attribute__((always_inline)) INLINE static float rt_limiter_vanLeer(
    const float dQi, const float dQj) {
  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float absr = fabs(r);
  return (r + absr) / (1.f + absr);
}

/**
 * the superbee limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
 * @return factor to slope limit the slope dQi
 */
__attribute__((always_inline)) INLINE static float rt_limiter_superbee(
    const float dQi, const float dQj) {
  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float minterm1 = min(1.f, 2.f * r);
  const float minterm2 = min(2.f, r);
  return max3(0.f, minterm1, minterm2);
}

/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Qi RT quantities of particle i (energy + fluxes in 3 dim)
 * @param Qj RT quantities of particle j (energy + fluxes in 3 dim)
 * @param dQi Difference between the RT quantities of particle i at the
 * position of particle i and at the interface position.
 * @param dQj Difference between the RT quantities of particle j at the
 * position of particle j and at the interface position.
 * @param dx Comoving distance vector between the particles (dx = pi->x -
 * pj->x).
 * @param r Comoving distance between particle i and particle j.
 * @param xij_i Position of the "interface" w.r.t. position of particle i
 * @param xij_j Position of the "interface" w.r.t. position of particle j
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_face(
    const float Qi[4], const float Qj[4], float dQi[4], float dQj[4],
    const float* dx, const float r, const float xij_i[3],
    const float xij_j[3]) {

  /* In 1D advection tests, any limiter works with the
   * GLF solver. If you also activate (uncomment) the
   * cell slope limiting, you should avoid the superbee
   * and MC limiters. The Gizmo-style clever limiter
   * that is used for hydro leads to a more diffusive
   * behaviour.
   * For the HLL solver, avoid superbee and MC limiters,
   * as they lead to unstable results.
   */

  /* const float xij_i_norm = */
  /*     sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);
   */
  /*  */
  /* const float xij_j_norm = */
  /*     sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);
   */
  /*  */
  /* const float r_inv = 1.f / r; */

  for (int i = 0; i < 4; i++) {

    rt_limiter_minmod(&dQi[i], &dQj[i]);

    /* const float alphai = rt_limiter_mc(dQi[i], dQj[i]); */
    /* const float alphaj = rt_limiter_mc(dQj[i], dQi[i]); */
    /* dQi[i] *= alphai; */
    /* dQj[i] *= alphaj; */

    /* const float alphai = rt_limiter_superbee(dQi[i], dQj[i]); */
    /* const float alphaj = rt_limiter_superbee(dQj[i], dQi[i]); */
    /* dQi[i] *= alphai; */
    /* dQj[i] *= alphaj; */

    /* const float alphai = rt_limiter_vanLeer(dQi[i], dQj[i]); */
    /* const float alphaj = rt_limiter_vanLeer(dQj[i], dQi[i]); */
    /* dQi[i] *= alphai; */
    /* dQj[i] *= alphaj; */

    /* dQi[i] = rt_slope_limit_face_quantity(Qi[i], Qj[i], Qi[i] + dQi[i], */
    /*                                          xij_i_norm, r_inv); */
    /*  */
    /* dQj[i] = rt_slope_limit_face_quantity(Qj[i], Qi[i], Qj[i] + dQj[i], */
    /*                                        xij_j_norm, r_inv); */
  }
}

#endif /* SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H */
