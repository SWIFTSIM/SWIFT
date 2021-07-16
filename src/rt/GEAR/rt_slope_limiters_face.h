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
 * limiter for the GEAR RT scheme. (= second slope limiting
 * step done during actual particle interactions.) */

/**
 * @brief Slope limit a single quantity at the interface
 *
 * @param phi_i Value of the quantity at the particle position.
 * @param phi_j Value of the quantity at the neighbouring particle position.
 * @param phi_mid0 Extrapolated value of the quantity at the interface position.
 * @param xij_norm Distance between the particle position and the interface
 * position.
 * @param r Distance between the particle and its neighbour.
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
__attribute__((always_inline)) INLINE static float rt_limiter_minmod(
    const float dQi, const float dQj) {

  if (dQi * dQj > 0) {
    if (fabsf(dQi) < fabsf(dQj)) {
      return dQi;
    } else {
      return dQj;
    }
  } else {
    return 0.f;
  }
}

/**
 * the monotonized cenetral limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
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
 * @param dQi Difference between the RT quantities of particle i at the
 * position of particle i and at the interface position.
 * @param dQj Difference between the RT quantities of particle j at the
 * position of particle j and at the interface position.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_face(
    float dQi[4], float dQj[4]) {

  for (int i = 0; i < 4; i++) {
    /* Minmod and monotinized central difference limiters
     * give stable results in advection tests, superbee
     * and vanLeer limiters are unstable. */

    /* const float alpha = rt_limiter_minmod(dQi[i], dQj[i]); */
    const float alpha = rt_limiter_mc(dQi[i], dQj[i]);

    dQi[i] *= alpha;
    dQj[i] *= alpha;
  }
}

#endif /* SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H */
