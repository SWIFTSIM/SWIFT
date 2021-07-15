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
__attribute__((always_inline)) INLINE static float
rt_slope_limit_face_quantity(float phi_i, float phi_j, float phi_mid0,
                                float xij_norm, float r_inv) {

  const float dphi = phi_i - phi_j;
  /* TODO: determine good factors for the limiter depending
   * on different problems and reduction factors for the
   * speed of light. For now, choose stability over diffusivity. */

  /* const float gamma1 = 1.0f; */
  /* const float gamma2 = 0.5f; */
  const float gamma1 = 0.5f;
  const float gamma2 = 0.25f;
  /* const float gamma1 = 0.1f; */
  /* const float gamma2 = 0.05f; */
  /* const float gamma1 = 0.f; */
  /* const float gamma2 = 0.f; */

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

__attribute__((always_inline)) INLINE static float
rt_minmod(const float a, const float b) {

  if (a * b > 0) {
    if (fabsf(a) < fabsf(b)) {
      return a;
    } else {
      return b;
    }
  } else {
    return 0.f;
  }
}

/* monotonized cenetral limiter */
__attribute__((always_inline)) INLINE static float rt_mc(const float dQi, const float dQj) {

  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float minterm = min3(0.5 * (1. + r), 2.f, 2.f * r);
  return max(0.f, minterm);
}

/* van Leer limiter */
__attribute__((always_inline)) INLINE static float rt_vanLeer(const float dQi, const float dQj) {
  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float absr = fabs(r);
  return (r + absr) / (1.f + absr);
}


/* superbee */
__attribute__((always_inline)) INLINE static float rt_superbee(const float dQi, const float dQj) {
  const float r = dQj == 0.f ? dQi * 1e6 : dQi / dQj;
  const float minterm1 = min(1.f, 2.f*r);
  const float minterm2 = min(2.f, r);
  return max3(0.f, minterm1, minterm2);
}








/**
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Qi RT quantities (photon energies and flux) of particle i.
 * @param Qj RT quantities (photon energies and flux) of particle j.
 * @param dQi Difference between the RT quantities of particle i at the
 * position of particle i and at the interface position.
 * @param dQj Difference between the RT quantities of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void rt_slope_limit_face(
    const float Qi[4], const float Qj[4], float dQi[4], float dQj[4],
    const float xij_i[3], const float xij_j[3], const float r) {

  /* const float xij_i_norm = */
  /*     sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]); */
  /*  */
  /* const float xij_j_norm = */
  /*     sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]); */
  /*  */
  /* const float r_inv = 1.f / r; */

  for (int i = 0; i < 4; i++) {
    /* dQi[i] = rt_slope_limit_face_quantity(Qi[i], Qj[i], Qi[i] + dQi[i], xij_i_norm, r_inv); */
    /* dQj[i] = rt_slope_limit_face_quantity(Qj[i], Qi[i], Qj[i] + dQj[i], xij_j_norm, r_inv); */

    /* const float alpha = rt_minmod(dQi[i], dQj[i]); */
    const float alpha = rt_mc(dQi[i], dQj[i]);
    /* const float alpha = rt_vanLeer(dQi[i], dQj[i]); */
    /* const float alpha = rt_superbee(dQi[i], dQj[i]); */

    dQi[i] *= alpha;
    dQj[i] *= alpha;
  }

}

#endif /* SWIFT_RT_SLOPE_LIMITERS_FACE_GEAR_H */
