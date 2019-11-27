/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_MFM_SLOPE_LIMITER_FACE_H
#define SWIFT_GIZMO_MFM_SLOPE_LIMITER_FACE_H

/* Some standard headers. */
#include <float.h>

/* Local headers. */
#include "minmax.h"
#include "sign.h"

__attribute__((always_inline)) INLINE static float
hydro_slope_limit_face_quantity(float phi_i, float phi_j, float phi_mid0,
                                float xij_norm, float r_inv) {

  const float psi1 = 0.5f;
  const float psi2 = 0.25f;

  const float delta1 = psi1 * fabsf(phi_i - phi_j);
  const float delta2 = psi2 * fabsf(phi_i - phi_j);

  const float phimin = min(phi_i, phi_j);
  const float phimax = max(phi_i, phi_j);

  const float phibar = phi_i + xij_norm * r_inv * (phi_j - phi_i);

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
 * @brief Slope limit the slopes at the interface between two particles
 *
 * @param Wi Hydrodynamic variables of particle i.
 * @param Wj Hydrodynamic variables of particle j.
 * @param dWi Difference between the hydrodynamic variables of particle i at the
 * position of particle i and at the interface position.
 * @param dWj Difference between the hydrodynamic variables of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float *Wi, float *Wj, float *dWi, float *dWj, const float *xij_i,
    const float *xij_j, float r) {

  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = 1.f / r;

  dWi[0] = hydro_slope_limit_face_quantity(Wi[0], Wj[0], Wi[0] + dWi[0],
                                           xij_i_norm, r_inv);
  dWi[1] = hydro_slope_limit_face_quantity(Wi[1], Wj[1], Wi[1] + dWi[1],
                                           xij_i_norm, r_inv);
  dWi[2] = hydro_slope_limit_face_quantity(Wi[2], Wj[2], Wi[2] + dWi[2],
                                           xij_i_norm, r_inv);
  dWi[3] = hydro_slope_limit_face_quantity(Wi[3], Wj[3], Wi[3] + dWi[3],
                                           xij_i_norm, r_inv);
  dWi[4] = hydro_slope_limit_face_quantity(Wi[4], Wj[4], Wi[4] + dWi[4],
                                           xij_i_norm, r_inv);

  dWj[0] = hydro_slope_limit_face_quantity(Wj[0], Wi[0], Wj[0] + dWj[0],
                                           xij_j_norm, r_inv);
  dWj[1] = hydro_slope_limit_face_quantity(Wj[1], Wi[1], Wj[1] + dWj[1],
                                           xij_j_norm, r_inv);
  dWj[2] = hydro_slope_limit_face_quantity(Wj[2], Wi[2], Wj[2] + dWj[2],
                                           xij_j_norm, r_inv);
  dWj[3] = hydro_slope_limit_face_quantity(Wj[3], Wi[3], Wj[3] + dWj[3],
                                           xij_j_norm, r_inv);
  dWj[4] = hydro_slope_limit_face_quantity(Wj[4], Wi[4], Wj[4] + dWj[4],
                                           xij_j_norm, r_inv);
}

#endif /* SWIFT_GIZMO_MFM_SLOPE_LIMITER_FACE_H */
