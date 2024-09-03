/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@epfl.ch)
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
#ifndef SWIFT_GEAR_CHEMISTRY_SLOPE_LIMITER_H
#define SWIFT_GEAR_CHEMISTRY_SLOPE_LIMITER_H


/* Some standard headers. */
#include <float.h>

/* Local headers. */
#include "minmax.h"
#include "sign.h"


/* Note: We do not need to slope limit the cell and the gradients, because we
   perform a first order reconstruction of nabla_otimes_q. If we were to use a
   first order reconstruction, then we would need the cell limiters. */

/**
 * The minmod slope limiter.
 *
 * @param a Left slope
 * @param b Right slope
 */
__attribute__((always_inline)) INLINE static double chemistry_limiter_minmod(double a, double b) {
  /* Write this more explicitely (taken from Gizmo) */
  return (a>0) ? ((b<0) ? 0 : min(a, b)) : ((b>=0) ? 0 : max(a, b));
}

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
__attribute__((always_inline)) INLINE static double
chemistry_slope_limit_face_quantity(double phi_i, double phi_j, double phi_mid0,
                                float xij_norm, float r_inv) {

  const double psi1 = 0.5;
  const double psi2 = 0.25;

  const double delta1 = psi1 * fabs(phi_i - phi_j);
  const double delta2 = psi2 * fabs(phi_i - phi_j);

  const double phimin = min(phi_i, phi_j);
  const double phimax = max(phi_i, phi_j);

  const double phibar = phi_i + xij_norm * r_inv * (phi_j - phi_i);

  double phiplus, phiminus, phi_mid;

  if (same_signf(phimax + delta1, phimax)) {
    phiplus = phimax + delta1;
  } else {
    phiplus =
        (phimax != 0.0f) ? phimax / (1.0f + delta1 / fabs(phimax)) : 0.0f;
  }

  if (same_signf(phimin - delta1, phimin)) {
    phiminus = phimin - delta1;
  } else {
    phiminus =
        (phimin != 0.0f) ? phimin / (1.0f + delta1 / fabs(phimin)) : 0.0f;
  }

  if (phi_i < phi_j) {
    const double temp = min(phibar + delta2, phi_mid0);
    phi_mid = max(phiminus, temp);
  } else {
    const double temp = max(phibar - delta2, phi_mid0);
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
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_face(
    double *Wi, double *Wj, double *dWi, double *dWj, const float xij_i[3],
    const float *xij_j, float r) {

  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  *dWi = chemistry_slope_limit_face_quantity(Wi[0], Wj[0], Wi[0] + dWi[0],
                                           xij_i_norm, r_inv);

  *dWj = chemistry_slope_limit_face_quantity(Wj[0], Wi[0], Wj[0] + dWj[0],
                                           xij_j_norm, r_inv);
}



#endif /*  SWIFT_GEAR_CHEMISTRY_SLOPE_LIMITER_H */
