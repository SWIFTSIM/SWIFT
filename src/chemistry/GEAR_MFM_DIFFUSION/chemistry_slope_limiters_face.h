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
#ifndef SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_FACE_H
#define SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_FACE_H

/* Some standard headers. */
#include <float.h>

/* Local headers. */
#include "minmax.h"
#include "sign.h"

/**
 * @file src/chemistry/GEAR_MFM_FIDDUSION/chemistry_slope_limiters_face.h
 * @brief File containing routines concerning the face slope
 * limiter for the GEAR MFM diffusion scheme, i.e. limiting the actual
 * interparticle flux.
 *
 * */

/**
 * The minmod limiter.
 *
 * @param a Left slope
 * @param b Right slope
 */
__attribute__((always_inline)) INLINE static double chemistry_minmod(double a,
                                                                     double b) {

  if (a>0 && b>0) {
    return min(a,b);
  } else if (a < 0 && b < 0) {
    return max(a,b);
  } else {
    return 0.0;
  }
}

/**
 * The minmod slope limiter.
 *
 * @param dQi left slope
 * @param dQj right slope
 */
__attribute__((always_inline)) INLINE static void chemistry_limiter_minmod(
    double *dQi, double *dQj) {

  if (*dQi * *dQj > 0.f) {
    if (fabs(*dQi) < fabs(*dQj)) {
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
 * @brief Slope limit a single quantity at the interface using Gizmo
 * slope-limiter.
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
    phiplus = (phimax != 0.0f) ? phimax / (1.0f + delta1 / fabs(phimax)) : 0.0f;
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
 * @param Ui Chemistry variables of particle i.
 * @param Uj Chemistry variables of particle j.
 * @param dUi Difference between the chemistry variables of particle i at the
 * position of particle i and at the interface position.
 * @param dUj Difference between the chemistry variables of particle j at the
 * position of particle j and at the interface position.
 * @param xij_i Relative position vector of the interface w.r.t. particle i.
 * @param xij_j Relative position vector of the interface w.r.t. partilce j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void chemistry_slope_limit_face(
    double *Ui, double *Uj, double *dUi, double *dUj, const float xij_i[3],
    const float *xij_j, float r) {

  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  *dUi = chemistry_slope_limit_face_quantity(Ui[0], Uj[0], Ui[0] + dUi[0],
                                             xij_i_norm, r_inv);

  *dUj = chemistry_slope_limit_face_quantity(Uj[0], Ui[0], Uj[0] + dUj[0],
                                             xij_j_norm, r_inv);

  /* Test to see whether it is better to use the minmode or the GIZMO style
     limiter */
  /* chemistry_limiter_minmod(dUi, dUj); */
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
__attribute__((always_inline)) INLINE static float
chemistry_slope_limit_face_quantity_float(float phi_i, float phi_j,
                                          float phi_mid0, float xij_norm,
                                          float r_inv) {

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
    phiplus =
        (phimax != 0.0f) ? phimax / (1.0f + delta1 / fabsf(phimax)) : 0.0f;
  }

  if (same_signf(phimin - delta1, phimin)) {
    phiminus = phimin - delta1;
  } else {
    phiminus =
        (phimin != 0.0f) ? phimin / (1.0f + delta1 / fabsf(phimin)) : 0.0f;
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
__attribute__((always_inline)) INLINE static void
chemistry_slope_limit_face_hydro(float *Wi, float *Wj, float *dvi,
                                    float *dvj, const float *xij_i,
                                    const float *xij_j, float r) {

  const float xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);

  const float xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  const float r_inv = (r > 0.0f) ? 1.0f / r : 0.0f;

  dvi[0] = chemistry_slope_limit_face_quantity_float(
      Wi[1], Wj[1], Wi[1] + dvi[0], xij_i_norm, r_inv);
  dvi[1] = chemistry_slope_limit_face_quantity_float(
      Wi[2], Wj[2], Wi[2] + dvi[1], xij_i_norm, r_inv);
  dvi[2] = chemistry_slope_limit_face_quantity_float(
      Wi[3], Wj[3], Wi[3] + dvi[2], xij_i_norm, r_inv);

  dvj[0] = chemistry_slope_limit_face_quantity_float(
      Wj[1], Wi[1], Wj[1] + dvj[0], xij_j_norm, r_inv);
  dvj[1] = chemistry_slope_limit_face_quantity_float(
      Wj[2], Wi[2], Wj[2] + dvj[1], xij_j_norm, r_inv);
  dvj[2] = chemistry_slope_limit_face_quantity_float(
      Wj[3], Wi[3], Wj[3] + dvj[2], xij_j_norm, r_inv);
}

#endif /* SWIFT_CHEMISTRY_GEAR_MFM_DIFFUSION_SLOPE_LIMITERS_FACE_H */
