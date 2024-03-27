//
// Created by yuyttenh on 20/04/22.
//

#ifndef SWIFTSIM_HYDRO_SLOPE_LIMITERS_PER_FACE_H
#define SWIFTSIM_HYDRO_SLOPE_LIMITERS_PER_FACE_H

/**
 * @brief Slope limit a single quantity at the interface
 *
 * @param phi_i Value of the quantity at the particle position.
 * @param phi_j Value of the quantity at the neighbouring particle position.
 * @param phi_mid0 Extrapolated value of the quantity at the interface position.
 * @param xij_norm_over_r Distance between the particle position and the
 * interface position divided by the distance between the particle and its
 * neighbour.
 * @return The slope limited difference between the quantity at the particle
 * position and the quantity at the interface position.
 */
__attribute__((always_inline)) INLINE static float
hydro_slope_limit_face_quantity(float phi_i, float phi_j, float phi_mid0,
                                float xij_norm_over_r) {

  float phi_mid, delta1, delta2, phimin, phimax, phibar, phiplus, phiminus;
  const float psi1 = 0.5f;
  const float psi2 = 0.25f;

  if (phi_i == phi_j) {
    return 0.f;
  }

  delta1 = psi1 * fabsf(phi_i - phi_j);
  delta2 = psi2 * fabsf(phi_i - phi_j);

  phimin = fminf(phi_i, phi_j);
  phimax = fmaxf(phi_i, phi_j);

  phibar = phi_i + xij_norm_over_r * (phi_j - phi_i);

  /* if sign(phimax+delta1) == sign(phimax) */
  if ((phimax + delta1) * phimax >= 0.0f) {
    phiplus = phimax + delta1;
  } else {
    phiplus = phimax / (1.0f + delta1 / fabsf(phimax));
  }

  /* if sign(phimin-delta1) == sign(phimin) */
  if ((phimin - delta1) * phimin >= 0.0f) {
    phiminus = phimin - delta1;
  } else {
    phiminus = phimin / (1.0f + delta1 / fabsf(phimin));
  }

  if (phi_i < phi_j) {
    phi_mid = fmaxf(phiminus, fminf(phibar + delta2, phi_mid0));
  } else {
    phi_mid = fminf(phiplus, fmaxf(phibar - delta2, phi_mid0));
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

  float xij_i_norm_over_r, xij_j_norm_over_r;
  const float r_inv = 1.f / r;
  xij_i_norm_over_r =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]) *
      r_inv;
  xij_j_norm_over_r =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]) *
      r_inv;

  dWi[0] = hydro_slope_limit_face_quantity(Wi[0], Wj[0], Wi[0] + dWi[0],
                                           xij_i_norm_over_r);
  dWi[1] = hydro_slope_limit_face_quantity(Wi[1], Wj[1], Wi[1] + dWi[1],
                                           xij_i_norm_over_r);
  dWi[2] = hydro_slope_limit_face_quantity(Wi[2], Wj[2], Wi[2] + dWi[2],
                                           xij_i_norm_over_r);
  dWi[3] = hydro_slope_limit_face_quantity(Wi[3], Wj[3], Wi[3] + dWi[3],
                                           xij_i_norm_over_r);
  dWi[4] = hydro_slope_limit_face_quantity(Wi[4], Wj[4], Wi[4] + dWi[4],
                                           xij_i_norm_over_r);
  dWi[5] = hydro_slope_limit_face_quantity(Wi[5], Wj[5], Wi[5] + dWi[5],
                                           xij_i_norm_over_r);

  dWj[0] = hydro_slope_limit_face_quantity(Wj[0], Wi[0], Wj[0] + dWj[0],
                                           xij_j_norm_over_r);
  dWj[1] = hydro_slope_limit_face_quantity(Wj[1], Wi[1], Wj[1] + dWj[1],
                                           xij_j_norm_over_r);
  dWj[2] = hydro_slope_limit_face_quantity(Wj[2], Wi[2], Wj[2] + dWj[2],
                                           xij_j_norm_over_r);
  dWj[3] = hydro_slope_limit_face_quantity(Wj[3], Wi[3], Wj[3] + dWj[3],
                                           xij_j_norm_over_r);
  dWj[4] = hydro_slope_limit_face_quantity(Wj[4], Wi[4], Wj[4] + dWj[4],
                                           xij_j_norm_over_r);
  dWj[5] = hydro_slope_limit_face_quantity(Wj[5], Wi[5], Wj[5] + dWj[5],
                                           xij_j_norm_over_r);
}

#endif  // SWIFTSIM_HYDRO_SLOPE_LIMITERS_PER_FACE_H
