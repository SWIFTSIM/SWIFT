//
// Created by yuyttenh on 29/03/22.
//

#ifndef SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
#define SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H

#include "hydro_gradients.h"

#include <float.h>

#define HYDRO_SLOPE_LIMITER_FACE_IMPLEMENTATION \
  "GIZMO piecewise slope limiter (Hopkins 2015)"
#define HYDRO_SLOPE_LIMITER_CELL_IMPLEMENTATION \
  "GIZMO Exact cell wide slope limiter (Hopkins 2015)"

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limiter_prepare(
    struct part *p) {
  p->limiter.rho[0] = FLT_MAX;
  p->limiter.rho[1] = -FLT_MAX;
  p->limiter.v[0][0] = FLT_MAX;
  p->limiter.v[0][1] = -FLT_MAX;
  p->limiter.v[1][0] = FLT_MAX;
  p->limiter.v[1][1] = -FLT_MAX;
  p->limiter.v[2][0] = FLT_MAX;
  p->limiter.v[2][1] = -FLT_MAX;
  p->limiter.P[0] = FLT_MAX;
  p->limiter.P[1] = -FLT_MAX;

  p->limiter.extrapolations.rho[0] = FLT_MAX;
  p->limiter.extrapolations.rho[1] = -FLT_MAX;
  p->limiter.extrapolations.v[0][0] = FLT_MAX;
  p->limiter.extrapolations.v[0][1] = -FLT_MAX;
  p->limiter.extrapolations.v[1][0] = FLT_MAX;
  p->limiter.extrapolations.v[1][1] = -FLT_MAX;
  p->limiter.extrapolations.v[2][0] = FLT_MAX;
  p->limiter.extrapolations.v[2][1] = -FLT_MAX;
  p->limiter.extrapolations.P[0] = FLT_MAX;
  p->limiter.extrapolations.P[1] = -FLT_MAX;
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param dx vector pointing from pi to the centroid of the face between pi and
 * pj.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj,
                               const float *dx) {

  /* Collect the maximal and the minimal value for the primitive variables
   * among the ngbs */
  pi->limiter.rho[0] = fminf(pj->rho, pi->limiter.rho[0]);
  pi->limiter.rho[1] = fmaxf(pj->rho, pi->limiter.rho[1]);

  pi->limiter.v[0][0] = fminf(pj->fluid_v[0], pi->limiter.v[0][0]);
  pi->limiter.v[0][1] = fmaxf(pj->fluid_v[0], pi->limiter.v[0][1]);
  pi->limiter.v[1][0] = fminf(pj->fluid_v[1], pi->limiter.v[1][0]);
  pi->limiter.v[1][1] = fmaxf(pj->fluid_v[1], pi->limiter.v[1][1]);
  pi->limiter.v[2][0] = fminf(pj->fluid_v[2], pi->limiter.v[2][0]);
  pi->limiter.v[2][1] = fmaxf(pj->fluid_v[2], pi->limiter.v[2][1]);

  pi->limiter.P[0] = fminf(pj->P, pi->limiter.P[0]);
  pi->limiter.P[1] = fmaxf(pj->P, pi->limiter.P[1]);

  /* Collect maximal and minimal extrapolated values of the primitive
   * quantities */
  float dW[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
  hydro_gradients_extrapolate(pi, dx, dW);

  pi->limiter.extrapolations.rho[0] =
      fminf(pi->limiter.extrapolations.rho[0], dW[0]);
  pi->limiter.extrapolations.rho[1] =
      fmaxf(pi->limiter.extrapolations.rho[1], dW[0]);

  pi->limiter.extrapolations.v[0][0] =
      fminf(dW[1], pi->limiter.extrapolations.v[0][0]);
  pi->limiter.extrapolations.v[0][1] =
      fmaxf(dW[1], pi->limiter.extrapolations.v[0][1]);
  pi->limiter.extrapolations.v[1][0] =
      fminf(dW[2], pi->limiter.extrapolations.v[1][0]);
  pi->limiter.extrapolations.v[1][1] =
      fmaxf(dW[2], pi->limiter.extrapolations.v[1][1]);
  pi->limiter.extrapolations.v[2][0] =
      fminf(dW[3], pi->limiter.extrapolations.v[2][0]);
  pi->limiter.extrapolations.v[2][1] =
      fmaxf(dW[3], pi->limiter.extrapolations.v[2][1]);

  pi->limiter.extrapolations.P[0] =
      fminf(pi->limiter.extrapolations.P[0], dW[4]);
  pi->limiter.extrapolations.P[1] =
      fmaxf(pi->limiter.extrapolations.P[1], dW[4]);
}

/**
 * @brief Apply the cell wide slope limiter to the gradient of a single quantity
 *
 * This corresponds to equation (B2) in Hopkins (2015).
 *
 * @param grad Gradient to slope limit
 * @param qval Value of the quantity at the cell generator
 * @param qmin Minimal value of the quantity among all cell neighbours
 * @param qmax Maximal value of the quantity among all cell neighbours
 * @param emax Maximal extrapolated value of the quantity among all cell
 *             neighbours
 * @param emin Minimal extrapolated value of the quantity among all cell
 *             neighbours
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_quantity(float *grad, float qval, float qmin, float qmax,
                                float emin, float emax) {
  float delta_max = qmax - qval;
  float delta_min = qmin - qval;
  float alpha = 1.f;
  if (emin != 0 && emax != 0) {
    alpha = fminf(1.0f, fminf(delta_max / emax, delta_min / emin));
  } else if (emin != 0) {
    alpha = fminf(1.0f, delta_min / emin);
  } else if (emax != 0) {
    alpha = fminf(1.0f, delta_max / emax);
  }
  if (alpha != 1.f) {
    grad[0] *= alpha;
    grad[1] *= alpha;
    grad[2] *= alpha;
  }
}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part *p) {
  hydro_slope_limit_cell_quantity(
      p->gradients.rho, p->rho, p->limiter.rho[0], p->limiter.rho[1],
      p->limiter.extrapolations.rho[0], p->limiter.extrapolations.rho[1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.v[0], p->fluid_v[0], p->limiter.v[0][0], p->limiter.v[0][1],
      p->limiter.extrapolations.v[0][0], p->limiter.extrapolations.v[0][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[1], p->fluid_v[1], p->limiter.v[1][0], p->limiter.v[1][1],
      p->limiter.extrapolations.v[1][0], p->limiter.extrapolations.v[1][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[2], p->fluid_v[2], p->limiter.v[2][0], p->limiter.v[2][1],
      p->limiter.extrapolations.v[2][0], p->limiter.extrapolations.v[2][1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.P, p->P, p->limiter.P[0], p->limiter.P[1],
      p->limiter.extrapolations.P[0], p->limiter.extrapolations.P[1]);
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
hydro_slope_limit_face_quantity(float phi_i, float phi_j, float phi_mid0,
                                float xij_norm, float r) {

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

  phibar = phi_i + xij_norm / r * (phi_j - phi_i);

  /* if sign(phimax+delta1) == sign(phimax) */
  if ((phimax + delta1) * phimax > 0.0f) {
    phiplus = phimax + delta1;
  } else {
    phiplus = phimax / (1.0f + delta1 / fabsf(phimax));
  }

  /* if sign(phimin-delta1) == sign(phimin) */
  if ((phimin - delta1) * phimin > 0.0f) {
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

  float xij_i_norm, xij_j_norm;
  xij_i_norm =
      sqrtf(xij_i[0] * xij_i[0] + xij_i[1] * xij_i[1] + xij_i[2] * xij_i[2]);
  xij_j_norm =
      sqrtf(xij_j[0] * xij_j[0] + xij_j[1] * xij_j[1] + xij_j[2] * xij_j[2]);

  dWi[0] = hydro_slope_limit_face_quantity(Wi[0], Wj[0], Wi[0] + dWi[0],
                                           xij_i_norm, r);
  dWi[1] = hydro_slope_limit_face_quantity(Wi[1], Wj[1], Wi[1] + dWi[1],
                                           xij_i_norm, r);
  dWi[2] = hydro_slope_limit_face_quantity(Wi[2], Wj[2], Wi[2] + dWi[2],
                                           xij_i_norm, r);
  dWi[3] = hydro_slope_limit_face_quantity(Wi[3], Wj[3], Wi[3] + dWi[3],
                                           xij_i_norm, r);
  dWi[4] = hydro_slope_limit_face_quantity(Wi[4], Wj[4], Wi[4] + dWi[4],
                                           xij_i_norm, r);

  dWj[0] = hydro_slope_limit_face_quantity(Wj[0], Wi[0], Wj[0] + dWj[0],
                                           xij_j_norm, r);
  dWj[1] = hydro_slope_limit_face_quantity(Wj[1], Wi[1], Wj[1] + dWj[1],
                                           xij_j_norm, r);
  dWj[2] = hydro_slope_limit_face_quantity(Wj[2], Wi[2], Wj[2] + dWj[2],
                                           xij_j_norm, r);
  dWj[3] = hydro_slope_limit_face_quantity(Wj[3], Wi[3], Wj[3] + dWj[3],
                                           xij_j_norm, r);
  dWj[4] = hydro_slope_limit_face_quantity(Wj[4], Wi[4], Wj[4] + dWj[4],
                                           xij_j_norm, r);
}

#endif  // SWIFTSIM_SHADOWSWIFT_HYDRO_SLOPE_LIMITERS_H
