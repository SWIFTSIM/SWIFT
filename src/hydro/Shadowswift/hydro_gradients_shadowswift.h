//
// Created by yuyttenh on 20/04/22.
//

#ifndef SWIFTSIM_HYDRO_GRADIENTS_SHADOWSWIFT_H
#define SWIFTSIM_HYDRO_GRADIENTS_SHADOWSWIFT_H

#include "hydro_unphysical.h"
#define HYDRO_GRADIENT_IMPLEMENTATION "Gradients following Springel (2010)"

/* Forward declarations */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_face(
    float* Wi, float* Wj, float* dWi, float* dWj, const float* xij_i,
    const float* xij_j, float r);

/**
 * @brief Add the gradient estimate for a single quantity due to a particle pair
 * to the total gradient for that quantity
 *
 * This corresponds to one term of equation (21) in Springel (2010).
 *
 * @param qL Value of the quantity on the left.
 * @param qR Value of the quantity on the right.
 * @param cLR Vector pointing from the midpoint of the particle pair to the
 * geometrical centroid of the face in between the particles.
 * @param xLR Vector pointing from the right particle to the left particle.
 * @param rLR_inv Inverse distance between two particles.
 * @param A Surface area of the face in between the particles.
 * @param grad Current value of the gradient for the quantity (is updated).
 */
__attribute__((always_inline)) INLINE void hydro_gradients_single_quantity(
    float qL, float qR, const float* restrict cLR, const float* restrict xLR,
    float rLR_inv, float A, float* restrict grad) {

  grad[0] += A * ((qR - qL) * cLR[0] * rLR_inv - 0.5f * (qL + qR) * xLR[0] * rLR_inv);
  grad[1] += A * ((qR - qL) * cLR[1] * rLR_inv - 0.5f * (qL + qR) * xLR[1] * rLR_inv);
  grad[2] += A * ((qR - qL) * cLR[2] * rLR_inv - 0.5f * (qL + qR) * xLR[2] * rLR_inv);
}

/**
 * @brief Update the gradient estimation for a particle using a given
 * neighbouring particle.
 *
 * @param pi Particle we are updating.
 * @param pj Particle we are using to update pi.
 * @param cLR Vector pointing from the midpoint of the particle pair to the
 * geometrical centroid of the face in between the particles.
 * @param dx Vector pointing from pj to pi.
 * @param r Distance between pi and pj.
 * @param surface_area Surface area of the face between pi an pj.
 */
__attribute__((always_inline)) INLINE void hydro_slope_estimate_collect(
    struct part* restrict pi, const struct part* restrict pj,
    const float* restrict cLR, const float* restrict dx, float r,
    float surface_area) {

  const float r_inv = 1.f / r;
  hydro_gradients_single_quantity(pi->rho, pj->rho, cLR, dx, r_inv, surface_area,
                                  pi->gradients.rho);
  hydro_gradients_single_quantity(pi->v[0], pj->v[0], cLR, dx, r_inv, surface_area,
                                  pi->gradients.v[0]);
  hydro_gradients_single_quantity(pi->v[1], pj->v[1], cLR, dx, r_inv, surface_area,
                                  pi->gradients.v[1]);
  hydro_gradients_single_quantity(pi->v[2], pj->v[2], cLR, dx, r_inv, surface_area,
                                  pi->gradients.v[2]);
  hydro_gradients_single_quantity(pi->P, pj->P, cLR, dx, r_inv, surface_area,
                                  pi->gradients.P);
  hydro_gradients_single_quantity(pi->A, pj->A, cLR, dx, r_inv, surface_area,
                                  pi->gradients.A);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {
  const float volume_inv = 1.f / p->geometry.volume;

  p->gradients.rho[0] *= volume_inv;
  p->gradients.rho[1] *= volume_inv;
  p->gradients.rho[2] *= volume_inv;

  p->gradients.v[0][0] *= volume_inv;
  p->gradients.v[0][1] *= volume_inv;
  p->gradients.v[0][2] *= volume_inv;
  p->gradients.v[1][0] *= volume_inv;
  p->gradients.v[1][1] *= volume_inv;
  p->gradients.v[1][2] *= volume_inv;
  p->gradients.v[2][0] *= volume_inv;
  p->gradients.v[2][1] *= volume_inv;
  p->gradients.v[2][2] *= volume_inv;

  p->gradients.P[0] *= volume_inv;
  p->gradients.P[1] *= volume_inv;
  p->gradients.P[2] *= volume_inv;

  p->gradients.A[0] *= volume_inv;
  p->gradients.A[1] *= volume_inv;
  p->gradients.A[2] *= volume_inv;
}

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, const float* dx,
    float r, const float* xij_i, float dt, float* Wi, float* Wj) {

  /* perform gradient reconstruction in space and time */
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8) */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};

  float dWi[6], dWj[6];
  hydro_gradients_extrapolate(pi, xij_i, dWi);
  hydro_gradients_extrapolate(pj, xij_j, dWj);

#ifdef SHADOWSWIFT_EXTRAPOLATE_TIME
  /* Add the extrapolations in time, so they also get include in the slope
   * limiting. */
  dWi[0] += pi->dW_time[0];
  dWi[1] += pi->dW_time[1];
  dWi[2] += pi->dW_time[2];
  dWi[3] += pi->dW_time[3];
  dWi[4] += pi->dW_time[4];
  dWi[5] += pi->dW_time[5];

  dWj[0] += pj->dW_time[0];
  dWj[1] += pj->dW_time[1];
  dWj[2] += pj->dW_time[2];
  dWj[3] += pj->dW_time[3];
  dWj[4] += pj->dW_time[4];
  dWj[5] += pj->dW_time[5];
#endif

  /* Apply the slope limiter at this interface */
  hydro_slope_limit_face(Wi, Wj, dWi, dWj, xij_i, xij_j, r);

  /* Apply the slope limited extrapolations */
  hydro_gradients_apply_extrapolation(Wi, dWi);
  hydro_gradients_apply_extrapolation(Wj, dWj);
}

#endif  // SWIFTSIM_HYDRO_GRADIENTS_SHADOWSWIFT_H
