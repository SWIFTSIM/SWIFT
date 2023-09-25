//
// Created by yuyttenh on 8/06/23.
//

#ifndef SWIFTSIM_HYDRO_GRADIENTS_MESHLESS_H
#define SWIFTSIM_HYDRO_GRADIENTS_MESHLESS_H

#include "hydro_getters.h"
#include "hydro_part.h"
#include "hydro_slope_limiters.h"
#include "hydro_unphysical.h"
#define HYDRO_GRADIENT_IMPLEMENTATION "Meshless (Gizmo) gradients"

/**
 * @brief Add the gradient estimate for a single quantity due to a particle pair
 * to the total gradient for that quantity
 *
 * This corresponds to one term in the right hand side of eq. 26 in Pakmor 2014.
 *
 * @param qL Value of the quantity on the left.
 * @param qR Value of the quantity on the right.
 * @param w Relative weight assigned to this neighbour
 * @param cLR A vector pointing from the left centroid to the right centroid.
 * @param grad Current value of the gradient for the quantity (is updated).
 */
__attribute__((always_inline)) INLINE void hydro_gradients_single_quantity(
    float qL, float qR, float w, const float* restrict cLR,
    float* restrict grad) {

  float diff = qR - qL;
  grad[0] += w * diff * cLR[0];
  grad[1] += w * diff * cLR[1];
  grad[2] += w * diff * cLR[2];
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_collect(
    float r2, const float* dx, float hi, float hj, struct part* restrict pi,
    struct part* restrict pj) {

  /* Get ds (distance vector between centroids) */
  float ds[3];
  for (int i = 0; i < 3; i++)
    ds[i] = (float)-dx[i] + pj->geometry.centroid[i] - pi->geometry.centroid[i];

  /* Get weights */
  float wi, wj, wi_dx, wj_dx;
  const float r = sqrtf(r2);
  kernel_deval(r / hi, &wi, &wi_dx);
  kernel_deval(r / hj, &wj, &wj_dx);

  /* Update gradient estimates for pi */
  hydro_gradients_single_quantity(pi->rho, pj->rho, wi, ds, pi->gradients.rho);
  hydro_gradients_single_quantity(pi->v[0], pj->v[0], wi, ds,
                                  pi->gradients.v[0]);
  hydro_gradients_single_quantity(pi->v[1], pj->v[1], wi, ds,
                                  pi->gradients.v[1]);
  hydro_gradients_single_quantity(pi->v[2], pj->v[2], wi, ds,
                                  pi->gradients.v[2]);
  hydro_gradients_single_quantity(pi->P, pj->P, wi, ds, pi->gradients.P);
  hydro_gradients_single_quantity(pi->A, pj->A, wi, ds, pi->gradients.A);
  /* Update matrix and weight counter */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->gradients.matrix_wls[i][j] += wi * ds[i] * ds[j];
    }
  }

  /* Flip ds for pj */
  for (int i = 0; i < 3; i++) ds[i] = -ds[i];
  hydro_gradients_single_quantity(pj->rho, pi->rho, wj, ds, pj->gradients.rho);
  hydro_gradients_single_quantity(pj->v[0], pi->v[0], wj, ds,
                                  pj->gradients.v[0]);
  hydro_gradients_single_quantity(pj->v[1], pi->v[1], wj, ds,
                                  pj->gradients.v[1]);
  hydro_gradients_single_quantity(pj->v[2], pi->v[2], wj, ds,
                                  pj->gradients.v[2]);
  hydro_gradients_single_quantity(pj->P, pi->P, wj, ds, pj->gradients.P);
  hydro_gradients_single_quantity(pj->A, pi->A, wj, ds, pj->gradients.A);
  /* Update matrix and weight counter */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pj->gradients.matrix_wls[i][j] += wj * ds[i] * ds[j];
    }
  }

  /* Collect slope limiting info */
  hydro_slope_limit_cell_collect(pi, pj, r);
  hydro_slope_limit_cell_collect(pj, pi, r);
}

/**
 * @brief Gradient calculations done during the neighbour loop
 *
 * @param r2 Squared distance between the two particles.
 * @param dx Distance vector (pi->x - pj->x).
 * @param hi Smoothing length of particle i.
 * @param hj Smoothing length of particle j.
 * @param pi Particle i.
 * @param pj Particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_gradients_nonsym_collect(float r2, const float* dx, float hi, float hj,
                               struct part* restrict pi,
                               struct part* restrict pj) {

  /* Get ds (distance vector between centroids) */
  float ds[3];
  for (int i = 0; i < 3; i++)
    ds[i] = (float)-dx[i] + pj->geometry.centroid[i] - pi->geometry.centroid[i];

  /* Get weights */
  float wi, wi_dx;
  const float r = sqrtf(r2);
  kernel_deval(r / hi, &wi, &wi_dx);

  /* Update gradient estimates */
  hydro_gradients_single_quantity(pi->rho, pj->rho, wi, ds, pi->gradients.rho);
  hydro_gradients_single_quantity(pi->v[0], pj->v[0], wi, ds,
                                  pi->gradients.v[0]);
  hydro_gradients_single_quantity(pi->v[1], pj->v[1], wi, ds,
                                  pi->gradients.v[1]);
  hydro_gradients_single_quantity(pi->v[2], pj->v[2], wi, ds,
                                  pi->gradients.v[2]);
  hydro_gradients_single_quantity(pi->P, pj->P, wi, ds, pi->gradients.P);
  hydro_gradients_single_quantity(pi->A, pj->A, wi, ds, pi->gradients.A);
  /* Update matrix and weight counter */
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pi->gradients.matrix_wls[i][j] += wi * ds[i] * ds[j];
    }
  }

  /* Collect slope limiting info */
  hydro_slope_limit_cell_collect(pi, pj, r);
}

__attribute__((always_inline)) INLINE static void
hydro_gradients_finalize_single_quantity(float grad[3], const float B[3][3]) {
  float result[3];
  for (int i = 0; i < 3; i++) {
    result[i] = B[i][0] * grad[0] + B[i][1] * grad[1] + B[i][2] * grad[2];
  }
  grad[0] = result[0];
  grad[1] = result[1];
  grad[2] = result[2];
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part* p) {
  if (invert_dimension_by_dimension_matrix(p->gradients.matrix_wls) != 0) {
    warning(
        "Matrix inversion failed during gradient calculation! "
        "Falling back to first order for this particle!");
    hydro_gradients_init(p);
    return;
  }

  /* Now actually compute the gradient estimates */
  hydro_gradients_finalize_single_quantity(p->gradients.rho,
                                           p->gradients.matrix_wls);
  hydro_gradients_finalize_single_quantity(p->gradients.v[0],
                                           p->gradients.matrix_wls);
  hydro_gradients_finalize_single_quantity(p->gradients.v[1],
                                           p->gradients.matrix_wls);
  hydro_gradients_finalize_single_quantity(p->gradients.v[2],
                                           p->gradients.matrix_wls);
  hydro_gradients_finalize_single_quantity(p->gradients.P,
                                           p->gradients.matrix_wls);
  hydro_gradients_finalize_single_quantity(p->gradients.A,
                                           p->gradients.matrix_wls);

  /* Slope limit the gradients */
  hydro_slope_limit_cell(p);
}

/**
 * @brief Gradients reconstruction. Is the same for all gradient types (although
 * gradients_none does nothing, since all gradients are zero -- are they?).
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_predict(
    struct part* restrict pi, struct part* restrict pj, const float* dx,
    float r, const float* xij_i, float dt, float* Wi, float* Wj) {

  /* Compute the position of the face relative to the centroids of pi and pj */
  float dx_i[3] = {
      xij_i[0] - pi->geometry.centroid[0],
      xij_i[1] - pi->geometry.centroid[1],
      xij_i[2] - pi->geometry.centroid[2],
  };
  /* Compute interface position (relative to pj, since we don't need the actual
   * position) eqn. (8) */
  const float xij_j[3] = {xij_i[0] + dx[0], xij_i[1] + dx[1], xij_i[2] + dx[2]};
  float dx_j[3] = {
      xij_j[0] - pj->geometry.centroid[0],
      xij_j[1] - pj->geometry.centroid[1],
      xij_j[2] - pj->geometry.centroid[2],
  };

  float dWi[6], dWj[6];
  hydro_gradients_extrapolate(pi, dx_i, dWi);
  hydro_gradients_extrapolate(pj, dx_j, dWj);

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
  float ds[3] = {-dx[0] + pi->geometry.centroid[0] - pj->geometry.centroid[0],
                 -dx[1] + pi->geometry.centroid[1] - pj->geometry.centroid[1],
                 -dx[2] + pi->geometry.centroid[2] - pj->geometry.centroid[2]};
  float s = sqrtf(ds[0] * ds[0] + ds[1] * ds[1] + ds[2] * ds[2]);
  hydro_slope_limit_face(Wi, Wj, dWi, dWj, dx_i, dx_j, s);

  /* Apply the slope limited extrapolations */
  hydro_gradients_apply_extrapolation(Wi, dWi);
  hydro_gradients_apply_extrapolation(Wj, dWj);
}

#endif  // SWIFTSIM_HYDRO_GRADIENTS_MESHLESS_H
