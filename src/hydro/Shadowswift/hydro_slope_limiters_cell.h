//
// Created by yuyttenh on 20/04/22.
//

#ifndef SWIFTSIM_HYDRO_SLOPE_LIMITERS_CELL_WIDE_H
#define SWIFTSIM_HYDRO_SLOPE_LIMITERS_CELL_WIDE_H

__attribute__((always_inline)) INLINE static void hydro_gradients_extrapolate(
    const struct part *p, const float *dx, float *dW);

/**
 * @brief Collect information for the cell wide slope limiter for a single
 * quantity.
 *
 * @param qR The value of the quantity for the neighbouring particle.
 * @param eR The extrapolated difference of the quantity towards the \
 * neighbouring particle.
 * @param limiter This particles limiter for this quantity.
 * @param limiter_extrapolations This particles limiter for the extrapolations
 * of this quantity.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect_quantity(
    float qR, float eR, float *restrict limiter,
    float *restrict limiter_extrapolations) {

  /* Collect the maximal and the minimal value for the primitive variables
   * among the ngbs */
  limiter[0] = fminf(qR, limiter[0]);
  limiter[1] = fmaxf(qR, limiter[1]);

  /* Collect maximal and minimal extrapolated values of the primitive
   * quantities */
  limiter_extrapolations[0] = fminf(eR, limiter_extrapolations[0]);
  limiter_extrapolations[1] = fmaxf(eR, limiter_extrapolations[1]);
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
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj, float *dx) {

  /* Calculate extrapolations */
  float dW[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
#ifdef SHADOWSWIFT_GRADIENTS_WLS
  /* Extrapolate from the centroids of the cells */
  dx[0] -= pi->geometry.centroid[0];
  dx[1] -= pi->geometry.centroid[1];
  dx[2] -= pi->geometry.centroid[2];
#endif
  hydro_gradients_extrapolate(pi, dx, dW);

  hydro_slope_limit_cell_collect_quantity(pj->rho, dW[0], pi->limiter.rho,
                                          pi->limiter.extrapolations.rho);
  hydro_slope_limit_cell_collect_quantity(pj->v[0], dW[1], pi->limiter.v[0],
                                          pi->limiter.extrapolations.v[0]);
  hydro_slope_limit_cell_collect_quantity(pj->v[1], dW[2], pi->limiter.v[1],
                                          pi->limiter.extrapolations.v[1]);
  hydro_slope_limit_cell_collect_quantity(pj->v[2], dW[3], pi->limiter.v[2],
                                          pi->limiter.extrapolations.v[2]);
  hydro_slope_limit_cell_collect_quantity(pj->P, dW[4], pi->limiter.P,
                                          pi->limiter.extrapolations.P);
  hydro_slope_limit_cell_collect_quantity(pj->A, dW[5], pi->limiter.A,
                                          pi->limiter.extrapolations.A);
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
      p->gradients.v[0], p->v[0], p->limiter.v[0][0], p->limiter.v[0][1],
      p->limiter.extrapolations.v[0][0], p->limiter.extrapolations.v[0][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[1], p->v[1], p->limiter.v[1][0], p->limiter.v[1][1],
      p->limiter.extrapolations.v[1][0], p->limiter.extrapolations.v[1][1]);
  hydro_slope_limit_cell_quantity(
      p->gradients.v[2], p->v[2], p->limiter.v[2][0], p->limiter.v[2][1],
      p->limiter.extrapolations.v[2][0], p->limiter.extrapolations.v[2][1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.P, p->P, p->limiter.P[0], p->limiter.P[1],
      p->limiter.extrapolations.P[0], p->limiter.extrapolations.P[1]);

  hydro_slope_limit_cell_quantity(
      p->gradients.A, p->A, p->limiter.A[0], p->limiter.A[1],
      p->limiter.extrapolations.A[0], p->limiter.extrapolations.A[1]);
}

#endif  // SWIFTSIM_HYDRO_SLOPE_LIMITERS_CELL_WIDE_H
