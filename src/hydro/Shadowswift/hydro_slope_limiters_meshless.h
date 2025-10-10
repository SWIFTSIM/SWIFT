//
// Created by yuyttenh on 8/06/23.
//

#ifndef SWIFTSIM_HYDRO_SLOPE_LIMITERS_MESHLESS_H
#define SWIFTSIM_HYDRO_SLOPE_LIMITERS_MESHLESS_H

/**
 * @brief Collect information for the cell wide slope limiter for a single
 * quantity.
 *
 * @param qR The value of the quantity for the neighbouring particle.
 * @param limiter This particles limiter for this quantity.
 * @param limiter_extrapolations This particles limiter for the extrapolations
 * of this quantity.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect_single(float qR, float *restrict limiter) {

  /* Collect the maximal and the minimal value for the primitive variables
   * among the ngbs */
  limiter[0] = fminf(qR, limiter[0]);
  limiter[1] = fmaxf(qR, limiter[1]);
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r distance between the centroids of pi and pj.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part *pi, struct part *pj, float r) {
  hydro_slope_limit_cell_collect_single(pj->rho, pi->limiter.rho);
  hydro_slope_limit_cell_collect_single(pj->v[0], pi->limiter.v[0]);
  hydro_slope_limit_cell_collect_single(pj->v[1], pi->limiter.v[1]);
  hydro_slope_limit_cell_collect_single(pj->v[2], pi->limiter.v[2]);
  hydro_slope_limit_cell_collect_single(pj->P, pi->limiter.P);
  hydro_slope_limit_cell_collect_single(pj->A, pi->limiter.A);
  pi->limiter.r_max = fmaxf(pi->limiter.r_max, r);
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
 * @param r Maximal distance to any of the neighbours
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_single(
    float *grad, float qval, float qmin, float qmax, float r) {
  float delta_max = qmax - qval;
  float delta_min = qmin - qval;
  float grad_nrm =
      sqrtf(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
  float emax = 0.5f * grad_nrm * r;
  float emin = -0.5f * grad_nrm * r;
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
  float r = p->limiter.r_max;
  hydro_slope_limit_cell_single(p->gradients.rho, p->rho, p->limiter.rho[0],
                                p->limiter.rho[1], r);
  hydro_slope_limit_cell_single(p->gradients.v[0], p->v[0], p->limiter.v[0][0],
                                p->limiter.v[0][1], r);
  hydro_slope_limit_cell_single(p->gradients.v[1], p->v[1], p->limiter.v[1][0],
                                p->limiter.v[1][1], r);
  hydro_slope_limit_cell_single(p->gradients.v[2], p->v[2], p->limiter.v[2][0],
                                p->limiter.v[2][1], r);
  hydro_slope_limit_cell_single(p->gradients.P, p->P, p->limiter.P[0],
                                p->limiter.P[1], r);
  hydro_slope_limit_cell_single(p->gradients.A, p->A, p->limiter.A[0],
                                p->limiter.A[1], r);
}

#endif  // SWIFTSIM_HYDRO_SLOPE_LIMITERS_MESHLESS_H
