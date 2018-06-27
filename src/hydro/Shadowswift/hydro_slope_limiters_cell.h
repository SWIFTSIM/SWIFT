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

#include <float.h>

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part* p) {

  p->primitives.limiter.rho[0] = FLT_MAX;
  p->primitives.limiter.rho[1] = -FLT_MAX;
  p->primitives.limiter.v[0][0] = FLT_MAX;
  p->primitives.limiter.v[0][1] = -FLT_MAX;
  p->primitives.limiter.v[1][0] = FLT_MAX;
  p->primitives.limiter.v[1][1] = -FLT_MAX;
  p->primitives.limiter.v[2][0] = FLT_MAX;
  p->primitives.limiter.v[2][1] = -FLT_MAX;
  p->primitives.limiter.P[0] = FLT_MAX;
  p->primitives.limiter.P[1] = -FLT_MAX;

  p->primitives.limiter.maxr = -FLT_MAX;
}

/**
 * @brief Collect information for the cell wide slope limiter during the
 * neighbour loop
 *
 * @param pi Particle i.
 * @param pj Particle j.
 * @param r Distance between particle i and particle j.
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_collect(struct part* pi, const struct part* pj,
                               float r) {

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->primitives.limiter.rho[0] =
      fmin(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      fmax(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      fmin(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      fmax(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      fmin(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      fmax(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      fmin(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      fmax(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      fmin(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      fmax(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.maxr = fmax(r, pi->primitives.limiter.maxr);
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
 * @param maxr Maximal distance between the generator and all of its neighbours
 */
__attribute__((always_inline)) INLINE static void
hydro_slope_limit_cell_quantity(float* grad, float qval, float qmin, float qmax,
                                float maxr) {

  float gradtrue, gradmax, gradmin, alpha;

  gradtrue = sqrtf(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
  if (gradtrue) {
    gradtrue *= maxr;
    gradmax = qmax - qval;
    gradmin = qval - qmin;
    alpha = fmin(1.0f, fmin(gradmax / gradtrue, gradmin / gradtrue));
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
    struct part* p) {

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.rho, p->primitives.rho,
      p->primitives.limiter.rho[0], p->primitives.limiter.rho[1],
      p->primitives.limiter.maxr);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[0], p->primitives.v[0],
      p->primitives.limiter.v[0][0], p->primitives.limiter.v[0][1],
      p->primitives.limiter.maxr);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[1], p->primitives.v[1],
      p->primitives.limiter.v[1][0], p->primitives.limiter.v[1][1],
      p->primitives.limiter.maxr);
  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.v[2], p->primitives.v[2],
      p->primitives.limiter.v[2][0], p->primitives.limiter.v[2][1],
      p->primitives.limiter.maxr);

  hydro_slope_limit_cell_quantity(
      p->primitives.gradients.P, p->primitives.P, p->primitives.limiter.P[0],
      p->primitives.limiter.P[1], p->primitives.limiter.maxr);
}
