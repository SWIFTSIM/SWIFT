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
 * @brief Gradient calculations done during the density loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_density_loop(
    struct part *pi, struct part *pj, float wi_dx, float wj_dx, float *dx,
    float r, int mode) {

  /* very basic gradient estimate */
  pi->primitives.gradients.rho[0] -=
      wi_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[1] -=
      wi_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) / r;
  pi->primitives.gradients.rho[2] -=
      wi_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) / r;

  pi->primitives.gradients.v[0][0] -=
      wi_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][1] -=
      wi_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
  pi->primitives.gradients.v[0][2] -=
      wi_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;

  pi->primitives.gradients.v[1][0] -=
      wi_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][1] -=
      wi_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
  pi->primitives.gradients.v[1][2] -=
      wi_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;

  pi->primitives.gradients.v[2][0] -=
      wi_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][1] -=
      wi_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
  pi->primitives.gradients.v[2][2] -=
      wi_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;

  pi->primitives.gradients.P[0] -=
      wi_dx * dx[0] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[1] -=
      wi_dx * dx[1] * (pi->primitives.P - pj->primitives.P) / r;
  pi->primitives.gradients.P[2] -=
      wi_dx * dx[2] * (pi->primitives.P - pj->primitives.P) / r;

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

  if (mode == 1) {
    /* signs are the same as before, since we swap i and j twice */
    pj->primitives.gradients.rho[0] -=
        wj_dx * dx[0] * (pi->primitives.rho - pj->primitives.rho) / r;
    pj->primitives.gradients.rho[1] -=
        wj_dx * dx[1] * (pi->primitives.rho - pj->primitives.rho) / r;
    pj->primitives.gradients.rho[2] -=
        wj_dx * dx[2] * (pi->primitives.rho - pj->primitives.rho) / r;

    pj->primitives.gradients.v[0][0] -=
        wj_dx * dx[0] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
    pj->primitives.gradients.v[0][1] -=
        wj_dx * dx[1] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;
    pj->primitives.gradients.v[0][2] -=
        wj_dx * dx[2] * (pi->primitives.v[0] - pj->primitives.v[0]) / r;

    pj->primitives.gradients.v[1][0] -=
        wj_dx * dx[0] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
    pj->primitives.gradients.v[1][1] -=
        wj_dx * dx[1] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
    pj->primitives.gradients.v[1][2] -=
        wj_dx * dx[2] * (pi->primitives.v[1] - pj->primitives.v[1]) / r;
    pj->primitives.gradients.v[2][0] -=
        wj_dx * dx[0] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
    pj->primitives.gradients.v[2][1] -=
        wj_dx * dx[1] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;
    pj->primitives.gradients.v[2][2] -=
        wj_dx * dx[2] * (pi->primitives.v[2] - pj->primitives.v[2]) / r;

    pj->primitives.gradients.P[0] -=
        wj_dx * dx[0] * (pi->primitives.P - pj->primitives.P) / r;
    pj->primitives.gradients.P[1] -=
        wj_dx * dx[1] * (pi->primitives.P - pj->primitives.P) / r;
    pj->primitives.gradients.P[2] -=
        wj_dx * dx[2] * (pi->primitives.P - pj->primitives.P) / r;

    /* basic slope limiter: collect the maximal and the minimal value for the
     * primitive variables among the ngbs */
    pj->primitives.limiter.rho[0] =
        fmin(pi->primitives.rho, pj->primitives.limiter.rho[0]);
    pj->primitives.limiter.rho[1] =
        fmax(pi->primitives.rho, pj->primitives.limiter.rho[1]);
    pj->primitives.limiter.v[0][0] =
        fmin(pi->primitives.v[0], pj->primitives.limiter.v[0][0]);
    pj->primitives.limiter.v[0][1] =
        fmax(pi->primitives.v[0], pj->primitives.limiter.v[0][1]);
    pj->primitives.limiter.v[1][0] =
        fmin(pi->primitives.v[1], pj->primitives.limiter.v[1][0]);
    pj->primitives.limiter.v[1][1] =
        fmax(pi->primitives.v[1], pj->primitives.limiter.v[1][1]);
    pj->primitives.limiter.v[2][0] =
        fmin(pi->primitives.v[2], pj->primitives.limiter.v[2][0]);
    pj->primitives.limiter.v[2][1] =
        fmax(pi->primitives.v[2], pj->primitives.limiter.v[2][1]);

    pj->primitives.limiter.P[0] =
        fmin(pi->primitives.P, pj->primitives.limiter.P[0]);
    pj->primitives.limiter.P[1] =
        fmax(pi->primitives.P, pj->primitives.limiter.P[1]);

    pj->primitives.limiter.maxr = fmax(r, pj->primitives.limiter.maxr);
  }
}

/**
 * @brief Gradient calculations done during the gradient loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_gradient_loop(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {}
