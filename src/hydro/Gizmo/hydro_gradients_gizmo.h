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
    float r, int mode) {}

/**
 * @brief Gradient calculations done during the gradient loop
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_gradient_loop(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj,
    int mode) {

  float r = sqrtf(r2);
  float xi, xj;
  float hi_inv, hj_inv;
  float wi, wj, wi_dx, wj_dx;
  int k, l;
  float Bi[3][3];
  float Bj[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  Wi[0] = pi->primitives.rho;
  Wi[1] = pi->primitives.v[0];
  Wi[2] = pi->primitives.v[1];
  Wi[3] = pi->primitives.v[2];
  Wi[4] = pi->primitives.P;
  Wj[0] = pj->primitives.rho;
  Wj[1] = pj->primitives.v[0];
  Wj[2] = pj->primitives.v[1];
  Wj[3] = pj->primitives.v[2];
  Wj[4] = pj->primitives.P;

  /* Compute kernel of pi. */
  hi_inv = 1.0 / hi;
  xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  pi->primitives.gradients.rho[0] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.rho[1] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.rho[2] +=
      (Wi[0] - Wj[0]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.v[0][0] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[0][1] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[0][2] +=
      (Wi[1] - Wj[1]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[1][0] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[1][1] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[1][2] +=
      (Wi[2] - Wj[2]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  pi->primitives.gradients.v[2][0] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.v[2][1] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.v[2][2] +=
      (Wi[3] - Wj[3]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

  pi->primitives.gradients.P[0] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
  pi->primitives.gradients.P[1] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
  pi->primitives.gradients.P[2] +=
      (Wi[4] - Wj[4]) * wi *
      (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);

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
    /* Compute kernel of pj. */
    hj_inv = 1.0 / hj;
    xj = r * hj_inv;
    kernel_deval(xj, &wj, &wj_dx);

    /* Compute gradients for pj */
    /* there is no sign difference w.r.t. eqn. (6) because dx is now what we
     * want
     * it to be */
    pj->primitives.gradients.rho[0] +=
        (Wi[0] - Wj[0]) * wj *
        (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    pj->primitives.gradients.rho[1] +=
        (Wi[0] - Wj[0]) * wj *
        (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    pj->primitives.gradients.rho[2] +=
        (Wi[0] - Wj[0]) * wj *
        (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

    pj->primitives.gradients.v[0][0] +=
        (Wi[1] - Wj[1]) * wj *
        (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    pj->primitives.gradients.v[0][1] +=
        (Wi[1] - Wj[1]) * wj *
        (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    pj->primitives.gradients.v[0][2] +=
        (Wi[1] - Wj[1]) * wj *
        (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
    pj->primitives.gradients.v[1][0] +=
        (Wi[2] - Wj[2]) * wj *
        (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    pj->primitives.gradients.v[1][1] +=
        (Wi[2] - Wj[2]) * wj *
        (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    pj->primitives.gradients.v[1][2] +=
        (Wi[2] - Wj[2]) * wj *
        (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);
    pj->primitives.gradients.v[2][0] +=
        (Wi[3] - Wj[3]) * wj *
        (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    pj->primitives.gradients.v[2][1] +=
        (Wi[3] - Wj[3]) * wj *
        (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    pj->primitives.gradients.v[2][2] +=
        (Wi[3] - Wj[3]) * wj *
        (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

    pj->primitives.gradients.P[0] +=
        (Wi[4] - Wj[4]) * wj *
        (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    pj->primitives.gradients.P[1] +=
        (Wi[4] - Wj[4]) * wj *
        (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    pj->primitives.gradients.P[2] +=
        (Wi[4] - Wj[4]) * wj *
        (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

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
