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
 * @brief Initialize gradient variables
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  p->primitives.gradients.rho[0] = 0.0f;
  p->primitives.gradients.rho[1] = 0.0f;
  p->primitives.gradients.rho[2] = 0.0f;

  p->primitives.gradients.v[0][0] = 0.0f;
  p->primitives.gradients.v[0][1] = 0.0f;
  p->primitives.gradients.v[0][2] = 0.0f;

  p->primitives.gradients.v[1][0] = 0.0f;
  p->primitives.gradients.v[1][1] = 0.0f;
  p->primitives.gradients.v[1][2] = 0.0f;

  p->primitives.gradients.v[2][0] = 0.0f;
  p->primitives.gradients.v[2][1] = 0.0f;
  p->primitives.gradients.v[2][2] = 0.0f;

  p->primitives.gradients.P[0] = 0.0f;
  p->primitives.gradients.P[1] = 0.0f;
  p->primitives.gradients.P[2] = 0.0f;

  hydro_slope_limit_cell_init(p);
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
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

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

  if (pi->density.wcorr > const_gizmo_min_wcorr) {
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

  } else {
    /* The gradient matrix was not well-behaved, switch to SPH gradients */

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
  }

  hydro_slope_limit_cell_collect(pi, pj, r);

  /* Compute kernel of pj. */
  hj_inv = 1.0 / hj;
  xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  if (pj->density.wcorr > const_gizmo_min_wcorr) {
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

  } else {
    /* SPH gradients */

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
  }

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
hydro_gradients_nonsym_collect(float r2, float *dx, float hi, float hj,
                               struct part *pi, struct part *pj) {

  float r = sqrtf(r2);
  float xi;
  float hi_inv;
  float wi, wi_dx;
  int k, l;
  float Bi[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (k = 0; k < 3; k++) {
    for (l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
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

  if (pi->density.wcorr > const_gizmo_min_wcorr) {
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

  } else {
    /* Gradient matrix is not well-behaved, switch to SPH gradients */

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
  }

  hydro_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  float h, ih;

  /* add kernel normalization to gradients */
  h = p->h;
  ih = 1.0f / h;
  const float ihdim = pow_dimension(ih);

  if (p->density.wcorr > const_gizmo_min_wcorr) {
    p->primitives.gradients.rho[0] *= ihdim;
    p->primitives.gradients.rho[1] *= ihdim;
    p->primitives.gradients.rho[2] *= ihdim;

    p->primitives.gradients.v[0][0] *= ihdim;
    p->primitives.gradients.v[0][1] *= ihdim;
    p->primitives.gradients.v[0][2] *= ihdim;
    p->primitives.gradients.v[1][0] *= ihdim;
    p->primitives.gradients.v[1][1] *= ihdim;
    p->primitives.gradients.v[1][2] *= ihdim;
    p->primitives.gradients.v[2][0] *= ihdim;
    p->primitives.gradients.v[2][1] *= ihdim;
    p->primitives.gradients.v[2][2] *= ihdim;

    p->primitives.gradients.P[0] *= ihdim;
    p->primitives.gradients.P[1] *= ihdim;
    p->primitives.gradients.P[2] *= ihdim;

  } else {
    const float ihdimp1 = pow_dimension_plus_one(ih);

    float volume = p->geometry.volume;

    /* finalize gradients by multiplying with volume */
    p->primitives.gradients.rho[0] *= ihdimp1 * volume;
    p->primitives.gradients.rho[1] *= ihdimp1 * volume;
    p->primitives.gradients.rho[2] *= ihdimp1 * volume;

    p->primitives.gradients.v[0][0] *= ihdimp1 * volume;
    p->primitives.gradients.v[0][1] *= ihdimp1 * volume;
    p->primitives.gradients.v[0][2] *= ihdimp1 * volume;

    p->primitives.gradients.v[1][0] *= ihdimp1 * volume;
    p->primitives.gradients.v[1][1] *= ihdimp1 * volume;
    p->primitives.gradients.v[1][2] *= ihdimp1 * volume;
    p->primitives.gradients.v[2][0] *= ihdimp1 * volume;
    p->primitives.gradients.v[2][1] *= ihdimp1 * volume;
    p->primitives.gradients.v[2][2] *= ihdimp1 * volume;

    p->primitives.gradients.P[0] *= ihdimp1 * volume;
    p->primitives.gradients.P[1] *= ihdimp1 * volume;
    p->primitives.gradients.P[2] *= ihdimp1 * volume;
  }

  hydro_slope_limit_cell(p);
}
