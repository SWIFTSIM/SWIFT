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
#ifndef SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H
#define SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H

__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  p->gradients.rho[0] = 0.0f;
  p->gradients.rho[1] = 0.0f;
  p->gradients.rho[2] = 0.0f;

  p->gradients.v[0][0] = 0.0f;
  p->gradients.v[0][1] = 0.0f;
  p->gradients.v[0][2] = 0.0f;

  p->gradients.v[1][0] = 0.0f;
  p->gradients.v[1][1] = 0.0f;
  p->gradients.v[1][2] = 0.0f;

  p->gradients.v[2][0] = 0.0f;
  p->gradients.v[2][1] = 0.0f;
  p->gradients.v[2][2] = 0.0f;

  p->gradients.P[0] = 0.0f;
  p->gradients.P[1] = 0.0f;
  p->gradients.P[2] = 0.0f;

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
    float r2, const float *dx, float hi, float hj, struct part *restrict pi,
    struct part *restrict pj) {

  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  float wi, wj, wi_dx, wj_dx;
  float Bi[3][3];
  float Bj[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  Wi[0] = pi->rho;
  Wi[1] = pi->v[0];
  Wi[2] = pi->v[1];
  Wi[3] = pi->v[2];
  Wi[4] = pi->P;
  Wj[0] = pj->rho;
  Wj[1] = pj->v[0];
  Wj[2] = pj->v[1];
  Wj[3] = pj->v[2];
  Wj[4] = pj->P;

  /* Compute kernel of pi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  float wiBidx[3];
  if (pi->geometry.wcorr > const_gizmo_min_wcorr) {
    wiBidx[0] = wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    wiBidx[1] = wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    wiBidx[2] = wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    wiBidx[0] = norm * dx[0];
    wiBidx[1] = norm * dx[1];
    wiBidx[2] = norm * dx[2];
  }

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  pi->gradients.rho[0] += dW[0] * wiBidx[0];
  pi->gradients.rho[1] += dW[0] * wiBidx[1];
  pi->gradients.rho[2] += dW[0] * wiBidx[2];

  pi->gradients.v[0][0] += dW[1] * wiBidx[0];
  pi->gradients.v[0][1] += dW[1] * wiBidx[1];
  pi->gradients.v[0][2] += dW[1] * wiBidx[2];
  pi->gradients.v[1][0] += dW[2] * wiBidx[0];
  pi->gradients.v[1][1] += dW[2] * wiBidx[1];
  pi->gradients.v[1][2] += dW[2] * wiBidx[2];
  pi->gradients.v[2][0] += dW[3] * wiBidx[0];
  pi->gradients.v[2][1] += dW[3] * wiBidx[1];
  pi->gradients.v[2][2] += dW[3] * wiBidx[2];

  pi->gradients.P[0] += dW[4] * wiBidx[0];
  pi->gradients.P[1] += dW[4] * wiBidx[1];
  pi->gradients.P[2] += dW[4] * wiBidx[2];

  hydro_slope_limit_cell_collect(pi, pj, r);

  /* Compute kernel of pj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  float wjBjdx[3];
  if (pj->geometry.wcorr > const_gizmo_min_wcorr) {

    wjBjdx[0] = wj * (Bj[0][0] * dx[0] + Bj[0][1] * dx[1] + Bj[0][2] * dx[2]);
    wjBjdx[1] = wj * (Bj[1][0] * dx[0] + Bj[1][1] * dx[1] + Bj[1][2] * dx[2]);
    wjBjdx[2] = wj * (Bj[2][0] * dx[0] + Bj[2][1] * dx[1] + Bj[2][2] * dx[2]);

  } else {
    const float norm = -wj_dx * r_inv;
    wjBjdx[0] = norm * dx[0];
    wjBjdx[1] = norm * dx[1];
    wjBjdx[2] = norm * dx[2];
  }

  /* Compute gradients for pj */
  /* there is no sign difference w.r.t. eqn. (6) because dx is now what we
   * want it to be */
  pj->gradients.rho[0] += dW[0] * wjBjdx[0];
  pj->gradients.rho[1] += dW[0] * wjBjdx[1];
  pj->gradients.rho[2] += dW[0] * wjBjdx[2];

  pj->gradients.v[0][0] += dW[1] * wjBjdx[0];
  pj->gradients.v[0][1] += dW[1] * wjBjdx[1];
  pj->gradients.v[0][2] += dW[1] * wjBjdx[2];
  pj->gradients.v[1][0] += dW[2] * wjBjdx[0];
  pj->gradients.v[1][1] += dW[2] * wjBjdx[1];
  pj->gradients.v[1][2] += dW[2] * wjBjdx[2];
  pj->gradients.v[2][0] += dW[3] * wjBjdx[0];
  pj->gradients.v[2][1] += dW[3] * wjBjdx[1];
  pj->gradients.v[2][2] += dW[3] * wjBjdx[2];

  pj->gradients.P[0] += dW[4] * wjBjdx[0];
  pj->gradients.P[1] += dW[4] * wjBjdx[1];
  pj->gradients.P[2] += dW[4] * wjBjdx[2];

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
hydro_gradients_nonsym_collect(float r2, const float *dx, float hi, float hj,
                               struct part *restrict pi,
                               struct part *restrict pj) {

  const float r_inv = 1.0f / sqrtf(r2);
  const float r = r2 * r_inv;

  float Bi[3][3];
  float Wi[5], Wj[5];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
    }
  }
  Wi[0] = pi->rho;
  Wi[1] = pi->v[0];
  Wi[2] = pi->v[1];
  Wi[3] = pi->v[2];
  Wi[4] = pi->P;
  Wj[0] = pj->rho;
  Wj[1] = pj->v[0];
  Wj[2] = pj->v[1];
  Wj[3] = pj->v[2];
  Wj[4] = pj->P;

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  float wiBidx[3];
  if (pi->geometry.wcorr > const_gizmo_min_wcorr) {
    wiBidx[0] = wi * (Bi[0][0] * dx[0] + Bi[0][1] * dx[1] + Bi[0][2] * dx[2]);
    wiBidx[1] = wi * (Bi[1][0] * dx[0] + Bi[1][1] * dx[1] + Bi[1][2] * dx[2]);
    wiBidx[2] = wi * (Bi[2][0] * dx[0] + Bi[2][1] * dx[1] + Bi[2][2] * dx[2]);
  } else {
    const float norm = -wi_dx * r_inv;
    wiBidx[0] = norm * dx[0];
    wiBidx[1] = norm * dx[1];
    wiBidx[2] = norm * dx[2];
  }

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  pi->gradients.rho[0] += dW[0] * wiBidx[0];
  pi->gradients.rho[1] += dW[0] * wiBidx[1];
  pi->gradients.rho[2] += dW[0] * wiBidx[2];

  pi->gradients.v[0][0] += dW[1] * wiBidx[0];
  pi->gradients.v[0][1] += dW[1] * wiBidx[1];
  pi->gradients.v[0][2] += dW[1] * wiBidx[2];
  pi->gradients.v[1][0] += dW[2] * wiBidx[0];
  pi->gradients.v[1][1] += dW[2] * wiBidx[1];
  pi->gradients.v[1][2] += dW[2] * wiBidx[2];
  pi->gradients.v[2][0] += dW[3] * wiBidx[0];
  pi->gradients.v[2][1] += dW[3] * wiBidx[1];
  pi->gradients.v[2][2] += dW[3] * wiBidx[2];

  pi->gradients.P[0] += dW[4] * wiBidx[0];
  pi->gradients.P[1] += dW[4] * wiBidx[1];
  pi->gradients.P[2] += dW[4] * wiBidx[2];

  hydro_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  /* add kernel normalization to gradients */
  const float volume = p->geometry.volume;
  const float h = p->h;
  const float h_inv = 1.0f / h;
  const float ihdim = pow_dimension(h_inv);

  float norm;
  if (p->geometry.wcorr > const_gizmo_min_wcorr) {
    norm = ihdim;
  } else {
    const float ihdimp1 = pow_dimension_plus_one(h_inv);
    norm = ihdimp1 * volume;
  }

  p->gradients.rho[0] *= norm;
  p->gradients.rho[1] *= norm;
  p->gradients.rho[2] *= norm;

  p->gradients.v[0][0] *= norm;
  p->gradients.v[0][1] *= norm;
  p->gradients.v[0][2] *= norm;
  p->gradients.v[1][0] *= norm;
  p->gradients.v[1][1] *= norm;
  p->gradients.v[1][2] *= norm;
  p->gradients.v[2][0] *= norm;
  p->gradients.v[2][1] *= norm;
  p->gradients.v[2][2] *= norm;

  p->gradients.P[0] *= norm;
  p->gradients.P[1] *= norm;
  p->gradients.P[2] *= norm;

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H */
