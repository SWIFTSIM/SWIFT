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
#ifndef SWIFT_GIZMO_MFM_HYDRO_SPH_GRADIENTS_H
#define SWIFT_GIZMO_MFM_HYDRO_SPH_GRADIENTS_H

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

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {pi->rho - pj->rho, pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1], pi->v[2] - pj->v[2], pi->P - pj->P};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  pi->gradients.rho[0] -= dW[0] * nidx[0];
  pi->gradients.rho[1] -= dW[0] * nidx[1];
  pi->gradients.rho[2] -= dW[0] * nidx[2];

  pi->gradients.v[0][0] -= dW[1] * nidx[0];
  pi->gradients.v[0][1] -= dW[1] * nidx[1];
  pi->gradients.v[0][2] -= dW[1] * nidx[2];
  pi->gradients.v[1][0] -= dW[2] * nidx[0];
  pi->gradients.v[1][1] -= dW[2] * nidx[1];
  pi->gradients.v[1][2] -= dW[2] * nidx[2];
  pi->gradients.v[2][0] -= dW[3] * nidx[0];
  pi->gradients.v[2][1] -= dW[3] * nidx[1];
  pi->gradients.v[2][2] -= dW[3] * nidx[2];

  pi->gradients.P[0] -= dW[4] * nidx[0];
  pi->gradients.P[1] -= dW[4] * nidx[1];
  pi->gradients.P[2] -= dW[4] * nidx[2];

  hydro_slope_limit_cell_collect(pi, pj, r);

  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  const float normj = wj_dx * r_inv;
  const float njdx[3] = {normj * dx[0], normj * dx[1], normj * dx[2]};

  /* signs are the same as before, since we swap i and j twice */
  pj->gradients.rho[0] -= dW[0] * njdx[0];
  pj->gradients.rho[1] -= dW[0] * njdx[1];
  pj->gradients.rho[2] -= dW[0] * njdx[2];

  pj->gradients.v[0][0] -= dW[1] * njdx[0];
  pj->gradients.v[0][1] -= dW[1] * njdx[1];
  pj->gradients.v[0][2] -= dW[1] * njdx[2];
  pj->gradients.v[1][0] -= dW[2] * njdx[0];
  pj->gradients.v[1][1] -= dW[2] * njdx[1];
  pj->gradients.v[1][2] -= dW[2] * njdx[2];
  pj->gradients.v[2][0] -= dW[3] * njdx[0];
  pj->gradients.v[2][1] -= dW[3] * njdx[1];
  pj->gradients.v[2][2] -= dW[3] * njdx[2];

  pj->gradients.P[0] -= dW[4] * njdx[0];
  pj->gradients.P[1] -= dW[4] * njdx[1];
  pj->gradients.P[2] -= dW[4] * njdx[2];

  hydro_slope_limit_cell_collect(pj, pi, r);
}

/**
 * @brief Gradient calculations done during the neighbour loop: non-symmetric
 * version
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

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {pi->rho - pj->rho, pi->v[0] - pj->v[0],
                       pi->v[1] - pj->v[1], pi->v[2] - pj->v[2], pi->P - pj->P};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  pi->gradients.rho[0] -= dW[0] * nidx[0];
  pi->gradients.rho[1] -= dW[0] * nidx[1];
  pi->gradients.rho[2] -= dW[0] * nidx[2];

  pi->gradients.v[0][0] -= dW[1] * nidx[0];
  pi->gradients.v[0][1] -= dW[1] * nidx[1];
  pi->gradients.v[0][2] -= dW[1] * nidx[2];
  pi->gradients.v[1][0] -= dW[2] * nidx[0];
  pi->gradients.v[1][1] -= dW[2] * nidx[1];
  pi->gradients.v[1][2] -= dW[2] * nidx[2];
  pi->gradients.v[2][0] -= dW[3] * nidx[0];
  pi->gradients.v[2][1] -= dW[3] * nidx[1];
  pi->gradients.v[2][2] -= dW[3] * nidx[2];

  pi->gradients.P[0] -= dW[4] * nidx[0];
  pi->gradients.P[1] -= dW[4] * nidx[1];
  pi->gradients.P[2] -= dW[4] * nidx[2];

  hydro_slope_limit_cell_collect(pi, pj, r);
}

/**
 * @brief Finalize the gradient variables after all data have been collected
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_gradients_finalize(
    struct part *p) {

  const float h = p->h;
  const float ih = 1.0f / h;
  const float ihdimp1 = pow_dimension_plus_one(ih);
  const float volume = p->geometry.volume;

  const float norm = ihdimp1 * volume;

  /* finalize gradients by multiplying with volume */
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

#endif /* SWIFT_GIZMO_MFM_HYDRO_SPH_GRADIENTS_H */
