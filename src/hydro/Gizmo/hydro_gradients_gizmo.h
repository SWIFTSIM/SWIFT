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

#include "hydro_getters.h"
#include "hydro_setters.h"

__attribute__((always_inline)) INLINE static void hydro_gradients_init(
    struct part *p) {

  hydro_part_reset_gradients(p);

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
  float Wi[6], Wj[6];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
      Bj[k][l] = pj->geometry.matrix_E[k][l];
    }
  }
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* Compute kernel of pi. */
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[6] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4], Wi[5] - Wj[5]};

  float wiBidx[3];
  if (hydro_part_geometry_well_behaved(pi)) {
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

  float gradWi[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWi[i][j] = dW[i] * wiBidx[j];
    }
  }

  hydro_part_update_gradients(pi, gradWi[0], gradWi[1], gradWi[2], gradWi[3],
                              gradWi[4], gradWi[5]);

  hydro_slope_limit_cell_collect(pi, pj, r);

  /* Compute kernel of pj. */
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  float wjBjdx[3];
  if (hydro_part_geometry_well_behaved(pj)) {

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
  float gradWj[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWj[i][j] = dW[i] * wjBjdx[j];
    }
  }

  hydro_part_update_gradients(pj, gradWj[0], gradWj[1], gradWj[2], gradWj[3],
                              gradWj[4], gradWj[5]);

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
  float Wi[6], Wj[6];

  /* Initialize local variables */
  for (int k = 0; k < 3; k++) {
    for (int l = 0; l < 3; l++) {
      Bi[k][l] = pi->geometry.matrix_E[k][l];
    }
  }
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[6] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4], Wi[5] - Wj[5]};

  float wiBidx[3];
  if (hydro_part_geometry_well_behaved(pi)) {
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
  float gradWi[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWi[i][j] = dW[i] * wiBidx[j];
    }
  }

  hydro_part_update_gradients(pi, gradWi[0], gradWi[1], gradWi[2], gradWi[3],
                              gradWi[4], gradWi[5]);

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
  if (hydro_part_geometry_well_behaved(p)) {
    norm = ihdim;
  } else {
    const float ihdimp1 = pow_dimension_plus_one(h_inv);
    norm = ihdimp1 * volume;
  }

  hydro_part_normalise_gradients(p, norm);

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_GIZMO_MFM_HYDRO_GRADIENTS_H */
