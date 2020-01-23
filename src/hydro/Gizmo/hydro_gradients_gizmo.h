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
  float Wi[5], Wj[5];

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

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

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

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  drho_i[0] = dW[0] * wiBidx[0];
  drho_i[1] = dW[0] * wiBidx[1];
  drho_i[2] = dW[0] * wiBidx[2];

  dvx_i[0] = dW[1] * wiBidx[0];
  dvx_i[1] = dW[1] * wiBidx[1];
  dvx_i[2] = dW[1] * wiBidx[2];
  dvy_i[0] = dW[2] * wiBidx[0];
  dvy_i[1] = dW[2] * wiBidx[1];
  dvy_i[2] = dW[2] * wiBidx[2];
  dvz_i[0] = dW[3] * wiBidx[0];
  dvz_i[1] = dW[3] * wiBidx[1];
  dvz_i[2] = dW[3] * wiBidx[2];

  dP_i[0] = dW[4] * wiBidx[0];
  dP_i[1] = dW[4] * wiBidx[1];
  dP_i[2] = dW[4] * wiBidx[2];

  hydro_part_update_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);

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

  float drho_j[3], dvx_j[3], dvy_j[3], dvz_j[3], dP_j[3];

  /* Compute gradients for pj */
  /* there is no sign difference w.r.t. eqn. (6) because dx is now what we
   * want it to be */
  drho_j[0] = dW[0] * wjBjdx[0];
  drho_j[1] = dW[0] * wjBjdx[1];
  drho_j[2] = dW[0] * wjBjdx[2];

  dvx_j[0] = dW[1] * wjBjdx[0];
  dvx_j[1] = dW[1] * wjBjdx[1];
  dvx_j[2] = dW[1] * wjBjdx[2];
  dvy_j[0] = dW[2] * wjBjdx[0];
  dvy_j[1] = dW[2] * wjBjdx[1];
  dvy_j[2] = dW[2] * wjBjdx[2];
  dvz_j[0] = dW[3] * wjBjdx[0];
  dvz_j[1] = dW[3] * wjBjdx[1];
  dvz_j[2] = dW[3] * wjBjdx[2];

  dP_j[0] = dW[4] * wjBjdx[0];
  dP_j[1] = dW[4] * wjBjdx[1];
  dP_j[2] = dW[4] * wjBjdx[2];

  hydro_part_update_gradients(pj, drho_j, dvx_j, dvy_j, dvz_j, dP_j);

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
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  /* Compute kernel of pi. */
  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

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

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  /* Compute gradients for pi */
  /* there is a sign difference w.r.t. eqn. (6) because of the inverse
   * definition of dx */
  drho_i[0] = dW[0] * wiBidx[0];
  drho_i[1] = dW[0] * wiBidx[1];
  drho_i[2] = dW[0] * wiBidx[2];

  dvx_i[0] = dW[1] * wiBidx[0];
  dvx_i[1] = dW[1] * wiBidx[1];
  dvx_i[2] = dW[1] * wiBidx[2];
  dvy_i[0] = dW[2] * wiBidx[0];
  dvy_i[1] = dW[2] * wiBidx[1];
  dvy_i[2] = dW[2] * wiBidx[2];
  dvz_i[0] = dW[3] * wiBidx[0];
  dvz_i[1] = dW[3] * wiBidx[1];
  dvz_i[2] = dW[3] * wiBidx[2];

  dP_i[0] = dW[4] * wiBidx[0];
  dP_i[1] = dW[4] * wiBidx[1];
  dP_i[2] = dW[4] * wiBidx[2];

  hydro_part_update_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);

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
