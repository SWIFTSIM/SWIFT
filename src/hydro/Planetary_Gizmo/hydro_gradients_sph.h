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
#ifndef SWIFT_PLANETARY_GIZMO_HYDRO_SPH_GRADIENTS_H
#define SWIFT_PLANETARY_GIZMO_HYDRO_SPH_GRADIENTS_H

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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  float Wi[5], Wj[5];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  drho_i[0] = -dW[0] * nidx[0];
  drho_i[1] = -dW[0] * nidx[1];
  drho_i[2] = -dW[0] * nidx[2];

  dvx_i[0] = -dW[1] * nidx[0];
  dvx_i[1] = -dW[1] * nidx[1];
  dvx_i[2] = -dW[1] * nidx[2];
  dvy_i[0] = -dW[2] * nidx[0];
  dvy_i[1] = -dW[2] * nidx[1];
  dvy_i[2] = -dW[2] * nidx[2];
  dvz_i[0] = -dW[3] * nidx[0];
  dvz_i[1] = -dW[3] * nidx[1];
  dvz_i[2] = -dW[3] * nidx[2];

  dP_i[0] = -dW[4] * nidx[0];
  dP_i[1] = -dW[4] * nidx[1];
  dP_i[2] = -dW[4] * nidx[2];

  hydro_part_update_gradients(pi, drho_i, dvx_i, dvy_i, dvz_i, dP_i);

  hydro_slope_limit_cell_collect(pi, pj, r);

  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  const float normj = wj_dx * r_inv;
  const float njdx[3] = {normj * dx[0], normj * dx[1], normj * dx[2]};

  float drho_j[3], dvx_j[3], dvy_j[3], dvz_j[3], dP_j[3];

  /* signs are the same as before, since we swap i and j twice */
  drho_j[0] = -dW[0] * njdx[0];
  drho_j[1] = -dW[0] * njdx[1];
  drho_j[2] = -dW[0] * njdx[2];

  dvx_j[0] = -dW[1] * njdx[0];
  dvx_j[1] = -dW[1] * njdx[1];
  dvx_j[2] = -dW[1] * njdx[2];
  dvy_j[0] = -dW[2] * njdx[0];
  dvy_j[1] = -dW[2] * njdx[1];
  dvy_j[2] = -dW[2] * njdx[2];
  dvz_j[0] = -dW[3] * njdx[0];
  dvz_j[1] = -dW[3] * njdx[1];
  dvz_j[2] = -dW[3] * njdx[2];

  dP_j[0] = -dW[4] * njdx[0];
  dP_j[1] = -dW[4] * njdx[1];
  dP_j[2] = -dW[4] * njdx[2];

  hydro_part_update_gradients(pj, drho_j, dvx_j, dvy_j, dvz_j, dP_j);

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

  /* Get r and 1/r. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  float Wi[5], Wj[5];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  const float dW[5] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4]};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  float drho_i[3], dvx_i[3], dvy_i[3], dvz_i[3], dP_i[3];

  drho_i[0] = -dW[0] * nidx[0];
  drho_i[1] = -dW[0] * nidx[1];
  drho_i[2] = -dW[0] * nidx[2];

  dvx_i[0] = -dW[1] * nidx[0];
  dvx_i[1] = -dW[1] * nidx[1];
  dvx_i[2] = -dW[1] * nidx[2];
  dvy_i[0] = -dW[2] * nidx[0];
  dvy_i[1] = -dW[2] * nidx[1];
  dvy_i[2] = -dW[2] * nidx[2];
  dvz_i[0] = -dW[3] * nidx[0];
  dvz_i[1] = -dW[3] * nidx[1];
  dvz_i[2] = -dW[3] * nidx[2];

  dP_i[0] = -dW[4] * nidx[0];
  dP_i[1] = -dW[4] * nidx[1];
  dP_i[2] = -dW[4] * nidx[2];

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

  const float h = p->h;
  const float ih = 1.0f / h;
  const float ihdimp1 = pow_dimension_plus_one(ih);
  const float volume = p->geometry.volume;

  /* finalize gradients by multiplying with volume */
  hydro_part_normalise_gradients(p, ihdimp1 * volume);

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_PLANETARY_GIZMO_HYDRO_SPH_GRADIENTS_H */
