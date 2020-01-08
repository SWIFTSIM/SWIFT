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
#ifndef SWIFT_GIZMO_HYDRO_SPH_GRADIENTS_H
#define SWIFT_GIZMO_HYDRO_SPH_GRADIENTS_H

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

  float wi, wi_dx;
  const float hi_inv = 1.0f / hi;
  const float xi = r * hi_inv;
  kernel_deval(xi, &wi, &wi_dx);

  float Wi[6], Wj[6];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  const float dW[6] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4], Wi[5] - Wj[5]};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  float gradWi[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWi[i][j] = -dW[i] * nidx[j];
    }
  }

  hydro_part_update_gradients(pi, gradWi[0], gradWi[1], gradWi[2], gradWi[3],
                              gradWi[4], gradWi[5]);

  hydro_slope_limit_cell_collect(pi, pj, r);

  float wj, wj_dx;
  const float hj_inv = 1.0f / hj;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);

  const float normj = wj_dx * r_inv;
  const float njdx[3] = {normj * dx[0], normj * dx[1], normj * dx[2]};

  /* signs are the same as before, since we swap i and j twice */
  float gradWj[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWj[i][j] = -dW[i] * njdx[j];
    }
  }

  hydro_part_update_gradients(pj, gradWj[0], gradWj[1], gradWj[2], gradWj[3],
                              gradWj[4], gradWj[5]);

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

  float Wi[6], Wj[6];
  hydro_part_get_primitive_variables(pi, Wi);
  hydro_part_get_primitive_variables(pj, Wj);

  const float dW[6] = {Wi[0] - Wj[0], Wi[1] - Wj[1], Wi[2] - Wj[2],
                       Wi[3] - Wj[3], Wi[4] - Wj[4], Wi[5] - Wj[5]};

  const float normi = wi_dx * r_inv;
  const float nidx[3] = {normi * dx[0], normi * dx[1], normi * dx[2]};

  float gradWi[6][3];

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 3; ++j) {
      gradWi[i][j] = -dW[i] * nidx[j];
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

  const float h = p->h;
  const float ih = 1.0f / h;
  const float ihdimp1 = pow_dimension_plus_one(ih);
  const float volume = p->geometry.volume;

  /* finalize gradients by multiplying with volume */
  hydro_part_normalise_gradients(p, ihdimp1 * volume);

  hydro_slope_limit_cell(p);
}

#endif /* SWIFT_GIZMO_HYDRO_SPH_GRADIENTS_H */
