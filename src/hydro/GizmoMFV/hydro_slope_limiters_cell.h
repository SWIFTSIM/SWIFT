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
#ifndef SWIFT_GIZMO_MFV_SLOPE_LIMITER_CELL_H
#define SWIFT_GIZMO_MFV_SLOPE_LIMITER_CELL_H

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
  p->primitives.limiter.A[0] = FLT_MAX;
  p->primitives.limiter.A[1] = -FLT_MAX;

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
hydro_slope_limit_cell_collect(struct part* pi, struct part* pj, float r) {

  /* basic slope limiter: collect the maximal and the minimal value for the
   * primitive variables among the ngbs */
  pi->primitives.limiter.rho[0] =
      min(pj->primitives.rho, pi->primitives.limiter.rho[0]);
  pi->primitives.limiter.rho[1] =
      max(pj->primitives.rho, pi->primitives.limiter.rho[1]);

  pi->primitives.limiter.v[0][0] =
      min(pj->primitives.v[0], pi->primitives.limiter.v[0][0]);
  pi->primitives.limiter.v[0][1] =
      max(pj->primitives.v[0], pi->primitives.limiter.v[0][1]);
  pi->primitives.limiter.v[1][0] =
      min(pj->primitives.v[1], pi->primitives.limiter.v[1][0]);
  pi->primitives.limiter.v[1][1] =
      max(pj->primitives.v[1], pi->primitives.limiter.v[1][1]);
  pi->primitives.limiter.v[2][0] =
      min(pj->primitives.v[2], pi->primitives.limiter.v[2][0]);
  pi->primitives.limiter.v[2][1] =
      max(pj->primitives.v[2], pi->primitives.limiter.v[2][1]);

  pi->primitives.limiter.P[0] =
      min(pj->primitives.P, pi->primitives.limiter.P[0]);
  pi->primitives.limiter.P[1] =
      max(pj->primitives.P, pi->primitives.limiter.P[1]);

  pi->primitives.limiter.A[0] =
      min(pj->primitives.A, pi->primitives.limiter.A[0]);
  pi->primitives.limiter.A[1] =
      max(pj->primitives.A, pi->primitives.limiter.A[1]);

  pi->primitives.limiter.maxr = max(r, pi->primitives.limiter.maxr);
}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part* p) {

  float gradtrue, gradrho[3], gradv[3][3], gradP[3], gradA[3];

  gradrho[0] = p->primitives.gradients.rho[0];
  gradrho[1] = p->primitives.gradients.rho[1];
  gradrho[2] = p->primitives.gradients.rho[2];

  gradv[0][0] = p->primitives.gradients.v[0][0];
  gradv[0][1] = p->primitives.gradients.v[0][1];
  gradv[0][2] = p->primitives.gradients.v[0][2];

  gradv[1][0] = p->primitives.gradients.v[1][0];
  gradv[1][1] = p->primitives.gradients.v[1][1];
  gradv[1][2] = p->primitives.gradients.v[1][2];

  gradv[2][0] = p->primitives.gradients.v[2][0];
  gradv[2][1] = p->primitives.gradients.v[2][1];
  gradv[2][2] = p->primitives.gradients.v[2][2];

  gradP[0] = p->primitives.gradients.P[0];
  gradP[1] = p->primitives.gradients.P[1];
  gradP[2] = p->primitives.gradients.P[2];

  gradA[0] = p->primitives.gradients.A[0];
  gradA[1] = p->primitives.gradients.A[1];
  gradA[2] = p->primitives.gradients.A[2];

  gradtrue = sqrtf(gradrho[0] * gradrho[0] + gradrho[1] * gradrho[1] +
                   gradrho[2] * gradrho[2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.rho[1] - p->primitives.rho;
    const float gradmin = p->primitives.rho - p->primitives.limiter.rho[0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.rho[0] *= alpha;
    p->primitives.gradients.rho[1] *= alpha;
    p->primitives.gradients.rho[2] *= alpha;
  }

  gradtrue = sqrtf(gradv[0][0] * gradv[0][0] + gradv[0][1] * gradv[0][1] +
                   gradv[0][2] * gradv[0][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.v[0][1] - p->primitives.v[0];
    const float gradmin = p->primitives.v[0] - p->primitives.limiter.v[0][0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.v[0][0] *= alpha;
    p->primitives.gradients.v[0][1] *= alpha;
    p->primitives.gradients.v[0][2] *= alpha;
  }

  gradtrue = sqrtf(gradv[1][0] * gradv[1][0] + gradv[1][1] * gradv[1][1] +
                   gradv[1][2] * gradv[1][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.v[1][1] - p->primitives.v[1];
    const float gradmin = p->primitives.v[1] - p->primitives.limiter.v[1][0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.v[1][0] *= alpha;
    p->primitives.gradients.v[1][1] *= alpha;
    p->primitives.gradients.v[1][2] *= alpha;
  }

  gradtrue = sqrtf(gradv[2][0] * gradv[2][0] + gradv[2][1] * gradv[2][1] +
                   gradv[2][2] * gradv[2][2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.v[2][1] - p->primitives.v[2];
    const float gradmin = p->primitives.v[2] - p->primitives.limiter.v[2][0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.v[2][0] *= alpha;
    p->primitives.gradients.v[2][1] *= alpha;
    p->primitives.gradients.v[2][2] *= alpha;
  }

  gradtrue =
      sqrtf(gradP[0] * gradP[0] + gradP[1] * gradP[1] + gradP[2] * gradP[2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.P[1] - p->primitives.P;
    const float gradmin = p->primitives.P - p->primitives.limiter.P[0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.P[0] *= alpha;
    p->primitives.gradients.P[1] *= alpha;
    p->primitives.gradients.P[2] *= alpha;
  }

  gradtrue =
      sqrtf(gradA[0] * gradA[0] + gradA[1] * gradA[1] + gradA[2] * gradA[2]);
  if (gradtrue) {
    gradtrue *= p->primitives.limiter.maxr;
    const float gradmax = p->primitives.limiter.A[1] - p->primitives.A;
    const float gradmin = p->primitives.A - p->primitives.limiter.A[0];
    const float gradtrue_inv = 1.f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    p->primitives.gradients.A[0] *= alpha;
    p->primitives.gradients.A[1] *= alpha;
    p->primitives.gradients.A[2] *= alpha;
  }
}

#endif /* SWIFT_GIZMO_MFV_SLOPE_LIMITER_CELL_H */
