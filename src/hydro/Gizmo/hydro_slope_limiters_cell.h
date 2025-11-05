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
#ifndef SWIFT_GIZMO_SLOPE_LIMITER_CELL_H
#define SWIFT_GIZMO_SLOPE_LIMITER_CELL_H

#include "hydro_getters.h"
#include "hydro_setters.h"

#include <float.h>

/**
 * @brief Initialize variables for the cell wide slope limiter
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell_init(
    struct part* p) {

  p->limiter.rho[0] = FLT_MAX;
  p->limiter.rho[1] = -FLT_MAX;
  p->limiter.v[0][0] = FLT_MAX;
  p->limiter.v[0][1] = -FLT_MAX;
  p->limiter.v[1][0] = FLT_MAX;
  p->limiter.v[1][1] = -FLT_MAX;
  p->limiter.v[2][0] = FLT_MAX;
  p->limiter.v[2][1] = -FLT_MAX;
  p->limiter.P[0] = FLT_MAX;
  p->limiter.P[1] = -FLT_MAX;

  p->limiter.maxr = -FLT_MAX;
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
  pi->limiter.rho[0] = min(pj->rho, pi->limiter.rho[0]);
  pi->limiter.rho[1] = max(pj->rho, pi->limiter.rho[1]);

  pi->limiter.v[0][0] = min(pj->fluid_v[0], pi->limiter.v[0][0]);
  pi->limiter.v[0][1] = max(pj->fluid_v[0], pi->limiter.v[0][1]);
  pi->limiter.v[1][0] = min(pj->fluid_v[1], pi->limiter.v[1][0]);
  pi->limiter.v[1][1] = max(pj->fluid_v[1], pi->limiter.v[1][1]);
  pi->limiter.v[2][0] = min(pj->fluid_v[2], pi->limiter.v[2][0]);
  pi->limiter.v[2][1] = max(pj->fluid_v[2], pi->limiter.v[2][1]);

  pi->limiter.P[0] = min(pj->P, pi->limiter.P[0]);
  pi->limiter.P[1] = max(pj->P, pi->limiter.P[1]);

  pi->limiter.maxr = max(r, pi->limiter.maxr);
}

/**
 * @brief Slope-limit the given quantity.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_quantity(
    float* gradient, const float maxr, const float value, const float valmin,
    const float valmax) {

  float gradtrue = sqrtf(gradient[0] * gradient[0] + gradient[1] * gradient[1] +
                         gradient[2] * gradient[2]);
  if (gradtrue != 0.0f) {
    gradtrue *= maxr;
    const float gradmax = valmax - value;
    const float gradmin = value - valmin;
    const float gradtrue_inv = 1.0f / gradtrue;
    const float alpha =
        min3(1.0f, gradmax * gradtrue_inv, gradmin * gradtrue_inv);
    gradient[0] *= alpha;
    gradient[1] *= alpha;
    gradient[2] *= alpha;
  }
}

/**
 * @brief Slope limit cell gradients
 *
 * @param p Particle.
 */
__attribute__((always_inline)) INLINE static void hydro_slope_limit_cell(
    struct part* p) {

  float W[5];
  float gradrho[3], gradvx[3], gradvy[3], gradvz[3], gradP[3];
  float rholim[2], vxlim[2], vylim[2], vzlim[2], Plim[2], maxr;

  hydro_part_get_primitive_variables(p, W);
  hydro_part_get_gradients(p, gradrho, gradvx, gradvy, gradvz, gradP);
  hydro_part_get_slope_limiter(p, rholim, vxlim, vylim, vzlim, Plim, &maxr);

  hydro_slope_limit_quantity(gradrho, maxr, W[0], rholim[0], rholim[1]);
  hydro_slope_limit_quantity(gradvx, maxr, W[1], vxlim[0], vxlim[1]);
  hydro_slope_limit_quantity(gradvy, maxr, W[2], vylim[0], vylim[1]);
  hydro_slope_limit_quantity(gradvz, maxr, W[3], vzlim[0], vzlim[1]);
  hydro_slope_limit_quantity(gradP, maxr, W[4], Plim[0], Plim[1]);

  hydro_part_set_gradients(p, gradrho, gradvx, gradvy, gradvz, gradP);
}

#endif /* SWIFT_GIZMO_SLOPE_LIMITER_CELL_H */
