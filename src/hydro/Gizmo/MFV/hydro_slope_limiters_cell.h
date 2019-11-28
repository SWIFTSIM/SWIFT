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

  pi->primitives.limiter.maxr = max(r, pi->primitives.limiter.maxr);
}

#endif /* SWIFT_GIZMO_MFV_SLOPE_LIMITER_CELL_H */
