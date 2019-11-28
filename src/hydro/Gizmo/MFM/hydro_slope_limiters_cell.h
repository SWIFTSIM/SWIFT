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
#ifndef SWIFT_GIZMO_MFM_SLOPE_LIMITER_CELL_H
#define SWIFT_GIZMO_MFM_SLOPE_LIMITER_CELL_H

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

  pi->limiter.v[0][0] = min(pj->v[0], pi->limiter.v[0][0]);
  pi->limiter.v[0][1] = max(pj->v[0], pi->limiter.v[0][1]);
  pi->limiter.v[1][0] = min(pj->v[1], pi->limiter.v[1][0]);
  pi->limiter.v[1][1] = max(pj->v[1], pi->limiter.v[1][1]);
  pi->limiter.v[2][0] = min(pj->v[2], pi->limiter.v[2][0]);
  pi->limiter.v[2][1] = max(pj->v[2], pi->limiter.v[2][1]);

  pi->limiter.P[0] = min(pj->P, pi->limiter.P[0]);
  pi->limiter.P[1] = max(pj->P, pi->limiter.P[1]);

  pi->limiter.maxr = max(r, pi->limiter.maxr);
}

#endif /* SWIFT_GIZMO_MFM_SLOPE_LIMITER_CELL_H */
