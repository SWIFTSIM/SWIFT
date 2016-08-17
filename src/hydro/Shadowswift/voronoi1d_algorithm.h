/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#ifndef SWIFT_VORONOI1D_ALGORITHM_H
#define SWIFT_VORONOI1D_ALGORITHM_H

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "inline.h"
#include "voronoi1d_cell.h"

__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, double *x) {
  cell->x = x[0];
  cell->xL = -DBL_MAX;
  cell->xR = DBL_MAX;
  cell->idL = 0;
  cell->idR = 0;
  cell->volume = 0.0f;
  cell->centroid = 0.0f;
}

__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, float *dx, unsigned long long id) {

  /* Check for stupidity */
  if (dx[0] == 0.0f) {
    error("Cannot interact a Voronoi cell generator with itself!");
  }

  if (-dx[0] < 0.0f) {
    /* New left neighbour? */
    if (-dx[0] > cell->xL) {
      cell->xL = -dx[0];
      cell->idL = id;
    }
  } else {
    /* New right neighbour? */
    if (-dx[0] < cell->xR) {
      cell->xR = -dx[0];
      cell->idR = id;
    }
  }
}

__attribute__((always_inline)) INLINE float voronoi_cell_finalize(
    struct voronoi_cell *cell) {

  float xL, xR;
  float max_radius;

  max_radius = fmax(-cell->xL, cell->xR);
  cell->xL = xL = 0.5f * cell->xL;
  cell->xR = xR = 0.5f * cell->xR;

  message("xL: %g, xR: %g", xL, xR);

  cell->volume = xR - xL;
  cell->centroid = cell->x + 0.5f * (xL + xR);

  return max_radius;
}

#endif  // SWIFT_VORONOI1D_ALGORITHM_H
