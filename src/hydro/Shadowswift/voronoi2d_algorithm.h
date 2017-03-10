/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com).
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

#ifndef SWIFT_VORONOIXD_ALGORITHM_H
#define SWIFT_VORONOIXD_ALGORITHM_H

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "voronoi2d_cell.h"

#define VORONOI2D_BOX_LEFT 18446744073709551602llu
#define VORONOI2D_BOX_RIGHT 18446744073709551603llu
#define VORONOI2D_BOX_TOP 18446744073709551604llu
#define VORONOI2D_BOX_BOTTOM 18446744073709551605llu

extern float global_voronoi_box_anchor[2];
extern float global_voronoi_box_side[2];

#define VORONOI_DECLARE_GLOBAL_VARIABLES() \
  float global_voronoi_box_anchor[2];      \
  float global_voronoi_box_side[2];

#define VORONOI2D_BOX_ANCHOR_X global_voronoi_box_anchor[0]
#define VORONOI2D_BOX_ANCHOR_Y global_voronoi_box_anchor[1]
#define VORONOI2D_BOX_SIDE_X global_voronoi_box_side[0]
#define VORONOI2D_BOX_SIDE_Y global_voronoi_box_side[1]

__attribute__((always_inline)) INLINE static void voronoi_set_box(float *anchor,
                                                                  float *side) {

  global_voronoi_box_anchor[0] = anchor[0];
  global_voronoi_box_anchor[1] = anchor[1];

  global_voronoi_box_side[0] = side[0];
  global_voronoi_box_side[1] = side[1];
}

/**
 * @brief Initialize a 2D Voronoi cell
 *
 * @param cell 2D Voronoi cell to initialize.
 * @param x Position of the generator of the cell.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, double *x) {

  cell->x[0] = x[0];
  cell->x[1] = x[1];

  cell->nvert = 4;

  cell->vertices[0][0] = VORONOI2D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[0][1] = VORONOI2D_BOX_ANCHOR_Y - cell->x[1];

  cell->vertices[1][0] = VORONOI2D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[1][1] =
      VORONOI2D_BOX_ANCHOR_Y + VORONOI2D_BOX_SIDE_Y - cell->x[1];

  cell->vertices[2][0] =
      VORONOI2D_BOX_ANCHOR_X + VORONOI2D_BOX_SIDE_X - cell->x[0];
  cell->vertices[2][1] =
      VORONOI2D_BOX_ANCHOR_Y + VORONOI2D_BOX_SIDE_Y - cell->x[1];

  cell->vertices[3][0] =
      VORONOI2D_BOX_ANCHOR_X + VORONOI2D_BOX_SIDE_X - cell->x[0];
  cell->vertices[3][1] = VORONOI2D_BOX_ANCHOR_Y - cell->x[1];

  cell->ngbs[0] = VORONOI2D_BOX_LEFT;
  cell->ngbs[1] = VORONOI2D_BOX_TOP;
  cell->ngbs[2] = VORONOI2D_BOX_RIGHT;
  cell->ngbs[3] = VORONOI2D_BOX_BOTTOM;

  cell->volume = 0.0f;
  cell->centroid[0] = 0.0f;
  cell->centroid[1] = 0.0f;
}

/**
 * @brief Interact a 2D Voronoi cell with a particle with given relative
 * position and ID
 *
 * @param cell 2D Voronoi cell.
 * @param dx Relative position of the interacting generator w.r.t. the cell
 * generator (in fact: dx = generator - neighbour).
 * @param id ID of the interacting neighbour.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, float *dx, unsigned long long id) {}

/**
 * @brief Finalize a 2D Voronoi cell
 *
 * @param cell 2D Voronoi cell.
 * @return Maximal radius that could still change the structure of the cell.
 */
__attribute__((always_inline)) INLINE float voronoi_cell_finalize(
    struct voronoi_cell *cell) {

  int i;
  float vertices[VORONOI2D_MAXNUMVERT][2];
  float A, x[2], y[2], r2, r2max;

  /* make a copy of the vertices (they are overwritten when the face midpoints
     are computed */
  for (i = 0; i < cell->nvert; ++i) {
    vertices[i][0] = cell->vertices[i][0];
    vertices[i][1] = cell->vertices[i][1];
  }

  r2max = 0.0f;
  for (i = 0; i < cell->nvert; ++i) {
    if (i < cell->nvert - 1) {
      x[0] = vertices[i][0];
      y[0] = vertices[i][1];
      x[1] = vertices[i + 1][0];
      y[1] = vertices[i + 1][1];
    } else {
      x[0] = vertices[i][0];
      y[0] = vertices[i][1];
      x[1] = vertices[0][0];
      y[1] = vertices[0][1];
    }
    A = x[1] * y[0] - x[0] * y[1];
    cell->volume += A;
    cell->centroid[0] += (x[0] + x[1]) * A;
    cell->centroid[1] += (y[0] + y[1]) * A;

    cell->face_midpoints[i][0] = 0.5f * (x[0] + x[1]) + cell->x[0];
    cell->face_midpoints[i][1] = 0.5f * (y[0] + y[1]) + cell->x[1];

    r2 = x[0] * x[0] + y[0] * y[0];
    r2max = max(r2max, r2);
  }

  cell->volume *= 0.5f;
  A = 6 * cell->volume;
  cell->centroid[0] /= A;
  cell->centroid[1] /= A;

  cell->centroid[0] += cell->x[0];
  cell->centroid[1] += cell->x[1];

  return 2.0f * sqrtf(r2max);
}

/**
 * @brief Get the oriented surface area and midpoint of the face between a
 * 2D Voronoi cell and the given neighbour
 *
 * @param cell 2D Voronoi cell.
 * @param ngb ID of a particle that is possibly a neighbour of this cell.
 * @param midpoint Array to store the relative position of the face in.
 * @return 0 if the given neighbour is not a neighbour, surface area otherwise.
 */
__attribute__((always_inline)) INLINE float voronoi_get_face(
    struct voronoi_cell *cell, unsigned long long ngb, float *midpoint) {

  /* look up the neighbour */
  int i = 0;
  while (i < cell->nvert && cell->ngbs[i] != ngb) {
    ++i;
  }

  if (i == cell->nvert) {
    /* The given cell is not a neighbour. */
    return 0.0f;
  }

  midpoint[0] = cell->face_midpoints[i][0];
  midpoint[1] = cell->face_midpoints[i][1];

  return cell->face_lengths[i];
}

/**
 * @brief Get the centroid of a 2D Voronoi cell
 *
 * @param cell 2D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    struct voronoi_cell *cell, float *centroid) {

  centroid[0] = cell->centroid[0];
  centroid[1] = cell->centroid[1];
}

#endif  // SWIFT_VORONOIXD_ALGORITHM_H
