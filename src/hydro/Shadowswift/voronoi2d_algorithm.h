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

/* IDs used to keep track of cells neighbouring walls of the simulation box
   This will only work if these IDs are never used for actual particles (which
   in practice means you want to have less than 2^63-4 (~9e18) particles in your
   simulation) */
#define VORONOI2D_BOX_LEFT 18446744073709551602llu
#define VORONOI2D_BOX_RIGHT 18446744073709551603llu
#define VORONOI2D_BOX_TOP 18446744073709551604llu
#define VORONOI2D_BOX_BOTTOM 18446744073709551605llu

/* Global variables used to store the extent of the simulation box */
extern float global_voronoi_box_anchor[2];
extern float global_voronoi_box_side[2];

/* Macro used to declare the global variables in any piece of code that uses
   this header file */
#define VORONOI_DECLARE_GLOBAL_VARIABLES() \
  float global_voronoi_box_anchor[2];      \
  float global_voronoi_box_side[2];

/* Macros that make it easier to change the box extents that are used throughout
   the code */
#define VORONOI2D_BOX_ANCHOR_X global_voronoi_box_anchor[0]
#define VORONOI2D_BOX_ANCHOR_Y global_voronoi_box_anchor[1]
#define VORONOI2D_BOX_SIDE_X global_voronoi_box_side[0]
#define VORONOI2D_BOX_SIDE_Y global_voronoi_box_side[1]

/**
 * @brief Set the values of the global variables that store the simulation box
 * extents.
 *
 * @param anchor Corner of the box with the lowest coordinate values.
 * @param side Side lengths of the box.
 */
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

  /* Set the position of the generator of the cell (for reference) */
  cell->x[0] = x[0];
  cell->x[1] = x[1];

  /* Initialize the cell as a box with the same extents as the simulation box
     (note: all vertex coordinates are relative w.r.t. the cell generator) */
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

  /* The neighbours are ordered such that neighbour i shares the face in between
     vertices i and i+1 (with last vertex + 1 = first vertex)
     We walk around the cell in clockwise direction */
  cell->ngbs[0] = VORONOI2D_BOX_LEFT;
  cell->ngbs[1] = VORONOI2D_BOX_TOP;
  cell->ngbs[2] = VORONOI2D_BOX_RIGHT;
  cell->ngbs[3] = VORONOI2D_BOX_BOTTOM;

  /* These variables are initialized to zero, we will compute them after the
     neighbour iteration has finished */
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
    struct voronoi_cell *cell, float *dx, unsigned long long id) {

  float half_dx[2];
  float r2;
  float test;
  int i;
  int index_above1, index_above2;
  int index_below1, index_below2;
  int increment;

  /* The process of cutting the current cell with the midline of the generator
     and the given relative neighbour position proceeds in two steps:
      - first we need to locate an edge of the current cell that is intersected
        by this midline. Such an edge does not necessarily exist; in this case
        the given neighbour is not an actual neighbour of this cell
      - Once we have an intersected edge, we create a new edge starting at the
        intersection point. We follow the edges connected to the intersected
        edge until we find another intersected edge, and use its intersection
        point as end point of the new edge. */

  /* First, we set up some variables that are used to check if a vertex is above
     or below the midplane. */

  /* we need a vector with half the size of the vector joining generator and
     neighbour, pointing to the neighbour */
  half_dx[0] = -0.5f * dx[0];
  half_dx[1] = -0.5f * dx[1];

  /* we need the squared length of this vector */
  r2 = half_dx[0] * half_dx[0] + half_dx[1] * half_dx[1];

  /* a vertex v = (vx, vy) is above the midline if
       vx*half_dx[0] + vy*half_dx[1] > r2
     i.e., if the length of the projected vertex position is longer than the
     length of the vector pointing to the closest point on the midline (both
     vectors originate at the position of the generator)
     the vertex is below the midline if the projected position vector is shorter
     if the projected position vector has the same length, the vertex is on the
     midline */

  /* start testing a random vertex: the first one */
  test = cell->vertices[0][0] * half_dx[0] + cell->vertices[0][1] * half_dx[1] -
         r2;
  if (test < 0.) {
    /* vertex is below midline */

    /* move on until we find a vertex that is above the midline */
    i = 1;
    test = cell->vertices[i][0] * half_dx[0] +
           cell->vertices[i][1] * half_dx[1] - r2;
    while (i < cell->nvert && test < 0.) {
      ++i;
      test = cell->vertices[i][0] * half_dx[0] +
             cell->vertices[i][1] * half_dx[1] - r2;
    }

    /* loop finished, there are two possibilities:
        - i == cell->nvert, all vertices lie below the midline and the given
          neighbour is not an actual neighbour of this cell
        - test >= 0., we found a vertex above (or on) the midline */
    if (i == cell->nvert) {
      /* the given neighbour is not an actual neighbour: exit the routine */
      return;
    }

    /* we have found an intersected edge: i-1 -> i
       we store the index of the vertex above the midline, and set the increment
       for the consecutive vertex search to +1 */
    index_below1 = i - 1;
    index_above1 = i;
    increment = 1;
  } else {
    /* vertex is above midline */

    /* move on until we find a vertex that is below the midline */
    i = 1;
    test = cell->vertices[i][0] * half_dx[0] +
           cell->vertices[i][1] * half_dx[1] - r2;
    while (i < cell->nvert && test > 0.) {
      ++i;
      test = cell->vertices[i][0] * half_dx[0] +
             cell->vertices[i][1] * half_dx[1] - r2;
    }

    /* loop finished, there are two possibilities:
        - i == cell->nvert, all vertices lie above the midline. This should
          never happen.
        - test <= 0., we found a vertex below (or on) the midline */
    if (i == cell->nvert) {
      /* fatal error! */
      error("Could not find a vertex below the midline!");
    }

    /* we have found an intersected edge: i-1 -> i
       we store the index of the vertex above the midline, and set the increment
       for the consecutive vertex search to -1 */
    index_below1 = i;
    index_above1 = i - 1;
    increment = -1;
  }

  /* now we need to find the second intersected edge
     we start from the vertex above the midline and search in the direction
     opposite to the intersected edge direction until we find a vertex below the
     midline */
  i = index_above1 + increment;
  if (i < 0) {
    i = cell->nvert - 1;
  }
  if (i == cell->nvert) {
    i = 0;
  }
  test = cell->vertices[i][0] * half_dx[0] + cell->vertices[i][1] * half_dx[1] -
         r2;
  /* this loop can never deadlock, as we know there is at least 1 vertex below
     the midline */
  while (test > 0.) {
    i += increment;
    if (i < 0) {
      i = cell->nvert - 1;
    }
    if (i == cell->nvert) {
      i = 0;
    }
    test = cell->vertices[i][0] * half_dx[0] +
           cell->vertices[i][1] * half_dx[1] - r2;
  }

  if (i == index_below1) {
    /* we only found 1 intersected edge. This is impossible. */
    error("Only 1 intersected edge found!");
  }

  index_below2 = i;
  index_above2 = i - increment;
  if (index_above2 < 0) {
    index_above2 = cell->nvert - 1;
  }
  if (index_above2 == cell->nvert) {
    index_above2 = 0;
  }

  if (increment < 0) {
    /* interchange index_below and above 1 and 2 */
    i = index_below1;
    index_below1 = index_below2;
    index_below2 = i;
    i = index_above1;
    index_above1 = index_above2;
    index_above2 = i;
  }

  /* index_above1 now holds the first vertex to remove, and index_above2 the
     last, in clockwise order
     they can be the same */
}

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
