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

#include "error.h"
#include "inline.h"
#include "minmax.h"
#include "voronoi2d_cell.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

/* Check if the number of vertices exceeds the maximal allowed number */
#define VORONOI_CHECK_SIZE()          \
  if (nvert > VORONOI2D_MAXNUMVERT) { \
    error("Too many vertices!");      \
  }

/* IDs used to keep track of cells neighbouring walls of the simulation box
   This will only work if these IDs are never used for actual particles (which
   in practice means you want to have less than 2^63-4 (~9e18) particles in your
   simulation) */
#define VORONOI2D_BOX_LEFT 18446744073709551602llu
#define VORONOI2D_BOX_RIGHT 18446744073709551603llu
#define VORONOI2D_BOX_TOP 18446744073709551604llu
#define VORONOI2D_BOX_BOTTOM 18446744073709551605llu

#define VORONOI2D_TOLERANCE 1.e-6f

/**
 * @brief Initialize a 2D Voronoi cell.
 *
 * @param cell 2D Voronoi cell to initialize.
 * @param x Position of the generator of the cell.
 * @param anchor Anchor of the simulation box containing all particles.
 * @param side Side lengths of the simulation box containing all particles.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, const double *x, const double *anchor,
    const double *side) {

  /* Set the position of the generator of the cell (for reference) */
  cell->x[0] = x[0];
  cell->x[1] = x[1];

  /* Initialize the cell as a box with the same extents as the simulation box
     (note: all vertex coordinates are relative w.r.t. the cell generator) */
  cell->nvert = 4;

  cell->vertices[0][0] = anchor[0] - cell->x[0];
  cell->vertices[0][1] = anchor[1] - cell->x[1];

  cell->vertices[1][0] = anchor[0] - cell->x[0];
  cell->vertices[1][1] = anchor[1] + side[1] - cell->x[1];

  cell->vertices[2][0] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[2][1] = anchor[1] + side[1] - cell->x[1];

  cell->vertices[3][0] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[3][1] = anchor[1] - cell->x[1];

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
 * position and ID.
 *
 * @param cell 2D Voronoi cell.
 * @param dx Relative position of the interacting generator w.r.t. the cell
 * generator (in fact: dx = generator - neighbour).
 * @param id ID of the interacting neighbour.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, const float *dx, unsigned long long id) {

  /* variables used for geometrical tests */
  float half_dx[2];
  float r2;
  /* variables used to store test results */
  float test, b1, b2, a1, a2;
  /* general loop index */
  int i;
  /* variables used to store indices of intersected edges */
  int index_above1, index_above2;
  int index_below1, index_below2;
  /* variable used to store directionality in edge traversal */
  int increment;
  /* new number of vertices and new vertex coordinates */
  int nvert;
  float vertices[VORONOI2D_MAXNUMVERT][2];
  unsigned long long ngbs[VORONOI2D_MAXNUMVERT];

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
  if (test < -VORONOI2D_TOLERANCE) {
/* vertex is below midline */
#ifdef VORONOI_VERBOSE
    message("First vertex is below midline (%g %g --> %g)!",
            cell->vertices[0][0] + cell->x[0],
            cell->vertices[0][1] + cell->x[1], test);
#endif

    /* store the test result; we might need it to compute the intersection
       coordinates */
    b1 = test;

    /* move on until we find a vertex that is above or on the midline */
    i = 1;
    test = cell->vertices[i][0] * half_dx[0] +
           cell->vertices[i][1] * half_dx[1] - r2;
    while (i < cell->nvert && test < 0.) {
      /* make sure we always store the latest test result */
      b1 = test;
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
#ifdef VORONOI_VERBOSE
      message("Not a neighbour!");
#endif
      return;
    }

    /* we have found an intersected edge: i-1 -> i
       we store the index of the vertices above and below the midline, make sure
       we store the test result for later intersection computation, and set the
       increment to positive, so that we look for the other intersected edge in
       clockwise direction */
    index_below1 = i - 1;
    index_above1 = i;
    a1 = test;
    increment = 1;
  } else {
/* vertex is above or on midline
   in the case where it is on the midline, we count that as above as well:
   the vertex will be removed, and a new vertex will be created at the same
   position */
#ifdef VORONOI_VERBOSE
    message("First vertex is above midline (%g %g --> %g)!",
            cell->vertices[0][0] + cell->x[0],
            cell->vertices[0][1] + cell->x[1], test);
#endif

    /* store the test result */
    a1 = test;

    /* move on until we find a vertex that is below the midline */
    i = 1;
    test = cell->vertices[i][0] * half_dx[0] +
           cell->vertices[i][1] * half_dx[1] - r2;
    while (i < cell->nvert && test > -VORONOI2D_TOLERANCE) {
      /* make sure we always store the most recent test result */
      a1 = test;
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
       we store the index of the vertices above and below the midline, make sure
       we store the test result for later intersection computation, and set the
       increment to negative, so that we look for the other intersected edge in
       counterclockwise direction */
    index_below1 = i;
    index_above1 = i - 1;
    increment = -1;
    b1 = test;
  }

#ifdef VORONOI_VERBOSE
  message("First intersected edge: %g %g --> %g %g (%i --> %i)",
          cell->vertices[index_below1][0] + cell->x[0],
          cell->vertices[index_below1][1] + cell->x[1],
          cell->vertices[index_above1][0] + cell->x[0],
          cell->vertices[index_above1][1] + cell->x[1], index_below1,
          index_above1);
#endif

  /* now we need to find the second intersected edge
     we start from the vertex above (or on) the midline and search in the
     direction opposite to the intersected edge direction until we find a vertex
     below the midline */

  /* we make sure we store the test result for the second vertex above the
     midline as well, since we need this for intersection point computations
     the second vertex can be equal to the first */
  a2 = a1;
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
  while (test > -VORONOI2D_TOLERANCE) {
    /* make sure we always store the most recent test result */
    a2 = test;
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

  index_below2 = i;
  index_above2 = i - increment;
  if (index_above2 < 0) {
    index_above2 = cell->nvert - 1;
  }
  if (index_above2 == cell->nvert) {
    index_above2 = 0;
  }
  /* we also store the test result for the second vertex below the midline */
  b2 = test;

  if (index_above1 == index_above2 && index_below1 == index_below2) {
    /* There can be only 1 vertex above or below the midline, but we need 2
       intersected edges, so if the vertices above the midline are the same, the
       ones below need to be different and vice versa */
    error("Only 1 intersected edge found!");
  }

  /* there is exactly one degenerate case we have not addressed yet: the case
     where index_above1 and index_above2 are the same and are on the midline.
     In this case we don't want to create 2 new vertices. Instead, we just keep
     index_above1, which basically means nothing happens at all and we can just
     return */
  if (index_above1 == index_above2 && a1 == 0.) {
    return;
  }

  /* to make the code below more clear, we make sure index_above1 always holds
     the first vertex to remove, and index_above2 the last one, in clockwise
     order
     This means we need to interchange 1 and 2 if we were searching in counter-
     clockwise direction above */
  if (increment < 0) {
    i = index_below1;
    index_below1 = index_below2;
    index_below2 = i;
    i = index_above1;
    index_above1 = index_above2;
    index_above2 = i;
    test = b1;
    b1 = b2;
    b2 = test;
    test = a1;
    a1 = a2;
    a2 = test;
  }

#ifdef VORONOI_VERBOSE
  message("First vertex below: %g %g (%i, %g)",
          cell->vertices[index_below1][0] + cell->x[0],
          cell->vertices[index_below1][1] + cell->x[1], index_below1, b1);
  message("First vertex above: %g %g (%i, %g)",
          cell->vertices[index_above1][0] + cell->x[0],
          cell->vertices[index_above1][1] + cell->x[1], index_above1, a1);
  message("Second vertex below: %g %g (%i, %g)",
          cell->vertices[index_below2][0] + cell->x[0],
          cell->vertices[index_below2][1] + cell->x[1], index_below2, b2);
  message("Second vertex above: %g %g (%i, %g)",
          cell->vertices[index_above2][0] + cell->x[0],
          cell->vertices[index_above2][1] + cell->x[1], index_above2, a2);
#endif

  if (b1 == 0. || b2 == 0.) {
    error("Vertex below midline is on midline!");
  }

  /* convert the test results (which correspond to the projected distance
     between the vertex and the midline) to the fractions of the intersected
     edges above and below the midline */
  test = a1 / (a1 - b1);
  a1 = test;
  b1 = 1.0f - test;

  test = a2 / (a2 - b2);
  a2 = test;
  b2 = 1.0f - test;

  /* remove the vertices above the midline, and insert two new vertices,
     corresponding to the intersection points of the intersected edges and the
     midline
     In practice, we just copy all remaining vertices, starting from the first
     vertex below the midline (in clockwise order) */
  nvert = 0;
  i = index_below2;
  while (i != index_above1) {
    vertices[nvert][0] = cell->vertices[i][0];
    vertices[nvert][1] = cell->vertices[i][1];
    ngbs[nvert] = cell->ngbs[i];
    ++nvert;
    VORONOI_CHECK_SIZE();
    ++i;
    if (i == cell->nvert) {
      i = 0;
    }
  }
  /* now add the new vertices, they are always last */
  vertices[nvert][0] = a1 * cell->vertices[index_below1][0] +
                       b1 * cell->vertices[index_above1][0];
  vertices[nvert][1] = a1 * cell->vertices[index_below1][1] +
                       b1 * cell->vertices[index_above1][1];
  ngbs[nvert] = id;
  ++nvert;
  VORONOI_CHECK_SIZE();
  vertices[nvert][0] = a2 * cell->vertices[index_below2][0] +
                       b2 * cell->vertices[index_above2][0];
  vertices[nvert][1] = a2 * cell->vertices[index_below2][1] +
                       b2 * cell->vertices[index_above2][1];
  ngbs[nvert] = cell->ngbs[index_above2];
  ++nvert;
  VORONOI_CHECK_SIZE();

  /* overwrite the original vertices */
  cell->nvert = nvert;
  for (i = 0; i < cell->nvert; ++i) {
    cell->vertices[i][0] = vertices[i][0];
    cell->vertices[i][1] = vertices[i][1];
    cell->ngbs[i] = ngbs[i];
  }
}

/**
 * @brief Finalize a 2D Voronoi cell.
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

    /* Note that we only need the RELATIVE positions of the midpoints */
    cell->face_midpoints[i][0] = 0.5f * (x[0] + x[1]);
    cell->face_midpoints[i][1] = 0.5f * (y[0] + y[1]);

    r2 = x[0] * x[0] + y[0] * y[0];
    r2max = max(r2max, r2);

    x[0] -= x[1];
    y[0] -= y[1];
    cell->face_lengths[i] = sqrtf(x[0] * x[0] + y[0] * y[0]);
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
 * 2D Voronoi cell and the given neighbour.
 *
 * @param cell 2D Voronoi cell.
 * @param ngb ID of a particle that is possibly a neighbour of this cell.
 * @param midpoint Array to store the relative position of the face in.
 * @return 0 if the given neighbour is not a neighbour, surface area otherwise.
 */
__attribute__((always_inline)) INLINE float voronoi_get_face(
    const struct voronoi_cell *cell, unsigned long long ngb, float *midpoint) {

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
  midpoint[2] = 0.0f;

  return cell->face_lengths[i];
}

/**
 * @brief Get the centroid of a 2D Voronoi cell.
 *
 * @param cell 2D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    const struct voronoi_cell *cell, float *centroid) {

  centroid[0] = cell->centroid[0];
  centroid[1] = cell->centroid[1];
  centroid[2] = 0.0f;
}

/*******************************************************************************
 ** EXTRA FUNCTIONS USED FOR DEBUGGING *****************************************
 ******************************************************************************/

/**
 * @brief Print the given cell to the stdout in a format that can be plotted
 * using gnuplot.
 *
 * @param cell voronoi_cell to print.
 */
__attribute__((always_inline)) INLINE void voronoi_print_cell(
    const struct voronoi_cell *cell) {

  int i, ip1;

  /* print cell generator */
  printf("%g %g\n\n", cell->x[0], cell->x[1]);

  /* print cell vertices */
  for (i = 0; i < cell->nvert; ++i) {
    ip1 = i + 1;
    if (ip1 == cell->nvert) {
      ip1 = 0;
    }
    printf("%g %g\n%g %g\n\n", cell->vertices[i][0] + cell->x[0],
           cell->vertices[i][1] + cell->x[1],
           cell->vertices[ip1][0] + cell->x[0],
           cell->vertices[ip1][1] + cell->x[1]);
  }
}

#endif  // SWIFT_VORONOIXD_ALGORITHM_H
