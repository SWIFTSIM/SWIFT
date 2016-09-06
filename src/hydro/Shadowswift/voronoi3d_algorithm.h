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

#ifndef SWIFT_VORONOI3D_ALGORITHM_H
#define SWIFT_VORONOI3D_ALGORITHM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "inline.h"
#include "voronoi3d_cell.h"

/* Tolerance parameter used to decide when to use more precise geometric
   criteria */
#define VORONOI3D_TOLERANCE 1.e-6f

/* Box boundary flags used to signal cells neighbouring the box boundary
   These values correspond to the top range of possible 64-bit integers, and
   we make the strong assumption that there will never be a particle that has
   one of these as particle ID. */
#define VORONOI3D_BOX_FRONT 18446744073709551600llu
#define VORONOI3D_BOX_BACK 18446744073709551601llu
#define VORONOI3D_BOX_TOP 18446744073709551602llu
#define VORONOI3D_BOX_BOTTOM 18446744073709551603llu
#define VORONOI3D_BOX_LEFT 18446744073709551604llu
#define VORONOI3D_BOX_RIGHT 18446744073709551605llu

/* Bottom front left corner and side lengths of the large box that contains all
   particles and is used as initial cell at the start of the construction */
/* We should make sure that this box is either so large a particle can never
   fall outside (by using FLT_MAX if that works), or is initialized to be larger
   than the (periodic) simulation box */
#define VORONOI3D_BOX_ANCHOR_X -2.0f
#define VORONOI3D_BOX_ANCHOR_Y -2.0f
#define VORONOI3D_BOX_ANCHOR_Z -2.0f
#define VORONOI3D_BOX_SIDE_X 6.0f
#define VORONOI3D_BOX_SIDE_Y 6.0f
#define VORONOI3D_BOX_SIDE_Z 6.0f

__attribute__((always_inline)) INLINE static float voronoi_get_box_volume() {
  return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Y * VORONOI3D_BOX_SIDE_Z;
}

__attribute__((always_inline)) INLINE static void voronoi_get_box_centroid(
    float *box_centroid) {
  box_centroid[0] = 0.5f * VORONOI3D_BOX_SIDE_X + VORONOI3D_BOX_ANCHOR_X;
  box_centroid[1] = 0.5f * VORONOI3D_BOX_SIDE_Y + VORONOI3D_BOX_ANCHOR_Y;
  box_centroid[2] = 0.5f * VORONOI3D_BOX_SIDE_Z + VORONOI3D_BOX_ANCHOR_Z;
}

__attribute__((always_inline)) INLINE static float voronoi_get_box_face(
    unsigned long long id, float *face_midpoint) {

  if (id == VORONOI3D_BOX_FRONT) {
    face_midpoint[0] = 0.5f * VORONOI3D_BOX_SIDE_X + VORONOI3D_BOX_ANCHOR_X;
    face_midpoint[1] = VORONOI3D_BOX_ANCHOR_Y;
    face_midpoint[2] = 0.5f * VORONOI3D_BOX_SIDE_Z + VORONOI3D_BOX_ANCHOR_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Z;
  }
  if (id == VORONOI3D_BOX_BACK) {
    face_midpoint[0] = 0.5f * VORONOI3D_BOX_SIDE_X + VORONOI3D_BOX_ANCHOR_X;
    face_midpoint[1] = VORONOI3D_BOX_ANCHOR_Y + VORONOI3D_BOX_SIDE_Y;
    face_midpoint[2] = 0.5f * VORONOI3D_BOX_SIDE_Z + VORONOI3D_BOX_ANCHOR_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Z;
  }

  if (id == VORONOI3D_BOX_BOTTOM) {
    face_midpoint[0] = 0.5f * VORONOI3D_BOX_SIDE_X + VORONOI3D_BOX_ANCHOR_X;
    face_midpoint[1] = 0.5f * VORONOI3D_BOX_SIDE_Y + VORONOI3D_BOX_ANCHOR_Y;
    face_midpoint[2] = VORONOI3D_BOX_ANCHOR_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Y;
  }
  if (id == VORONOI3D_BOX_TOP) {
    face_midpoint[0] = 0.5f * VORONOI3D_BOX_SIDE_X + VORONOI3D_BOX_ANCHOR_X;
    face_midpoint[1] = 0.5f * VORONOI3D_BOX_SIDE_Y + VORONOI3D_BOX_ANCHOR_Y;
    face_midpoint[2] = VORONOI3D_BOX_ANCHOR_Z + VORONOI3D_BOX_SIDE_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Y;
  }

  if (id == VORONOI3D_BOX_LEFT) {
    face_midpoint[0] = VORONOI3D_BOX_ANCHOR_X;
    face_midpoint[1] = 0.5f * VORONOI3D_BOX_SIDE_Y + VORONOI3D_BOX_ANCHOR_Y;
    face_midpoint[2] = 0.5f * VORONOI3D_BOX_SIDE_Z + VORONOI3D_BOX_ANCHOR_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Y;
  }
  if (id == VORONOI3D_BOX_RIGHT) {
    face_midpoint[0] = VORONOI3D_BOX_ANCHOR_X + VORONOI3D_BOX_SIDE_X;
    face_midpoint[1] = 0.5f * VORONOI3D_BOX_SIDE_Y + VORONOI3D_BOX_ANCHOR_Y;
    face_midpoint[2] = 0.5f * VORONOI3D_BOX_SIDE_Z + VORONOI3D_BOX_ANCHOR_Z;
    return VORONOI3D_BOX_SIDE_X * VORONOI3D_BOX_SIDE_Y;
  }

  return 0.0f;
}

/*******************************************************************************
 * 3D specific methods
 *
 * Most of these methods are based on the source code of voro++:
 *  http://math.lbl.gov/voro++/
 ******************************************************************************/

__attribute__((always_inline)) INLINE void voronoi_print_gnuplot_c(
    struct voronoi_cell *c) {

  int i, j, v;
  double *x = c->x;

  fprintf(stderr, "%g\t%g\t%g\n\n", x[0], x[1], x[2]);

  for (i = 0; i < c->nvert; i++) {
    for (j = 0; j < c->orders[i]; j++) {
      v = c->edges[c->offsets[i] + j];
      fprintf(stderr, "%g\t%g\t%g\n", c->vertices[3 * i + 0] + x[0],
              c->vertices[3 * i + 1] + x[1], c->vertices[3 * i + 2] + x[2]);
      fprintf(stderr, "%g\t%g\t%g\n\n", c->vertices[3 * v + 0] + x[0],
              c->vertices[3 * v + 1] + x[1], c->vertices[3 * v + 2] + x[2]);
    }
  }
  fprintf(stderr, "\n");
}

/**
 * @brief Print the contents of a 3D Voronoi cell
 *
 * @param cell 3D Voronoi cell
 */
__attribute__((always_inline)) INLINE void voronoi_print_cell(
    struct voronoi_cell *cell) {

  int i, j;

  for (i = 0; i < cell->nvert; i++) {
    message("%i: %g %g %g (%i)", i, cell->vertices[3 * i],
            cell->vertices[3 * i + 1], cell->vertices[3 * i + 2],
            cell->orders[i]);
    for (j = 0; j < cell->orders[i]; j++) {
      message("%i (%i)", cell->edges[cell->offsets[i] + j],
              cell->edgeindices[cell->offsets[i] + j]);
    }
  }
  message("\n");
}

/**
 * @brief Get the index of the vertex pointed to by the given edge of the given
 *  vertex
 */
__attribute__((always_inline)) INLINE int voronoi_get_edge(
    struct voronoi_cell *c, int vertex, int edge) {
  return c->edges[c->offsets[vertex] + edge];
}

/**
 * @brief Get the index of vertex in the edge list of the vertex pointed to by
 *  the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE int voronoi_get_edgeindex(
    struct voronoi_cell *c, int vertex, int edge) {
  return c->edgeindices[c->offsets[vertex] + edge];
}

/**
 * @brief Set the index of the vertex pointed to by the given edge of the given
 *  vertex
 */
__attribute__((always_inline)) INLINE void voronoi_set_edge(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->edges[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Set the index of vertex in the edge list of the vertex pointed to by
 *  the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE void voronoi_set_edgeindex(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->edgeindices[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Get the neighbour for the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE int voronoi_get_ngb(
    struct voronoi_cell *c, int vertex, int edge) {
  return c->ngbs[c->offsets[vertex] + edge];
}

/**
 * @brief Set the neighbour for the given edge of the given vertex
 */
__attribute__((always_inline)) INLINE void voronoi_set_ngb(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->ngbs[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Check if the given vertex is above, below or on the cutting plane
 */
__attribute__((always_inline)) INLINE int voronoi_test_vertex(
    float *v, float *dx, float r2, float *test, float *teststack,
    int *teststack_size) {

  *test = v[0] * dx[0] + v[1] * dx[1] + v[2] * dx[2] - r2;

  teststack[*teststack_size] = *test;
  *teststack_size = *teststack_size + 1;
  if (*teststack_size == 2 * VORONOI3D_MAXNUMVERT) {
    *teststack_size = 0;
  }

  if (*test < -VORONOI3D_TOLERANCE) {
    return -1;
  }
  if (*test > VORONOI3D_TOLERANCE) {
    return 1;
  }
  return 0;
}

/**
 * @brief Initialize the cell as a cube with side 2*h, centered on the particle
 *  position
 */
__attribute__((always_inline)) INLINE void voronoi_initialize(
    struct voronoi_cell *cell) {

  cell->nvert = 8;

  /* (0, 0, 0) -- 0 */
  cell->vertices[0] = VORONOI3D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[1] = VORONOI3D_BOX_ANCHOR_Y - cell->x[1];
  cell->vertices[2] = VORONOI3D_BOX_ANCHOR_Z - cell->x[2];

  /* (0, 0, 1)-- 1 */
  cell->vertices[3] = VORONOI3D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[4] = VORONOI3D_BOX_ANCHOR_Y - cell->x[1];
  cell->vertices[5] =
      VORONOI3D_BOX_ANCHOR_Z + VORONOI3D_BOX_SIDE_Z - cell->x[2];

  /* (0, 1, 0) -- 2 */
  cell->vertices[6] = VORONOI3D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[7] =
      VORONOI3D_BOX_ANCHOR_Y + VORONOI3D_BOX_SIDE_Y - cell->x[1];
  cell->vertices[8] = VORONOI3D_BOX_ANCHOR_Z - cell->x[2];

  /* (0, 1, 1) -- 3 */
  cell->vertices[9] = VORONOI3D_BOX_ANCHOR_X - cell->x[0];
  cell->vertices[10] =
      VORONOI3D_BOX_ANCHOR_Y + VORONOI3D_BOX_SIDE_Y - cell->x[1];
  cell->vertices[11] =
      VORONOI3D_BOX_ANCHOR_Z + VORONOI3D_BOX_SIDE_Z - cell->x[2];

  /* (1, 0, 0) -- 4 */
  cell->vertices[12] =
      VORONOI3D_BOX_ANCHOR_X + VORONOI3D_BOX_SIDE_X - cell->x[0];
  cell->vertices[13] = VORONOI3D_BOX_ANCHOR_Y - cell->x[1];
  cell->vertices[14] = VORONOI3D_BOX_ANCHOR_Z - cell->x[2];

  /* (1, 0, 1) -- 5 */
  cell->vertices[15] =
      VORONOI3D_BOX_ANCHOR_X + VORONOI3D_BOX_SIDE_X - cell->x[0];
  cell->vertices[16] = VORONOI3D_BOX_ANCHOR_Y - cell->x[1];
  cell->vertices[17] =
      VORONOI3D_BOX_ANCHOR_Z + VORONOI3D_BOX_SIDE_Z - cell->x[2];

  /* (1, 1, 0) -- 6 */
  cell->vertices[18] =
      VORONOI3D_BOX_ANCHOR_X + VORONOI3D_BOX_SIDE_X - cell->x[0];
  cell->vertices[19] =
      VORONOI3D_BOX_ANCHOR_Y + VORONOI3D_BOX_SIDE_Y - cell->x[1];
  cell->vertices[20] = VORONOI3D_BOX_ANCHOR_Z - cell->x[2];

  /* (1, 1, 1) -- 7 */
  cell->vertices[21] =
      VORONOI3D_BOX_ANCHOR_X + VORONOI3D_BOX_SIDE_X - cell->x[0];
  cell->vertices[22] =
      VORONOI3D_BOX_ANCHOR_Y + VORONOI3D_BOX_SIDE_Y - cell->x[1];
  cell->vertices[23] =
      VORONOI3D_BOX_ANCHOR_Z + VORONOI3D_BOX_SIDE_Z - cell->x[2];

  cell->orders[0] = 3;
  cell->orders[1] = 3;
  cell->orders[2] = 3;
  cell->orders[3] = 3;
  cell->orders[4] = 3;
  cell->orders[5] = 3;
  cell->orders[6] = 3;
  cell->orders[7] = 3;

  /* edges are ordered counterclockwise w.r.t. a vector pointing from the
     cell generator to the vertex
     (0, 0, 0) corner */
  cell->offsets[0] = 0;
  cell->edges[0] = 1;
  cell->edges[1] = 2;
  cell->edges[2] = 4;
  cell->edgeindices[0] = 0;
  cell->edgeindices[1] = 2;
  cell->edgeindices[2] = 0;

  /* (0, 0, 1) corner */
  cell->offsets[1] = 3;
  cell->edges[3] = 0;
  cell->edges[4] = 5;
  cell->edges[5] = 3;
  cell->edgeindices[3] = 0;
  cell->edgeindices[4] = 2;
  cell->edgeindices[5] = 1;

  /* (0, 1, 0) corner */
  cell->offsets[2] = 6;
  cell->edges[6] = 3;
  cell->edges[7] = 6;
  cell->edges[8] = 0;
  cell->edgeindices[6] = 0;
  cell->edgeindices[7] = 0;
  cell->edgeindices[8] = 1;

  /* (0, 1, 1) corner */
  cell->offsets[3] = 9;
  cell->edges[9] = 2;
  cell->edges[10] = 1;
  cell->edges[11] = 7;
  cell->edgeindices[9] = 0;
  cell->edgeindices[10] = 2;
  cell->edgeindices[11] = 0;

  /* (1, 0, 0) corner */
  cell->offsets[4] = 12;
  cell->edges[12] = 0;
  cell->edges[13] = 6;
  cell->edges[14] = 5;
  cell->edgeindices[12] = 2;
  cell->edgeindices[13] = 2;
  cell->edgeindices[14] = 0;

  /* (1, 0, 1) corner */
  cell->offsets[5] = 15;
  cell->edges[15] = 4;
  cell->edges[16] = 7;
  cell->edges[17] = 1;
  cell->edgeindices[15] = 2;
  cell->edgeindices[16] = 1;
  cell->edgeindices[17] = 1;

  /* (1, 1, 0) corner */
  cell->offsets[6] = 18;
  cell->edges[18] = 2;
  cell->edges[19] = 7;
  cell->edges[20] = 4;
  cell->edgeindices[18] = 1;
  cell->edgeindices[19] = 2;
  cell->edgeindices[20] = 1;

  /* (1, 1, 1) corner */
  cell->offsets[7] = 21;
  cell->edges[21] = 3;
  cell->edges[22] = 5;
  cell->edges[23] = 6;
  cell->edgeindices[21] = 2;
  cell->edgeindices[22] = 1;
  cell->edgeindices[23] = 1;

  /* ngbs[3*i+j] is the neighbour corresponding to the plane clockwise of
     edge j of vertex i (when going from edge j to vertex i)
     we set them to a ridiculously large value to be able to track faces without
     neighbour */
  cell->ngbs[0] = VORONOI3D_BOX_FRONT;  /* (000) - (001) */
  cell->ngbs[1] = VORONOI3D_BOX_LEFT;   /* (000) - (010) */
  cell->ngbs[2] = VORONOI3D_BOX_BOTTOM; /* (000) - (100) */

  cell->ngbs[3] = VORONOI3D_BOX_LEFT;  /* (001) - (000) */
  cell->ngbs[4] = VORONOI3D_BOX_FRONT; /* (001) - (101) */
  cell->ngbs[5] = VORONOI3D_BOX_TOP;   /* (001) - (011) */

  cell->ngbs[6] = VORONOI3D_BOX_LEFT;   /* (010) - (011) */
  cell->ngbs[7] = VORONOI3D_BOX_BACK;   /* (010) - (110) */
  cell->ngbs[8] = VORONOI3D_BOX_BOTTOM; /* (010) - (000) */

  cell->ngbs[9] = VORONOI3D_BOX_BACK;  /* (011) - (010) */
  cell->ngbs[10] = VORONOI3D_BOX_LEFT; /* (011) - (001) */
  cell->ngbs[11] = VORONOI3D_BOX_TOP;  /* (011) - (111) */

  cell->ngbs[12] = VORONOI3D_BOX_FRONT;  /* (100) - (000) */
  cell->ngbs[13] = VORONOI3D_BOX_BOTTOM; /* (100) - (110) */
  cell->ngbs[14] = VORONOI3D_BOX_RIGHT;  /* (100) - (101) */

  cell->ngbs[15] = VORONOI3D_BOX_FRONT; /* (101) - (100) */
  cell->ngbs[16] = VORONOI3D_BOX_RIGHT; /* (101) - (111) */
  cell->ngbs[17] = VORONOI3D_BOX_TOP;   /* (101) - (001) */

  cell->ngbs[18] = VORONOI3D_BOX_BOTTOM; /* (110) - (010) */
  cell->ngbs[19] = VORONOI3D_BOX_BACK;   /* (110) - (111) */
  cell->ngbs[20] = VORONOI3D_BOX_RIGHT;  /* (110) - (100) */

  cell->ngbs[21] = VORONOI3D_BOX_BACK;  /* (111) - (011) */
  cell->ngbs[22] = VORONOI3D_BOX_TOP;   /* (111) - (101) */
  cell->ngbs[23] = VORONOI3D_BOX_RIGHT; /* (111) - (110) */
}

/**
 * @brief Find an edge of the voronoi_cell that intersects the cutting plane
 *
 * @param c voronoi_cell
 * @param dx vector pointing from the midpoint of the line segment between pi
 * and pj to pj
 * @param r2 Squared norm of dx
 * @param u Distance between the plane and the closest vertex above the plane
 * @param up Index of the closest vertex above the plane
 * @param us Index of the edge of vertex up that intersects the plane
 * @param uw Result of the last test_vertex call for vertex up
 * @param l Distance between the plane and the closest vertex below the plane
 * @param lp Index of the closest vertex below the plane
 * @param ls Index of the edge of vertex lp that intersects the plane
 * @param lw Result of the last test_vertex call for vertex lp
 * @param q Distance between the plane and a testing vertex
 * @param qp Index of the testing vertex
 * @param qs Index of the edge of the testing vertex that is connected to up
 * @param qw Result of the last test_vertex call involving qp
 * @return A negative value if an error occurred, 0 if the plane does not
 * intersect the cell, 1 if nothing special happened and 2 if we have a
 * complicated setup
 */
__attribute__((always_inline)) INLINE int voronoi_intersect_find_closest_vertex(
    struct voronoi_cell *c, float *dx, float r2, float *u, int *up, int *us,
    int *uw, float *l, int *lp, int *ls, int *lw, float *q, int *qp, int *qs,
    int *qw) {

  // stack to store all vertices that have already been tested (debugging only)
  float teststack[2 * VORONOI3D_MAXNUMVERT];
  // size of the used part of the stack
  int teststack_size = 0;
  int complicated;

  // test the first vertex: uw = -1 if it is below the plane, 1 if it is above
  // 0 if it is very close to the plane, and things become complicated...
  *uw = voronoi_test_vertex(&c->vertices[0], dx, r2, u, teststack,
                            &teststack_size);
  *up = 0;
  complicated = 0;
  if ((*uw) == 0) {

    /* PATH 0 */
    complicated = 1;

  } else {

    // two options: either the vertex is above or below the plane

    if ((*uw) == 1) {

      /* PATH 1 */

      // above: try to find a vertex below
      // we test all edges of the current vertex stored in up (vertex 0) until
      // we either find one below the plane or closer to the plane
      *lp = voronoi_get_edge(c, (*up), 0);
      *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l, teststack,
                                &teststack_size);
      *us = 1;
      /* Not in while: PATH 1.0 */
      /* somewhere in while: PATH 1.1 */
      /* last valid option of while: PATH 1.2 */
      while ((*us) < c->orders[(*up)] && (*l) >= (*u)) {
        *lp = voronoi_get_edge(c, (*up), (*us));
        *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l, teststack,
                                  &teststack_size);
        (*us)++;
      }
      // we increased us too much, correct this
      (*us)--;
      if ((*l) >= (*u)) {
        /* PATH 1.3 */
        // up is the closest vertex to the plane, but is above the plane
        // since the entire cell is convex, up is the closest vertex of all
        // vertices of the cell
        // this means the entire cell is supposedly above the plane, which is
        // impossible
        message(
            "Cell completely gone! This should not happen. (l >= u, l = %g, u "
            "= %g)",
            (*l), (*u));
        return -1;
      }
      // we know that lp is closer to the plane or below the plane
      // now find the index of the edge up-lp in the edge list of lp
      *ls = voronoi_get_edgeindex(c, (*up), (*us));

      // if lp is also above the plane, replace up by lp and repeat the process
      // until lp is below the plane
      while ((*lw) == 1) {
        /* PATH 1.4 */
        *u = (*l);
        *up = (*lp);
        *us = 0;
        /* no while: PATH 1.4.0 */
        /* somewhere in while: PATH 1.4.1 */
        /* last valid option of while: PATH 1.4.2 */
        while ((*us) < (*ls) && (*l) >= (*u)) {
          *lp = voronoi_get_edge(c, (*up), (*us));
          *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l,
                                    teststack, &teststack_size);
          (*us)++;
        }
        if ((*l) >= (*u)) {
          (*us)++;
          /* no while: PATH 1.4.3 */
          /* somewhere in while: PATH 1.4.4 */
          /* last valid option of while: PATH 1.4.5 */
          while ((*us) < c->orders[(*up)] && (*l) >= (*u)) {
            *lp = voronoi_get_edge(c, (*up), (*us));
            *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l,
                                      teststack, &teststack_size);
            (*us)++;
          }
          if ((*l) >= (*u)) {
            /* PATH 1.4.6 */
            message(
                "Cell completely gone! This should not happen. (l >= u, l = "
                "%g, u = %g)",
                (*l), (*u));
            return -1;
          }
        }
        (*us)--;
        *ls = voronoi_get_edgeindex(c, (*up), (*us));
      }
      // if lp is too close to the plane, replace up by lp and proceed to
      // complicated setup
      if ((*lw) == 0) {
        /* PATH 1.5 */
        *up = (*lp);
        complicated = 1;
      }

    } else { /* if(uw == 1) */

      /* PATH 2 */

      // below: try to find a vertex above
      // we test all edges of the current vertex stored in up (vertex 0) until
      // we either find one above the plane or closer to the plane

      *qp = voronoi_get_edge(c, (*up), 0);
      *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q, teststack,
                                &teststack_size);
      *us = 1;
      /* not in while: PATH 2.0 */
      /* somewhere in while: PATH 2.1 */
      /* last valid option of while: PATH 2.2 */
      while ((*us) < c->orders[(*up)] && (*u) >= (*q)) {
        *qp = voronoi_get_edge(c, (*up), (*us));
        *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q, teststack,
                                  &teststack_size);
        (*us)++;
      }
      if ((*u) >= (*q)) {
        /* PATH 2.3 */
        // up is the closest vertex to the plane and is below the plane
        // since the cell is convex, up is the closest vertex of all vertices of
        // the cell
        // this means that the entire cell is below the plane
        /* cell unaltered */
        return 0;
      } else {
        // the last increase in the loop pushed us too far, correct this
        (*us)--;
      }

      // repeat the above process until qp is closer or above the plane
      while ((*qw) == -1) {
        /* PATH 2.4 */
        *qs = voronoi_get_edgeindex(c, (*up), (*us));
        *u = (*q);
        *up = (*qp);
        *us = 0;
        /* no while: PATH 2.4.0 */
        /* somewhere in while: PATH 2.4.1 */
        /* last valid option of while: 2.4.2 */
        while ((*us) < (*qs) && (*u) >= (*q)) {
          *qp = voronoi_get_edge(c, (*up), (*us));
          *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q,
                                    teststack, &teststack_size);
          (*us)++;
        }
        if ((*u) >= (*q)) {
          (*us)++;
          /* no while: PATH 2.4.3 */
          /* somewhere in while: PATH 2.4.4 */
          /* last valid option of while: PATH 2.4.5 */
          while ((*us) < c->orders[(*up)] && (*u) >= (*q)) {
            *qp = voronoi_get_edge(c, (*up), (*us));
            *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q,
                                      teststack, &teststack_size);
            (*us)++;
          }
          if ((*u) >= (*q)) {
            /* PATH 2.4.6 */
            /* cell unaltered */
            return 0;
          }
        }
        (*us)--;
      }
      if ((*qw) == 1) {
        // qp is above the plane: initialize lp to up and replace up by qp
        *lp = (*up);
        *ls = (*us);
        *l = (*u);
        *up = (*qp);
        *us = voronoi_get_edgeindex(c, (*lp), (*ls));
        *u = (*q);
      } else {
        /* PATH 2.5 */
        // too close to call: go to complicated setup
        *up = (*qp);
        complicated = 1;
      }

    } /* if(uw == 1) */

  } /* if(uw == 0) */

  if (complicated) {
    return 2;
  } else {
    return 1;
  }
}

/**
 * @brief Intersect particle pi with particle pj and adapt its Voronoi cell
 *  structure
 *
 * odx = x_i - x_j!!!
 */
__attribute__((always_inline)) INLINE void voronoi_intersect(
    float *odx, struct voronoi_cell *c, unsigned long long ngb) {

  // vector pointing from the midpoint of the line segment between pi and pj to
  // pj.
  float dx[3];
  // squared norm of dx
  float r2;
  // u: distance between the plane and the closest vertex above the plane (up)
  // l: distance between the plane and the closest vertex below the plane (low)
  // q: distance between the plane and the vertex that is currently being tested
  float u = 0.0f, l = 0.0f, q = 0.0f;
  // up: index of the closest vertex above the plane
  // us: index of the edge of vertex up that intersects the plane
  // uw: result of the last orientation test involving vertex u
  // same naming used for vertex l and vertex q
  int up = -1, us = -1, uw = -1, lp = -1, ls = -1, lw = -1, qp = -1, qs = -1,
      qw = -1;
  // auxiliary flag used to capture degeneracies
  int complicated = -1;

  // stack to store all vertices that have already been tested (debugging only)
  float teststack[2 * VORONOI3D_MAXNUMVERT];
  // size of the used part of the stack
  int teststack_size = 0;

  dx[0] = -0.5f * odx[0];
  dx[1] = -0.5f * odx[1];
  dx[2] = -0.5f * odx[2];
  r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  int result = voronoi_intersect_find_closest_vertex(
      c, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
  if (result < 0) {
    voronoi_print_gnuplot_c(c);
    error("Error while searching intersected edge!");
  }
  if (!result) {
    /* no intersection */
    return;
  }
  if (result == 2) {
    complicated = 1;
  } else {
    complicated = 0;
  }

  int vindex = -1;
  int visitflags[VORONOI3D_MAXNUMVERT];
  int dstack[2 * VORONOI3D_MAXNUMVERT];
  int dstack_size = 1;
  float r = 0.0f;
  int cs = -1, rp = -1;
  int double_edge = 0;
  int i = -1, j = -1, k = -1;

  /* initialize visitflags */
  for (i = 0; i < VORONOI3D_MAXNUMVERT; i++) {
    visitflags[i] = 0;
  }

  if (complicated) {

    // the search routine detected a vertex very close to the plane
    // the index of this vertex is stored in up
    // we proceed by checking the edges of this vertex

    lp = voronoi_get_edge(c, up, 0);
    lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                             &teststack_size);

    // the first edge can be below, above or on the plane
    if (lw != -1) {

      // above or on the plane: we try to find one below the plane

      rp = lw;
      i = 1;
      lp = voronoi_get_edge(c, up, i);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      while (lw != -1) {
        i++;
        if (i == c->orders[up]) {
          // none of the edges of up is below the plane. Since the cell is
          // supposed to be convex, this means the entire cell is above or on
          // the plane. This should not happen...
          voronoi_print_gnuplot_c(c);
          error(
              "Cell completely gone! This should not happen. (i == "
              "c->order[up], i = %d, c->orders[up] = %d, up = %d)\n"
              "dx: [%g %g %g]\nv[up]: [%g %g %g]\nx: [%g %g %g]",
              i, c->orders[up], up, dx[0], dx[1], dx[2], c->vertices[3 * up],
              c->vertices[3 * up + 1], c->vertices[3 * up + 2], c->x[0],
              c->x[1], c->x[2]);
        }
        lp = voronoi_get_edge(c, up, i);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      j = i + 1;
      while (j < c->orders[up]) {
        lp = voronoi_get_edge(c, up, j);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
        if (lw != -1) {
          break;
        }
        j++;
      }

      if (j == c->orders[up] && i == 1 && rp == 0) {
        k = c->orders[up];
        double_edge = 1;
      } else {
        k = j - i + 2;
      }

      /* create new order k vertex */
      vindex = c->nvert;
      c->nvert++;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }

      k = 1;

      visitflags[vindex] = 0;
      c->vertices[3 * vindex + 0] = c->vertices[3 * up + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * up + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * up + 2];

      us = i - 1;
      if (i < 0) {
        i = c->orders[up] - 1;
      }
      while (i < j) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        voronoi_set_edge(c, up, i, -1);
        i++;
        k++;
      }
      if (i == c->orders[up]) {
        qs = 0;
      } else {
        qs = i;
      }

    } else { /* if(lw != -1) */

      i = c->orders[up] - 1;
      lp = voronoi_get_edge(c, up, i);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      while (lw == -1) {
        i--;
        if (i == 0) {
          /* cell unaltered */
          return;
        }
        lp = voronoi_get_edge(c, up, i);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      j = 1;
      qp = voronoi_get_edge(c, up, j);
      qw = voronoi_test_vertex(&c->vertices[3 * qp], dx, r2, &q, teststack,
                               &teststack_size);
      while (qw == -1) {
        j++;
        qp = voronoi_get_edge(c, up, j);
        qw = voronoi_test_vertex(&c->vertices[3 * qp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      if (i == j && qw == 0) {
        double_edge = 1;
        k = c->orders[up];
      } else {
        k = c->orders[up] - i + j + 1;
      }

      /* create new order k vertex */
      vindex = c->nvert;
      c->nvert++;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }
      k = 1;

      visitflags[vindex] = 0;
      c->vertices[3 * vindex + 0] = c->vertices[3 * up + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * up + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * up + 2];

      us = i;
      i++;
      while (i < c->orders[up]) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        voronoi_set_edge(c, up, i, -1);
        i++;
        k++;
      }
      i = 0;
      while (i < j) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        voronoi_set_edge(c, up, i, -1);
        i++;
        k++;
      }
      qs = j;
    }

    if (!double_edge) {
      voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, qs));
      voronoi_set_ngb(c, vindex, 0, ngb);
    } else {
      voronoi_set_ngb(c, vindex, 0, voronoi_get_ngb(c, up, qs));
    }

    dstack[0] = up;

    cs = k;
    qp = up;
    q = u;
    i = voronoi_get_edge(c, up, us);
    us = voronoi_get_edgeindex(c, up, us);
    up = i;
    visitflags[qp] = -vindex;

  } else { /* if(complicated) */

    r = u / (u - l);
    l = 1.0f - r;

    /* create new order 3 vertex */
    vindex = c->nvert;
    c->nvert++;
    if (c->nvert == VORONOI3D_MAXNUMVERT) {
      error("Too many vertices!");
    }
    c->orders[vindex] = 3;
    c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
    if (c->offsets[vindex] + 3 >= VORONOI3D_MAXNUMEDGE) {
      error("Too many edges!");
    }

    visitflags[vindex] = 0;
    c->vertices[3 * vindex + 0] =
        c->vertices[3 * lp + 0] * r + c->vertices[3 * up + 0] * l;
    c->vertices[3 * vindex + 1] =
        c->vertices[3 * lp + 1] * r + c->vertices[3 * up + 1] * l;
    c->vertices[3 * vindex + 2] =
        c->vertices[3 * lp + 2] * r + c->vertices[3 * up + 2] * l;

    dstack[0] = up;

    voronoi_set_edge(c, vindex, 1, lp);
    voronoi_set_edgeindex(c, vindex, 1, ls);
    voronoi_set_edge(c, lp, ls, vindex);
    voronoi_set_edgeindex(c, lp, ls, 1);
    voronoi_set_edge(c, up, us, -1);

    voronoi_set_ngb(c, vindex, 0, ngb);
    voronoi_set_ngb(c, vindex, 1, voronoi_get_ngb(c, up, us));
    voronoi_set_ngb(c, vindex, 2, voronoi_get_ngb(c, lp, ls));

    qs = us + 1;
    if (qs == c->orders[up]) {
      qs = 0;
    }
    qp = up;
    q = u;

    cs = 2;

  } /* if(complicated) */

  int cp = -1;
  int iqs = -1;
  int new_double_edge = -1;

  cp = vindex;
  rp = vindex;
  while (qp != up || qs != us) {

    lp = voronoi_get_edge(c, qp, qs);
    lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                             &teststack_size);
    if (lw == 0) {

      k = 1;
      if (double_edge) {
        k = 0;
      }
      qs = voronoi_get_edgeindex(c, qp, qs);
      qp = lp;
      iqs = qs;

      k++;
      qs++;
      if (qs == c->orders[qp]) {
        qs = 0;
      }
      lp = voronoi_get_edge(c, qp, qs);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      while (lw == -1) {
        k++;
        qs++;
        if (qs == c->orders[qp]) {
          qs = 0;
        }
        lp = voronoi_get_edge(c, qp, qs);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      j = visitflags[qp];
      if (qp == up && qs == us) {
        new_double_edge = 0;
        if (j > 0) {
          k += c->orders[j];
        }
      } else {
        if (j > 0) {
          k += c->orders[j];
          if (lw == 0) {
            i = -visitflags[lp];
            if (i > 0) {
              if (voronoi_get_edge(c, i, c->orders[i] - 1) == j) {
                new_double_edge = 1;
                k--;
              } else {
                new_double_edge = 0;
              }
            } else {
              if (j == rp && lp == up && voronoi_get_edge(c, qp, qs) == us) {
                new_double_edge = 1;
                k--;
              } else {
                new_double_edge = 0;
              }
            }
          } else {
            new_double_edge = 0;
          }
        } else {
          if (lw == 0) {
            i = -visitflags[lp];
            if (i == cp) {
              new_double_edge = 1;
              k--;
            } else {
              new_double_edge = 0;
            }
          } else {
            new_double_edge = 0;
          }
        }
      }

      if (j > 0) {
        error("Case not handled!");
      }

      /* create new order k vertex */
      vindex = c->nvert;
      c->nvert++;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        fprintf(stderr, "dx: %g %g %g\n", dx[0], dx[1], dx[2]);
        fprintf(stderr, "r2: %g\n", r2);
        fprintf(stderr, "u: %g, l: %g, q: %g\n", u, l, q);
        fprintf(stderr, "up: %i, us: %i, uw: %i\n", up, us, uw);
        fprintf(stderr, "lp: %i, ls: %i, lw: %i\n", lp, ls, lw);
        fprintf(stderr, "qp: %i, qs: %i, qw: %i\n", qp, qs, qw);
        fprintf(stderr, "complicated: %i\n", complicated);
        fprintf(stderr, "vindex: %i\n", vindex);
        for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
          fprintf(stderr, "visitflags[%i]: %i\n", ic, visitflags[ic]);
        }
        for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
          fprintf(stderr, "dstack[%i]: %i\n", ic, dstack[ic]);
        }
        fprintf(stderr, "dstack_size: %i\n", dstack_size);
        for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
          fprintf(stderr, "teststack[%i]: %g\n", ic, teststack[ic]);
        }
        fprintf(stderr, "teststack_size: %i\n", teststack_size);
        fprintf(stderr, "r: %g\n", r);
        fprintf(stderr, "cs: %i, rp: %i\n", cs, rp);
        fprintf(stderr, "double_edge: %i\n", double_edge);
        fprintf(stderr, "i: %i, j: %i, k: %i\n", i, j, k);
        fprintf(stderr, "cp: %i, iqs: %i\n", cp, iqs);
        fprintf(stderr, "new_double_edge: %i\n", new_double_edge);
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }

      visitflags[vindex] = 0;
      c->vertices[3 * vindex + 0] = c->vertices[3 * qp + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * qp + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * qp + 2];
      visitflags[qp] = -vindex;
      dstack[dstack_size] = qp;
      dstack_size++;
      j = vindex;
      i = 0;

      if (!double_edge) {
        voronoi_set_ngb(c, j, i, ngb);
        voronoi_set_edge(c, j, i, cp);
        voronoi_set_edgeindex(c, j, i, cs);
        voronoi_set_edge(c, cp, cs, j);
        voronoi_set_edgeindex(c, cp, cs, i);
        i++;
      }

      qs = iqs;
      iqs = k - 1;
      if (new_double_edge) {
        iqs = k;
      }
      while (i < iqs) {
        qs++;
        if (qs == c->orders[qp]) {
          qs = 0;
        }
        lp = voronoi_get_edge(c, qp, qs);
        ls = voronoi_get_edgeindex(c, qp, qs);
        voronoi_set_ngb(c, j, i, voronoi_get_ngb(c, qp, qs));
        voronoi_set_edge(c, j, i, lp);
        voronoi_set_edgeindex(c, j, i, ls);
        voronoi_set_edge(c, lp, ls, j);
        voronoi_set_edgeindex(c, lp, ls, i);
        voronoi_set_edge(c, qp, qs, -1);
        i++;
      }
      qs++;
      if (qs == c->orders[qp]) {
        qs = 0;
      }
      cs = i;
      cp = j;

      if (new_double_edge) {
        voronoi_set_ngb(c, j, 0, voronoi_get_ngb(c, qp, qs));
      } else {
        voronoi_set_ngb(c, j, cs, voronoi_get_ngb(c, qp, qs));
      }

      double_edge = new_double_edge;

    } else { /* if(lw == 0) */

      if (lw == 1) {
        qs = voronoi_get_edgeindex(c, qp, qs) + 1;
        if (qs == c->orders[lp]) {
          qs = 0;
        }
        qp = lp;
        q = l;
        dstack[dstack_size] = qp;
        dstack_size++;
      } else {

        r = q / (q - l);
        l = 1.0f - r;

        /* create new order 3 vertex */
        vindex = c->nvert;
        c->nvert++;
        if (c->nvert == VORONOI3D_MAXNUMVERT) {
          fprintf(stderr, "dx: %g %g %g\n", dx[0], dx[1], dx[2]);
          fprintf(stderr, "r2: %g\n", r2);
          fprintf(stderr, "u: %g, l: %g, q: %g\n", u, l, q);
          fprintf(stderr, "up: %i, us: %i, uw: %i\n", up, us, uw);
          fprintf(stderr, "lp: %i, ls: %i, lw: %i\n", lp, ls, lw);
          fprintf(stderr, "qp: %i, qs: %i, qw: %i\n", qp, qs, qw);
          fprintf(stderr, "complicated: %i\n", complicated);
          fprintf(stderr, "vindex: %i\n", vindex);
          for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
            fprintf(stderr, "visitflags[%i]: %i\n", ic, visitflags[ic]);
          }
          for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
            fprintf(stderr, "dstack[%i]: %i\n", ic, dstack[ic]);
          }
          fprintf(stderr, "dstack_size: %i\n", dstack_size);
          for (int ic = 0; ic < VORONOI3D_MAXNUMVERT; ic++) {
            fprintf(stderr, "teststack[%i]: %g\n", ic, teststack[ic]);
          }
          fprintf(stderr, "teststack_size: %i\n", teststack_size);
          fprintf(stderr, "r: %g\n", r);
          fprintf(stderr, "cs: %i, rp: %i\n", cs, rp);
          fprintf(stderr, "double_edge: %i\n", double_edge);
          fprintf(stderr, "i: %i, j: %i, k: %i\n", i, j, k);
          fprintf(stderr, "cp: %i, iqs: %i\n", cp, iqs);
          fprintf(stderr, "new_double_edge: %i\n", new_double_edge);
          error("Too many vertices!");
        }
        c->orders[vindex] = 3;
        c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
        if (c->offsets[vindex] + 3 >= VORONOI3D_MAXNUMEDGE) {
          error("Too many edges!");
        }

        c->vertices[3 * vindex + 0] =
            c->vertices[3 * lp + 0] * r + c->vertices[3 * qp + 0] * l;
        c->vertices[3 * vindex + 1] =
            c->vertices[3 * lp + 1] * r + c->vertices[3 * qp + 1] * l;
        c->vertices[3 * vindex + 2] =
            c->vertices[3 * lp + 2] * r + c->vertices[3 * qp + 2] * l;

        ls = voronoi_get_edgeindex(c, qp, qs);
        voronoi_set_edge(c, vindex, 0, cp);
        voronoi_set_edge(c, vindex, 1, lp);
        voronoi_set_edgeindex(c, vindex, 0, cs);
        voronoi_set_edgeindex(c, vindex, 1, ls);
        voronoi_set_edge(c, lp, ls, vindex);
        voronoi_set_edgeindex(c, lp, ls, 1);
        voronoi_set_edge(c, cp, cs, vindex);
        voronoi_set_edgeindex(c, cp, cs, 0);
        voronoi_set_edge(c, qp, qs, -1);

        voronoi_set_ngb(c, vindex, 0, ngb);
        voronoi_set_ngb(c, vindex, 1, voronoi_get_ngb(c, qp, qs));
        voronoi_set_ngb(c, vindex, 2, voronoi_get_ngb(c, lp, ls));

        qs++;
        if (qs == c->orders[qp]) {
          qs = 0;
        }
        cp = vindex;
        cs = 2;
      } /* if(lw == 1) */

    } /* if(lw == 0) */

  } /* while() */

  voronoi_set_edge(c, cp, cs, rp);
  voronoi_set_edge(c, rp, 0, cp);
  voronoi_set_edgeindex(c, cp, cs, 0);
  voronoi_set_edgeindex(c, rp, 0, cs);

  for (i = 0; i < dstack_size; i++) {
    for (j = 0; j < c->orders[dstack[i]]; j++) {
      if (voronoi_get_edge(c, dstack[i], j) >= 0) {
        dstack[dstack_size] = voronoi_get_edge(c, dstack[i], j);
        dstack_size++;
        voronoi_set_edge(c, dstack[i], j, -1);
        voronoi_set_edgeindex(c, dstack[i], j, -1);
      }
    }
  }

  /* remove deleted vertices from all arrays */
  struct voronoi_cell new_cell;
  // make sure the contents of the new cell are the same as for the old cell
  memcpy(&new_cell, c, sizeof(struct voronoi_cell));
  int m, n;
  for (vindex = 0; vindex < c->nvert; vindex++) {
    j = vindex;
    /* find next edge that is not deleted */
    while (j < c->nvert && voronoi_get_edge(c, j, 0) < 0) {
      j++;
    }

    if (j == c->nvert) {
      /* ready */
      break;
    }

    /* copy vertices */
    new_cell.vertices[3 * vindex + 0] = c->vertices[3 * j + 0];
    new_cell.vertices[3 * vindex + 1] = c->vertices[3 * j + 1];
    new_cell.vertices[3 * vindex + 2] = c->vertices[3 * j + 2];

    /* copy order */
    new_cell.orders[vindex] = c->orders[j];

    /* set offset */
    if (vindex) {
      new_cell.offsets[vindex] =
          new_cell.offsets[vindex - 1] + new_cell.orders[vindex - 1];
    } else {
      new_cell.offsets[vindex] = 0;
    }

    /* copy edges, edgeindices and ngbs */
    for (k = 0; k < c->orders[j]; k++) {
      voronoi_set_edge(&new_cell, vindex, k, voronoi_get_edge(c, j, k));
      voronoi_set_edgeindex(&new_cell, vindex, k,
                            voronoi_get_edgeindex(c, j, k));
      voronoi_set_ngb(&new_cell, vindex, k, voronoi_get_ngb(c, j, k));
    }

    /* update other edges */
    for (k = 0; k < c->orders[j]; k++) {
      m = voronoi_get_edge(c, j, k);
      n = voronoi_get_edgeindex(c, j, k);
      if (m < vindex) {
        voronoi_set_edge(&new_cell, m, n, vindex);
      } else {
        voronoi_set_edge(c, m, n, vindex);
      }
    }

    /* deactivate edge */
    voronoi_set_edge(c, j, 0, -1);
  }
  new_cell.nvert = vindex;

  new_cell.x[0] = c->x[0];
  new_cell.x[1] = c->x[1];
  new_cell.x[2] = c->x[2];
  new_cell.centroid[0] = c->centroid[0];
  new_cell.centroid[1] = c->centroid[1];
  new_cell.centroid[2] = c->centroid[2];
  new_cell.volume = c->volume;
  new_cell.nface = c->nface;

  /* Update the cell values. */
  voronoi3d_cell_copy(&new_cell, c);
}

/**
 * @brief Get the volume of the tetrahedron made up by the four given vertices
 *
 * The vertices are expected to be oriented. No idea how though...
 *
 * @param v1 First vertex.
 * @param v2 Second vertex.
 * @param v3 Third vertex.
 * @param v4 Fourth vertex.
 * @return Volume of the tetrahedron.
 */
__attribute__((always_inline)) INLINE float voronoi_volume_tetrahedron(
    float *v1, float *v2, float *v3, float *v4) {

  float V;
  float r1[3], r2[3], r3[3];

  r1[0] = v2[0] - v1[0];
  r1[1] = v2[1] - v1[1];
  r1[2] = v2[2] - v1[2];
  r2[0] = v3[0] - v1[0];
  r2[1] = v3[1] - v1[1];
  r2[2] = v3[2] - v1[2];
  r3[0] = v4[0] - v1[0];
  r3[1] = v4[1] - v1[1];
  r3[2] = v4[2] - v1[2];
  V = fabs(r1[0] * r2[1] * r3[2] + r1[1] * r2[2] * r3[0] +
           r1[2] * r2[0] * r3[1] - r1[2] * r2[1] * r3[0] -
           r2[2] * r3[1] * r1[0] - r3[2] * r1[1] * r2[0]);
  V /= 6.;
  return V;
}

/**
 * @brief Get the centroid of the tetrahedron made up by the four given vertices
 *
 * This time, there is no need to orient the vertices.
 *
 * @param centroid Array to store the centroid in.
 * @param v1 First vertex.
 * @param v2 Second vertex.
 * @param v3 Third vertex.
 * @param v4 Fourth vertex.
 */
__attribute__((always_inline)) INLINE void voronoi_centroid_tetrahedron(
    float *centroid, float *v1, float *v2, float *v3, float *v4) {

  centroid[0] = 0.25f * (v1[0] + v2[0] + v3[0] + v4[0]);
  centroid[1] = 0.25f * (v1[1] + v2[1] + v3[1] + v4[1]);
  centroid[2] = 0.25f * (v1[2] + v2[2] + v3[2] + v4[2]);
}

/**
 * @brief Calculate the volume and centroid of a 3D Voronoi cell
 *
 * @param cell 3D Voronoi cell.
 */
__attribute__((always_inline)) INLINE void voronoi_calculate_cell(
    struct voronoi_cell *cell) {

  float v1[3], v2[3], v3[3], v4[3];
  int i, j, k, l, m, n;
  float tcentroid[3];
  float tvol;

  /* we need to calculate the volume of the tetrahedra formed by the first
     vertex and the triangles that make up the other faces
     since we do not store faces explicitly, this means keeping track of the
     edges that have been processed somehow
     we follow the method used in voro++ and "flip" processed edges to
     negative values
     this also means that we need to process all triangles corresponding to
     an edge at once */
  cell->volume = 0.0f;
  v1[0] = cell->vertices[0];
  v1[1] = cell->vertices[1];
  v1[2] = cell->vertices[2];
  cell->centroid[0] = 0.0f;
  cell->centroid[1] = 0.0f;
  cell->centroid[2] = 0.0f;

  /* loop over all vertices (except the first one) */
  for (i = 1; i < cell->nvert; i++) {

    v2[0] = cell->vertices[3 * i + 0];
    v2[1] = cell->vertices[3 * i + 1];
    v2[2] = cell->vertices[3 * i + 2];

    /*  loop over the edges of the vertex*/
    for (j = 0; j < cell->orders[i]; j++) {

      k = voronoi_get_edge(cell, i, j);

      if (k >= 0) {

        /* mark the edge as processed */
        voronoi_set_edge(cell, i, j, -k - 1);

        l = voronoi_get_edgeindex(cell, i, j) + 1;
        if (l == cell->orders[k]) {
          l = 0;
        }
        v3[0] = cell->vertices[3 * k + 0];
        v3[1] = cell->vertices[3 * k + 1];
        v3[2] = cell->vertices[3 * k + 2];
        m = voronoi_get_edge(cell, k, l);
        voronoi_set_edge(cell, k, l, -1 - m);

        while (m != i) {
          n = voronoi_get_edgeindex(cell, k, l) + 1;
          if (n == cell->orders[m]) {
            n = 0;
          }
          v4[0] = cell->vertices[3 * m + 0];
          v4[1] = cell->vertices[3 * m + 1];
          v4[2] = cell->vertices[3 * m + 2];
          tvol = voronoi_volume_tetrahedron(v1, v2, v3, v4);
          cell->volume += tvol;
          voronoi_centroid_tetrahedron(tcentroid, v1, v2, v3, v4);
          cell->centroid[0] += tcentroid[0] * tvol;
          cell->centroid[1] += tcentroid[1] * tvol;
          cell->centroid[2] += tcentroid[2] * tvol;
          k = m;
          l = n;
          v3[0] = v4[0];
          v3[1] = v4[1];
          v3[2] = v4[2];
          m = voronoi_get_edge(cell, k, l);
          voronoi_set_edge(cell, k, l, -1 - m);
        } /* while() */

      } /* if(k >= 0) */

    } /* for(j) */

  } /* for(i) */

  cell->centroid[0] /= cell->volume;
  cell->centroid[1] /= cell->volume;
  cell->centroid[2] /= cell->volume;

  /* centroid was calculated relative w.r.t. particle position */
  cell->centroid[0] += cell->x[0];
  cell->centroid[1] += cell->x[1];
  cell->centroid[2] += cell->x[2];

  /* Reset the edges: we still need them for the face calculation */
  for (i = 0; i < VORONOI3D_MAXNUMEDGE; ++i) {
    if (cell->edges[i] < 0) {
      cell->edges[i] = -1 - cell->edges[i];
    }
  }
}

/**
 * @brief Calculate the faces for a 3D Voronoi cell. This reorganizes the
 * internal variables of the cell, so no new neighbours can be added after
 * this method has been called!
 *
 * @param cell 3D Voronoi cell.
 */
__attribute__((always_inline)) INLINE void voronoi_calculate_faces(
    struct voronoi_cell *cell) {

  int i, j, k, l, m, n;
  float area;
  float midpoint[3];
  float u[3], v[3], w[3];
  float loc_area;
  unsigned long long newngbs[VORONOI3D_MAXNUMEDGE];

  cell->nface = 0;
  for (i = 0; i < cell->nvert; ++i) {

    for (j = 0; j < cell->orders[i]; ++j) {

      k = voronoi_get_edge(cell, i, j);

      if (k >= 0) {

        newngbs[cell->nface] = voronoi_get_ngb(cell, i, j);
        area = 0.;
        midpoint[0] = 0.;
        midpoint[1] = 0.;
        midpoint[2] = 0.;
        voronoi_set_edge(cell, i, j, -1 - k);
        l = voronoi_get_edgeindex(cell, i, j) + 1;
        if (l == cell->orders[k]) {
          l = 0;
        }
        m = voronoi_get_edge(cell, k, l);
        voronoi_set_edge(cell, k, l, -1 - m);

        while (m != i) {
          n = voronoi_get_edgeindex(cell, k, l) + 1;
          if (n == cell->orders[m]) {
            n = 0;
          }
          u[0] = cell->vertices[3 * k + 0] - cell->vertices[3 * i + 0];
          u[1] = cell->vertices[3 * k + 1] - cell->vertices[3 * i + 1];
          u[2] = cell->vertices[3 * k + 2] - cell->vertices[3 * i + 2];
          v[0] = cell->vertices[3 * m + 0] - cell->vertices[3 * i + 0];
          v[1] = cell->vertices[3 * m + 1] - cell->vertices[3 * i + 1];
          v[2] = cell->vertices[3 * m + 2] - cell->vertices[3 * i + 2];
          w[0] = u[1] * v[2] - u[2] * v[1];
          w[1] = u[2] * v[0] - u[0] * v[2];
          w[2] = u[0] * v[1] - u[1] * v[0];
          loc_area = sqrtf(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
          area += loc_area;
          midpoint[0] += loc_area * (cell->vertices[3 * k + 0] +
                                     cell->vertices[3 * i + 0] +
                                     cell->vertices[3 * m + 0]);
          midpoint[1] += loc_area * (cell->vertices[3 * k + 1] +
                                     cell->vertices[3 * i + 1] +
                                     cell->vertices[3 * m + 1]);
          midpoint[2] += loc_area * (cell->vertices[3 * k + 2] +
                                     cell->vertices[3 * i + 2] +
                                     cell->vertices[3 * m + 2]);
          k = m;
          l = n;
          m = voronoi_get_edge(cell, k, l);
          voronoi_set_edge(cell, k, l, -1 - m);
        }

        cell->face_areas[cell->nface] = 0.5f * area;
        cell->face_midpoints[cell->nface][0] = midpoint[0] / area / 3.0f;
        cell->face_midpoints[cell->nface][1] = midpoint[1] / area / 3.0f;
        cell->face_midpoints[cell->nface][2] = midpoint[2] / area / 3.0f;
        /* face midpoint was calculated relative to particle position */
        cell->face_midpoints[cell->nface][0] += cell->x[0];
        cell->face_midpoints[cell->nface][1] += cell->x[1];
        cell->face_midpoints[cell->nface][2] += cell->x[2];
        cell->nface++;

        if (cell->nface == VORONOI3D_MAXFACE) {
          error("Too many faces!");
        }

      } /* if(k >= 0) */

    } /* for(j) */

  } /* for(i) */

  /* Overwrite the old neighbour array. */
  for (i = 0; i < cell->nface; ++i) {
    cell->ngbs[i] = newngbs[i];
  }
}

/*******************************************************************************
 * voronoi_algorithm interface implementations
 *
 * If you change any function parameters below, you also have to change them in
 * the 1D and 2D algorithm!
 ******************************************************************************/

/**
 * @brief Initialize a 3D Voronoi cell
 *
 * @param cell 3D Voronoi cell to initialize.
 * @param x Position of the generator of the cell->
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, double *x) {

  cell->x[0] = x[0];
  cell->x[1] = x[1];
  cell->x[2] = x[2];

  voronoi_initialize(cell);

  cell->volume = 0.0f;
  cell->centroid[0] = 0.0f;
  cell->centroid[1] = 0.0f;
  cell->centroid[2] = 0.0f;
  cell->nface = 0;
}

/**
 * @brief Interact a #D Voronoi cell with a particle with given relative
 * position and ID
 *
 * @param cell 3D Voronoi cell.
 * @param dx Relative position of the interacting generator w.r.t. the cell
 * generator (in fact: dx = generator - neighbour).
 * @param id ID of the interacting neighbour.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, float *dx, unsigned long long id) {

  voronoi_intersect(dx, cell, id);
}

/**
 * @brief Finalize a 3D Voronoi cell
 *
 * @param cell 3D Voronoi cell.
 * @return Maximal radius that could still change the structure of the cell->
 */
__attribute__((always_inline)) INLINE float voronoi_cell_finalize(
    struct voronoi_cell *cell) {

  int i;
  float max_radius, v[3], v2;

  /* Calculate the volume and centroid of the cell. */
  voronoi_calculate_cell(cell);
  /* Calculate the faces. */
  voronoi_calculate_faces(cell);

  /* Loop over the vertices and calculate the maximum radius. */
  max_radius = 0.0f;
  for (i = 0; i < cell->nvert; ++i) {
    v[0] = cell->vertices[3 * i];
    v[1] = cell->vertices[3 * i + 1];
    v[2] = cell->vertices[3 * i + 2];
    v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    max_radius = fmaxf(max_radius, v2);
  }
  max_radius = sqrtf(max_radius);

  return 2.0f * max_radius;
}

/**
 * @brief Get the surface area and midpoint of the face between a 3D Voronoi
 * cell and the given neighbour
 *
 * @param cell 3D Voronoi cell.
 * @param ngb ID of a particle that is possibly a neighbour of this cell->
 * @param midpoint Array to store the relative position of the face in.
 * @return 0 if the given neighbour is not a neighbour, the surface area of
 * the face otherwise.
 */
__attribute__((always_inline)) INLINE float voronoi_get_face(
    struct voronoi_cell *cell, unsigned long long ngb, float *midpoint) {

  int i = 0;
  while (i < cell->nface && cell->ngbs[i] != ngb) {
    i++;
  }
  if (i == cell->nface) {
    /* Ngb not found */
    return 0.0f;
  }

  midpoint[0] = cell->face_midpoints[i][0];
  midpoint[1] = cell->face_midpoints[i][1];
  midpoint[2] = cell->face_midpoints[i][2];

  return cell->face_areas[i];
}

/**
 * @brief Get the centroid of a 3D Voronoi cell
 *
 * @param cell 3D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    struct voronoi_cell *cell, float *centroid) {

  centroid[0] = cell->centroid[0];
  centroid[1] = cell->centroid[1];
  centroid[2] = cell->centroid[2];
}

#endif  // SWIFT_VORONOI3D_ALGORITHM_H
