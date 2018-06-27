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

#ifndef SWIFT_VORONOIXD_ALGORITHM_H
#define SWIFT_VORONOIXD_ALGORITHM_H

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "inline.h"
#include "voronoi3d_cell.h"

/* For debugging purposes */
//#define LOOP_CHECK 1000

#ifdef LOOP_CHECK
/* We need to do the trickery below to get a unique counter for each call to the
   macro. This only works if the macro is never called twice on the same line.
 */
#define MERGE(a, b) a##b
#define LOOPCOUNTER_NAME(line) MERGE(loopcount, line)

/**
 * @brief Increase the given counter variable and check if it is still valid.
 *
 * @param counter Counter to increase.
 * @param line_number Line number where the while is called.
 * @return 1 if the counter is still valid, 0 otherwise.
 */
__attribute__((always_inline)) INLINE int check_counter(int *counter,
                                                        int line_number) {
  ++(*counter);
  if ((*counter) == LOOP_CHECK) {
    error("Number of iterations reached maximum (=%i) in while on line %i!",
          LOOP_CHECK, line_number);
  }
  return 1;
}

/* safewhile is a wrapper around a while that adds a unique counter variable to
   the loop that is increased by 1 for each time the loop is executed, and
   causes the code to crash if this number exceeds a given value.
   We use this to quickly enable or disable number of iterations checks for a
   large number of while loops */
#define safewhile(condition)          \
  int LOOPCOUNTER_NAME(__LINE__) = 0; \
  while (check_counter(&LOOPCOUNTER_NAME(__LINE__), __LINE__) && (condition))

#else /* LOOP_CHECK */

/* If LOOP_CHECK is not defined, safewhile and while are EXACTLY the same */
#define safewhile(condition) while (condition)

#endif /* LOOP_CHECK */

/* This flag activates a number of expensive geometrical checks that help
   finding bugs. */
//#define VORONOI3D_EXPENSIVE_CHECKS

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

/*******************************************************************************
 * 3D specific methods
 *
 * Most of these methods are based on the source code of voro++:
 *  http://math.lbl.gov/voro++/
 ******************************************************************************/

/**
 * @brief Print the given cell to the stderr in a format that can be easily
 * plotted using gnuplot.
 *
 * This method prints to the stderr instead of stdout to make it possible to use
 * it right before crashing the code.
 *
 * @param c Voronoi cell to print.
 */
__attribute__((always_inline)) INLINE void voronoi_print_gnuplot_c(
    const struct voronoi_cell *c) {

  int i, j, v;
  const double *x = c->x;

  fprintf(stderr, "%g\t%g\t%g\n\n", x[0], x[1], x[2]);

  for (i = 0; i < c->nvert; ++i) {
    for (j = 0; j < c->orders[i]; ++j) {
      v = c->edges[c->offsets[i] + j];
      if (v < 0) {
        v = -v - 1;
      }
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
    const struct voronoi_cell *cell) {

  int i, j;

  fprintf(stderr, "x: %g %g %g\n", cell->x[0], cell->x[1], cell->x[2]);
  fprintf(stderr, "nvert: %i\n", cell->nvert);

  for (i = 0; i < cell->nvert; ++i) {
    fprintf(stderr, "%i: %g %g %g (%i)\n", i, cell->vertices[3 * i],
            cell->vertices[3 * i + 1], cell->vertices[3 * i + 2],
            cell->orders[i]);
    for (j = 0; j < cell->orders[i]; ++j) {
      fprintf(stderr, "%i (%i)\n", cell->edges[cell->offsets[i] + j],
              cell->edgeindices[cell->offsets[i] + j]);
    }
  }
  fprintf(stderr, "\n");
}

/**
 * @brief Get the index of the vertex pointed to by the given edge of the given
 * vertex.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of the cell.
 * @param edge Edge of that vertex.
 * @return Index of the vertex on the other side of the edge.
 */
__attribute__((always_inline)) INLINE int voronoi_get_edge(
    const struct voronoi_cell *c, int vertex, int edge) {
  return c->edges[c->offsets[vertex] + edge];
}

/**
 * @brief Get the index of the given edge in the edge list of the vertex on the
 * other side of the edge of the given vertex.
 *
 * Suppose that the given vertex has edges [edge1, edge2, given_edge], and that
 * the vertex on the other side of given_edge has edges [edge1, given_edge,
 * edge2], then this method returns 1.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of the cell.
 * @param edge Edge of that vertex.
 * @return Index of that edge in the edge list of the vertex on the other side
 * of the edge.
 */
__attribute__((always_inline)) INLINE int voronoi_get_edgeindex(
    const struct voronoi_cell *c, int vertex, int edge) {
  return c->edgeindices[c->offsets[vertex] + edge];
}

/**
 * @brief Set the index of the vertex on the other side of the given edge of the
 * given vertex.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of the cell.
 * @param edge Edge of that vertex.
 * @param value Index of the vertex on the other side of that edge.
 */
__attribute__((always_inline)) INLINE void voronoi_set_edge(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->edges[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Set the index of the given edge in the edge list of the vertex on the
 * other side of the edge of the given vertex.
 *
 * Suppose that the given vertex has edges [edge1, edge2, given_edge], and we
 * want to tell this method that the vertex on the other side of given_edge has
 * edges [edge1, given_edge, edge2], then we need to pass on a value of 1 to
 * this method.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of that cell.
 * @param edge Edge of that vertex.
 * @param value Index of that edge in the edge list of the vertex on the other
 * side of the edge.
 */
__attribute__((always_inline)) INLINE void voronoi_set_edgeindex(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->edgeindices[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Get the neighbour for the given edge of the given vertex.
 *
 * An edge is shared by two faces, and each face has a neighbour. Luckily, each
 * edge also has two endpoints, so we can get away with storing only one
 * neighbour per endpoint of an edge. We have complete freedom in choosing which
 * neighbour to store in which endpoint, but we need to be consistent about it.
 * Here we use the following convention: if we take a vector pointing away from
 * the given vertex along the given edge direction, then we store the neighbour
 * that corresponds to the face to the right if looking to the cell from the
 * outside. This is the face that you encounter first when rotating counter-
 * clockwise around that vector, starting from outside the cell.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of that cell.
 * @param edge Edge of that vertex.
 * @return Index of the neighbour corresponding to that edge and vertex.
 */
__attribute__((always_inline)) INLINE int voronoi_get_ngb(
    const struct voronoi_cell *c, int vertex, int edge) {
  return c->ngbs[c->offsets[vertex] + edge];
}

/**
 * @brief Set the neighbour for the given edge of the given vertex.
 *
 * An edge is shared by two faces, and each face has a neighbour. Luckily, each
 * edge also has two endpoints, so we can get away with storing only one
 * neighbour per endpoint of an edge. We have complete freedom in choosing which
 * neighbour to store in which endpoint, but we need to be consistent about it.
 * Here we use the following convention: if we take a vector pointing away from
 * the given vertex along the given edge direction, then we store the neighbour
 * that corresponds to the face to the right if looking to the cell from the
 * outside. This is the face that you encounter first when rotating counter-
 * clockwise around that vector, starting from outside the cell.
 *
 * @param c 3D Voronoi cell.
 * @param vertex Index of a vertex of that cell.
 * @param edge Edge of that vertex.
 * @param value Index of the neighbour corresponding to that edge and vertex.
 */
__attribute__((always_inline)) INLINE void voronoi_set_ngb(
    struct voronoi_cell *c, int vertex, int edge, int value) {
  c->ngbs[c->offsets[vertex] + edge] = value;
}

/**
 * @brief Check if the 3D Voronoi cell is still consistent.
 *
 * A cell is consistent if its edges are consistent, i.e. if edge e of vertex v1
 * points to vertex v2, then v2 should have an edge that points to v1 as well,
 * and then the edge index of vertex v1 should contain the index of that edge
 * in the edge list of v2. We also check if all vertices have orders of at least
 * 3, and if all vertices are actually part of the vertex list.
 * Oh, and we check if the cell actually has vertices.
 *
 * @param cell 3D Voronoi cell to check
 */
__attribute__((always_inline)) INLINE void voronoi_check_cell_consistency(
    const struct voronoi_cell *c) {

  int i, j, e, l, m;

  if (c->nvert < 4) {
    error("Found cell with only %i vertices!", c->nvert);
  }

  for (i = 0; i < c->nvert; ++i) {
    if (c->orders[i] < 3) {
      voronoi_print_cell(c);
      error("Found cell with vertex of order %i!", c->orders[i]);
    }
    for (j = 0; j < c->orders[i]; ++j) {
      e = voronoi_get_edge(c, i, j);
      if (e >= c->nvert) {
        voronoi_print_cell(c);
        error("Found cell with edges that lead to non-existing vertices!");
      }
      if (e < 0) {
        continue;
      }
      l = voronoi_get_edgeindex(c, i, j);
      m = voronoi_get_edge(c, e, l);
      if (m != i) {
        /* voronoi_print_gnuplot_c(c); */
        voronoi_print_cell(c);
        fprintf(stderr, "i: %i, j: %i, e: %i, l: %i, m: %i\n", i, j, e, l, m);
        error("Cell inconsistency!");
      }
    }
  }
}

/**
 * @brief Check if the given vertex is above, below or on the cutting plane
 * defined by the given parameters.
 *
 * @param v Coordinates of a cell vertex, relative w.r.t. the position of the
 * generator of the cell.
 * @param dx Half of the relative distance vector between the position of the
 * generator of the cell and the position of the neighbouring cell that
 * generates the cutting plane, pointing from the generator position to the
 * cutting plane.
 * @param r2 Squared length of dx.
 * @param test Variable to store the result of the geometric test in, which
 * corresponds to the projected distance between the generator and the vertex
 * along dx.
 * @param teststack Stack to store the results of the N last tests in (for
 * debugging purposes only).
 * @param teststack_size Next available field in the teststack, is reset to 0 if
 * the teststack is full (so the N+1th results is overwritten; for debugging
 * purposes only).
 * @return Result of the test: -1 if the vertex is below the cutting plane, +1
 * if it is above, and 0 if it is on the cutting plane.
 */
__attribute__((always_inline)) INLINE int voronoi_test_vertex(
    const float *v, const float *dx, float r2, float *test, float *teststack,
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
 * @brief Initialize the cell as a cube that spans the entire simulation box.
 *
 * @param c 3D Voronoi cell to initialize.
 * @param anchor Anchor of the simulation box.
 * @param side Side lengths of the simulation box.
 */
__attribute__((always_inline)) INLINE void voronoi_initialize(
    struct voronoi_cell *cell, const double *anchor, const double *side) {

  cell->nvert = 8;

  /* (0, 0, 0) -- 0 */
  cell->vertices[0] = anchor[0] - cell->x[0];
  cell->vertices[1] = anchor[1] - cell->x[1];
  cell->vertices[2] = anchor[2] - cell->x[2];

  /* (0, 0, 1)-- 1 */
  cell->vertices[3] = anchor[0] - cell->x[0];
  cell->vertices[4] = anchor[1] - cell->x[1];
  cell->vertices[5] = anchor[2] + side[2] - cell->x[2];

  /* (0, 1, 0) -- 2 */
  cell->vertices[6] = anchor[0] - cell->x[0];
  cell->vertices[7] = anchor[1] + side[1] - cell->x[1];
  cell->vertices[8] = anchor[2] - cell->x[2];

  /* (0, 1, 1) -- 3 */
  cell->vertices[9] = anchor[0] - cell->x[0];
  cell->vertices[10] = anchor[1] + side[1] - cell->x[1];
  cell->vertices[11] = anchor[2] + side[2] - cell->x[2];

  /* (1, 0, 0) -- 4 */
  cell->vertices[12] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[13] = anchor[1] - cell->x[1];
  cell->vertices[14] = anchor[2] - cell->x[2];

  /* (1, 0, 1) -- 5 */
  cell->vertices[15] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[16] = anchor[1] - cell->x[1];
  cell->vertices[17] = anchor[2] + side[2] - cell->x[2];

  /* (1, 1, 0) -- 6 */
  cell->vertices[18] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[19] = anchor[1] + side[1] - cell->x[1];
  cell->vertices[20] = anchor[2] - cell->x[2];

  /* (1, 1, 1) -- 7 */
  cell->vertices[21] = anchor[0] + side[0] - cell->x[0];
  cell->vertices[22] = anchor[1] + side[1] - cell->x[1];
  cell->vertices[23] = anchor[2] + side[2] - cell->x[2];

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
 * @brief Find an edge of the voronoi_cell that intersects the cutting plane.
 *
 * There is a large number of possible paths through this method, each of which
 * is covered by a separate unit test in testVoronoi3D. Paths have been numbered
 * in the inline comments to help identify them.
 *
 * @param c 3D Voronoi cell.
 * @param dx Vector pointing from pj to the midpoint of the line segment between
 * pi and pj.
 * @param r2 Squared length of dx.
 * @param u Projected distance between the plane and the closest vertex above
 * the plane, along dx.
 * @param up Index of the closest vertex above the plane.
 * @param us Index of the edge of vertex up that intersects the plane.
 * @param uw Result of the last test_vertex call for vertex up.
 * @param l Projected distance between the plane and the closest vertex below
 * the plane, along dx.
 * @param lp Index of the closest vertex below the plane.
 * @param ls Index of the edge of vertex lp that intersects the plane.
 * @param lw Result of the last test_vertex call for vertex lp.
 * @param q Projected distance between the plane and a test vertex, along dx.
 * @param qp Index of the test vertex.
 * @param qs Index of the edge of the test vertex that is connected to up.
 * @param qw Result of the last test_vertex call involving qp.
 * @return A negative value if an error occurred, 0 if the plane does not
 * intersect the cell, 1 if nothing special happened and 2 if we have a
 * complicated setup.
 */
__attribute__((always_inline)) INLINE int voronoi_intersect_find_closest_vertex(
    struct voronoi_cell *c, const float *dx, float r2, float *u, int *up,
    int *us, int *uw, float *l, int *lp, int *ls, int *lw, float *q, int *qp,
    int *qs, int *qw) {

  /* stack to store all vertices that have already been tested (debugging
     only) */
  float teststack[2 * VORONOI3D_MAXNUMVERT];
  /* size of the used part of the stack */
  int teststack_size = 0;
  /* flag signalling a complicated setup */
  int complicated;

  /* test the first vertex: uw = -1 if it is below the plane, 1 if it is above
     0 if it is very close to the plane, and things become complicated... */
  *uw = voronoi_test_vertex(&c->vertices[0], dx, r2, u, teststack,
                            &teststack_size);
  *up = 0;
  complicated = 0;
  if ((*uw) == 0) {

    /* PATH 0 */
    complicated = 1;

  } else {

    /* two options: either the vertex is above or below the plane */

    if ((*uw) == 1) {

      /* PATH 1 */

      /* above: try to find a vertex below
         we test all edges of the current vertex stored in up (vertex 0) until
         we either find one below the plane or closer to the plane */
      *lp = voronoi_get_edge(c, (*up), 0);
      *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l, teststack,
                                &teststack_size);
      *us = 1;
      /* Not in while: PATH 1.0 */
      /* somewhere in while: PATH 1.1 */
      /* last valid option of while: PATH 1.2 */
      safewhile((*us) < c->orders[(*up)] && (*l) >= (*u)) {
        *lp = voronoi_get_edge(c, (*up), (*us));
        *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l, teststack,
                                  &teststack_size);
        ++(*us);
      }
      /* we increased us too much, correct this */
      --(*us);
      if ((*l) >= (*u)) {
        /* PATH 1.3 */
        /* up is the closest vertex to the plane, but is above the plane
           since the entire cell is convex, up is the closest vertex of all
           vertices of the cell
           this means the entire cell is supposedly above the plane, which is
           impossible */
        message(
            "Cell completely gone! This should not happen. (l >= u, l = %g, u "
            "= %g)",
            (*l), (*u));
        return -1;
      }
      /* we know that lp is closer to the plane or below the plane
         now find the index of the edge up-lp in the edge list of lp */
      *ls = voronoi_get_edgeindex(c, (*up), (*us));

      /* if lp is also above the plane, replace up by lp and repeat the process
         until lp is below the plane */
      safewhile((*lw) == 1) {
        /* PATH 1.4 */
        *u = (*l);
        *up = (*lp);
        *us = 0;
        /* no while: PATH 1.4.0 */
        /* somewhere in while: PATH 1.4.1 */
        /* last valid option of while: PATH 1.4.2 */
        safewhile((*us) < (*ls) && (*l) >= (*u)) {
          *lp = voronoi_get_edge(c, (*up), (*us));
          *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l,
                                    teststack, &teststack_size);
          ++(*us);
        }
        if ((*l) >= (*u)) {
          ++(*us);
          /* no while: PATH 1.4.3 */
          /* somewhere in while: PATH 1.4.4 */
          /* last valid option of while: PATH 1.4.5 */
          safewhile((*us) < c->orders[(*up)] && (*l) >= (*u)) {
            *lp = voronoi_get_edge(c, (*up), (*us));
            *lw = voronoi_test_vertex(&c->vertices[3 * (*lp)], dx, r2, l,
                                      teststack, &teststack_size);
            ++(*us);
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
        --(*us);
        *ls = voronoi_get_edgeindex(c, (*up), (*us));
      }
      /* if lp is too close to the plane, replace up by lp and proceed to
         complicated setup */
      if ((*lw) == 0) {
        /* PATH 1.5 */
        *up = (*lp);
        complicated = 1;
      }
    } else { /* if(uw == 1) */

      /* PATH 2 */

      /* below: try to find a vertex above
         we test all edges of the current vertex stored in up (vertex 0) until
         we either find one above the plane or closer to the plane */

      *qp = voronoi_get_edge(c, (*up), 0);
      *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q, teststack,
                                &teststack_size);
      *us = 1;
      /* not in while: PATH 2.0 */
      /* somewhere in while: PATH 2.1 */
      /* last valid option of while: PATH 2.2 */
      safewhile((*us) < c->orders[(*up)] && (*u) >= (*q)) {
        *qp = voronoi_get_edge(c, (*up), (*us));
        *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q, teststack,
                                  &teststack_size);
        ++(*us);
      }
      if ((*u) >= (*q)) {
        /* PATH 2.3 */
        /* up is the closest vertex to the plane and is below the plane
           since the cell is convex, up is the closest vertex of all vertices of
           the cell
           this means that the entire cell is below the plane
           The cell is unaltered. */
        return 0;
      } else {
        /* the last increase in the loop pushed us too far, correct this */
        --(*us);
      }

      /* repeat the above process until qp is closer or above the plane */
      safewhile((*qw) == -1) {
        /* PATH 2.4 */
        *qs = voronoi_get_edgeindex(c, (*up), (*us));
        *u = (*q);
        *up = (*qp);
        *us = 0;
        /* no while: PATH 2.4.0 */
        /* somewhere in while: PATH 2.4.1 */
        /* last valid option of while: 2.4.2 */
        safewhile((*us) < (*qs) && (*u) >= (*q)) {
          *qp = voronoi_get_edge(c, (*up), (*us));
          *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q,
                                    teststack, &teststack_size);
          ++(*us);
        }
        if ((*u) >= (*q)) {
          ++(*us);
          /* no while: PATH 2.4.3 */
          /* somewhere in while: PATH 2.4.4 */
          /* last valid option of while: PATH 2.4.5 */
          safewhile((*us) < c->orders[(*up)] && (*u) >= (*q)) {
            *qp = voronoi_get_edge(c, (*up), (*us));
            *qw = voronoi_test_vertex(&c->vertices[3 * (*qp)], dx, r2, q,
                                      teststack, &teststack_size);
            ++(*us);
          }
          if ((*u) >= (*q)) {
            /* PATH 2.4.6 */
            /* cell unaltered */
            return 0;
          }
        }
        --(*us);
      }
      if ((*qw) == 1) {
        /* qp is above the plane: initialize lp to up and replace up by qp */
        *lp = (*up);
        *ls = (*us);
        *l = (*u);
        *up = (*qp);
        *us = voronoi_get_edgeindex(c, (*lp), (*ls));
        *u = (*q);
      } else {
        /* PATH 2.5 */
        /* too close to call: go to complicated setup */
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
 * @brief Intersect the given cell with the midplane between the cell generator
 * and a neighbouring cell at the given relative position and with the given ID.
 *
 * This method is the core of the Voronoi algorithm. If anything goes wrong
 * geometrically, it most likely goes wrong somewhere within this method.
 *
 * @param c 3D Voronoi cell.
 * @param odx The original relative distance vector between the cell generator
 * and the intersecting neighbour, as it is passed on to runner_iact_density
 * (remember: odx = pi->x - pj->x).
 * @param ngb ID of the intersecting neighbour (pj->id in runner_iact_density).
 */
__attribute__((always_inline)) INLINE void voronoi_intersect(
    struct voronoi_cell *c, const float *odx, unsigned long long ngb) {

  /* vector pointing from pi to the midpoint of the line segment between pi and
     pj. This corresponds to -0.5*odx */
  float dx[3];
  /* squared norm of dx */
  float r2;
  /* u: distance between the plane and the closest vertex above the plane (up)
     l: distance between the plane and the closest vertex below the plane (lp)
     q: distance between the plane and the vertex that is currently being
     tested (qp) */
  float u = 0.0f, l = 0.0f, q = 0.0f;
  /* up: index of the closest vertex above the plane
     us: index of the edge of vertex up that intersects the plane
     uw: result of the last orientation test involving vertex u
     same naming used for vertex l and vertex q */
  int up = -1, us = -1, uw = -1, lp = -1, ls = -1, lw = -1, qp = -1, qs = -1,
      qw = -1;
  /* auxiliary flag used to capture degeneracies */
  int complicated = -1;

  /* stack to store all vertices that have already been tested (debugging
     only) */
  float teststack[2 * VORONOI3D_MAXNUMVERT];
  /* size of the used part of the stack */
  int teststack_size = 0;

#ifdef VORONOI3D_EXPENSIVE_CHECKS
  voronoi_check_cell_consistency(c);
#endif

  /* initialize dx and r2 */
  dx[0] = -0.5f * odx[0];
  dx[1] = -0.5f * odx[1];
  dx[2] = -0.5f * odx[2];
  r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* find an intersected edge of the cell */
  int result = voronoi_intersect_find_closest_vertex(
      c, dx, r2, &u, &up, &us, &uw, &l, &lp, &ls, &lw, &q, &qp, &qs, &qw);
  if (result < 0) {
    /* the closest_vertex test only found vertices above the intersecting plane
       this would mean that the entire cell lies above the midplane of the line
       segment connecting a point inside the cell (the generator) and a point
       that could be inside or outside the cell (the neighbour). This is
       geometrically absurd and should NEVER happen. */
    voronoi_print_gnuplot_c(c);
    error("Error while searching intersected edge!");
  }
  if (result == 0) {
    /* no intersection */
    return;
  }
  if (result == 2) {
    complicated = 1;
  } else {
    complicated = 0;
  }

  /* At this point:
      up contains a vertex above the plane
      lp contains a vertex below the plane
      us and ls contain the index of the edge that connects up and lp, this edge
      is intersected by the midplane
      u and l contain the projected distances of up and lp to the midplane,
      along dx
     IF complicated is 1, up contains a vertex that is considered to be on the
     plane. All other variables can be considered to be uninitialized in this
     case. */

  int vindex = -1;
  int visitflags[VORONOI3D_MAXNUMVERT];
  int dstack[2 * VORONOI3D_MAXNUMVERT];
  int dstack_size = 0;
  float r = 0.0f;
  int cs = -1, rp = -1;
  int double_edge = 0;
  int i = -1, j = -1, k = -1;

  /* initialize visitflags */
  for (i = 0; i < VORONOI3D_MAXNUMVERT; ++i) {
    visitflags[i] = 0;
  }

  if (complicated) {

    /* We've entered the complicated setup, which means that somewhere along the
       way we found a vertex that is on or very close to the midplane. The index
       of that vertex is stored in up, all other variables are meaningless at
       this point. */

    /* first of all, we need to find a vertex which has edges that extend below
       the plane (since the remainder of our algorithm depends on that). This is
       not necessarily the case: in principle a vertex can only have edges that
       extend inside or above the plane.
       we create a stack of vertices to test (we use dstack for this), and add
       vertex up. For each vertex on the stack, we then traverse its edges. If
       the edge extends above the plane, we ignore it. If it extends below, we
       stop. If the edge lies in the plane, we add the vertex on the other end
       to the stack.
       We make sure that up contains the index of a vertex extending beyond the
       plane on exit. */
    dstack[dstack_size] = up;
    ++dstack_size;
    lw = 0;
    j = 0;
    safewhile(j < dstack_size && lw != -1) {
      up = dstack[j];
      for (i = 0; i < c->orders[up]; ++i) {
        lp = voronoi_get_edge(c, up, i);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
        if (lw == -1) {
          /* jump out of the for loop */
          break;
        }
        if (lw == 0) {
          /* only add each vertex to the stack once */
          k = 0;
          safewhile(k < dstack_size && dstack[k] != lp) { ++k; }
          if (k == dstack_size) {
            dstack[dstack_size] = lp;
            ++dstack_size;
          }
        }
      }
      ++j;
    }

    /* we increased j after lw was calculated, so only the value of lw should be
       used to determine whether or not the loop was successful */
    if (lw != -1) {
      /* we did not find an edge that extends below the plane. There are two
         possible reasons for this: either all vertices of the cell lie above
         or inside the midplane of the segment connecting a point inside the
         cell (the generator) with a point inside or outside the cell (the
         neighbour). This is geometrically absurd.
         Another reason might be that somehow all vertices in the midplane only
         have edges that extend outwards. This is contradictory to the fact that
         a Voronoi cell is convex, and therefore also unacceptable.
         We conclude that we should NEVER end up here. */
      voronoi_print_cell(c);
      error("Unable to find a vertex below the midplane!");
    }
    /* reset the delete stack, we need it later on */
    dstack_size = 0;

    /* the search routine detected a vertex very close to or in the midplane
       the index of this vertex is stored in up
       we proceed by checking the edges of this vertex */

    lp = voronoi_get_edge(c, up, 0);
    lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                             &teststack_size);

    /* the first edge can be below, above or on the plane */
    if (lw != -1) {

      /* above or on the plane: we try to find one below the plane */

      rp = lw;
      i = 1;
      lp = voronoi_get_edge(c, up, i);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      safewhile(lw != -1) {
        ++i;
        if (i == c->orders[up]) {
          /* none of the edges of up is below the plane. Since the cell is
             supposed to be convex, this means the entire cell is above or on
             the plane. This should not happen...
             Furthermore, we should really NEVER end up here, as in this case
             an error should already have be thrown above. */
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

      /* lp, l and lw now contain values corresponding to an edge below the
         plane
         rp contains the result of test_vertex for the first edge of up, for
         reference */

      /* we go on to the next edge of up, and see if we can find an edge that
         does not extend below the plane */

      j = i + 1;
      safewhile(j < c->orders[up] && lw == -1) {
        lp = voronoi_get_edge(c, up, j);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
        ++j;
      }

      if (lw != -1) {
        /* the last iteration increased j by 1 too many, correct this */
        --j;
      }

      /* j-i now contains the number of edges below the plane. We will replace
         up by a new vertex of order this number + 2 (since 2 new edges will be
         created inside the plane)
         however, we do not do this if there is exactly one edge that lies in
         the plane, and all other edges lie below, because in this case we can
         just keep vertex up as is */

      if (j == c->orders[up] && i == 1 && rp == 0) {
        /* keep the order of up, and flag this event for later reference */
        k = c->orders[up];
        double_edge = 1;
      } else {
        /* general case: keep all edges below the plane, and create 2 new ones
           in the plane */
        k = j - i + 2;
      }

      /* create new order k vertex */
      vindex = c->nvert;
      ++c->nvert;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }

      visitflags[vindex] = -vindex;
      /* the new vertex adopts the coordinates of the old vertex */
      c->vertices[3 * vindex + 0] = c->vertices[3 * up + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * up + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * up + 2];

      /* us contains the index of the last edge NOT below the plane
         note that i is at least 1, so there is no need to wrap in this case */
      us = i - 1;

      /* copy all edges of up below the plane into the new vertex, starting from
         edge 1 (edge 0 is reserved to connect to a newly created vertex
         below) */
      k = 1;
      safewhile(i < j) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        /* disconnect up, since this vertex will be removed */
        voronoi_set_edge(c, up, i, -1);
        ++i;
        ++k;
      }

      /* store the index of the first edge not below the plane */
      if (i == c->orders[up]) {
        qs = 0;
      } else {
        qs = i;
      }
    } else { /* if(lw != -1) */

      /* the first edge lies below the plane, try to find one that does not */

      /* we first do a reverse search */
      i = c->orders[up] - 1;
      lp = voronoi_get_edge(c, up, i);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      safewhile(lw == -1) {
        --i;
        if (i == 0) {
          /* No edge above or in the plane found: the cell is unaltered */
          return;
        }
        lp = voronoi_get_edge(c, up, i);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      /* now we do a forward search */
      j = 1;
      qp = voronoi_get_edge(c, up, j);
      qw = voronoi_test_vertex(&c->vertices[3 * qp], dx, r2, &q, teststack,
                               &teststack_size);
      safewhile(qw == -1) {
        ++j;
        qp = voronoi_get_edge(c, up, j);
        qw = voronoi_test_vertex(&c->vertices[3 * qp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      /* at this point j contains the index of the first edge not below the
         plane, i the index of the last edge not below the plane
         we use this to compute the number of edges below the plane. up is
         replaced by a new vertex that has that number + 2 edges (since 2 new
         edges are created inside the plane). We again capture the special event
         where there is only one edge not below the plane, which lies inside the
         plane. In this case up is copied as is. */

      if (i == j && qw == 0) {
        /* we keep up as is, and flag this event */
        double_edge = 1;
        k = c->orders[up];
      } else {
        /* (c->orders[up]-1 - i) + j is the number of edges below the plane */
        k = c->orders[up] - i + j + 1;
      }

      /* create new order k vertex */
      vindex = c->nvert;
      ++c->nvert;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }

      visitflags[vindex] = -vindex;
      /* the new vertex is just a copy of vertex up */
      c->vertices[3 * vindex + 0] = c->vertices[3 * up + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * up + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * up + 2];

      /* as above, us stores the index of the last edge NOT below the plane */
      us = i;

      /* copy all edges below the plane into the new vertex, starting from edge
         1 (edge 0 will be connected to a newly created vertex below)
         We have to do this in two steps: first we copy the high index edges of
         up, then the low index ones (since the edges below the plane are not a
         continuous block of indices in this case) */
      k = 1;
      ++i;
      safewhile(i < c->orders[up]) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        /* disconnect up, it will be removed */
        voronoi_set_edge(c, up, i, -1);
        ++i;
        ++k;
      }
      i = 0;
      safewhile(i < j) {
        qp = voronoi_get_edge(c, up, i);
        qs = voronoi_get_edgeindex(c, up, i);
        voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, i));
        voronoi_set_edge(c, vindex, k, qp);
        voronoi_set_edgeindex(c, vindex, k, qs);
        voronoi_set_edge(c, qp, qs, vindex);
        voronoi_set_edgeindex(c, qp, qs, k);
        voronoi_set_edge(c, up, i, -1);
        ++i;
        ++k;
      }
      /* qs stores the index of the first edge not below the plane */
      qs = j;
    }

    /* at this point, we have created a new vertex that contains all edges of up
       below the plane, and two dangling edges: 0 and k
       Furthermore, us stores the index of the last edge not below the plane,
       qs the index of the first edge not below the plane */

    /* now set the neighbours for the dangling edge(s) */
    if (!double_edge) {
      /* the last edge has the same neighbour as the first edge not below the
         plane */
      voronoi_set_ngb(c, vindex, k, voronoi_get_ngb(c, up, qs));
      /* the first edge has the new neighbour as neighbour */
      voronoi_set_ngb(c, vindex, 0, ngb);
    } else {
      /* up is copied as is, so we also copy its last remaining neighbour */
      voronoi_set_ngb(c, vindex, 0, voronoi_get_ngb(c, up, qs));
    }

    /* add up to the delete stack */
    dstack[dstack_size] = up;
    ++dstack_size;

    /* make sure the variables below have the same meaning as they would have
       if we had the non complicated setup:
       cs contains the index of the last dangling edge of the new vertex
       qp and q correspond to the last vertex that has been deleted
       qs corresponds to the first edge not below the plane
       up and us correspond to the last edge not below the plane, i.e. the edge
       that will be the last one to connect to the new vertex
       note that the value of i is ignored below, it is just used to temporary
       store the new value of up */
    cs = k;
    qp = up;
    q = u;
    i = voronoi_get_edge(c, up, us);
    us = voronoi_get_edgeindex(c, up, us);
    up = i;
    /* we store the index of the newly created vertex in the visitflags of the
       last deleted vertex */
    visitflags[qp] = vindex;
  } else { /* if(complicated) */

    if (u == l) {
      error("Upper and lower vertex are the same!");
    }

    /* the line joining up and lp has general (vector) equation
         x = lp + (up-lp)*t,
       with t a parameter ranging from 0 to 1
       we can rewrite this as
         x = lp*(1-t) + up*t
       the value for t corresponding to the intersection of the line and the
       midplane can be found as the ratio of the projected distance between one
       of the vertices and the midplane, and the total projected distance
       between the two vertices: u-l (remember that u > 0 and l < 0) */
    r = u / (u - l);
    l = 1.0f - r;

    if (r > FLT_MAX || r < -FLT_MAX || l > FLT_MAX || l < -FLT_MAX) {
      error("Value overflow (r: %g, l: %g)", r, l);
    }

    /* create a new order 3 vertex */
    vindex = c->nvert;
    ++c->nvert;
    if (c->nvert == VORONOI3D_MAXNUMVERT) {
      error("Too many vertices!");
    }
    c->orders[vindex] = 3;
    c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
    if (c->offsets[vindex] + 3 >= VORONOI3D_MAXNUMEDGE) {
      error("Too many edges!");
    }

    visitflags[vindex] = -vindex;
    c->vertices[3 * vindex + 0] =
        c->vertices[3 * lp + 0] * r + c->vertices[3 * up + 0] * l;
    c->vertices[3 * vindex + 1] =
        c->vertices[3 * lp + 1] * r + c->vertices[3 * up + 1] * l;
    c->vertices[3 * vindex + 2] =
        c->vertices[3 * lp + 2] * r + c->vertices[3 * up + 2] * l;

    /* add vertex up to the delete stack */
    dstack[dstack_size] = up;
    ++dstack_size;

    /* connect the new vertex to lp (and update lp as well) */
    voronoi_set_edge(c, vindex, 1, lp);
    voronoi_set_edgeindex(c, vindex, 1, ls);
    voronoi_set_edge(c, lp, ls, vindex);
    voronoi_set_edgeindex(c, lp, ls, 1);
    /* disconnect vertex up, it will be deleted */
    voronoi_set_edge(c, up, us, -1);
    /* note that we do not connect edges 0 and 2: edge 2 will be connected to
       the next new vertex that we created, while edge 0 will be connected to
       the last new vertex */

    /* set neighbour relations for the new vertex:
        - edge 0 will be connected to the next intersection point (below), and
          hence has pj as ngb
        - edge 1 is connected to lp and has the original neighbour of the
          intersected edge corresponding to up as neighbour
        - edge 2 has the neighbour on the other side of the original intersected
          edge as neighbour, which is the same as the neighbour of the edge
          corresponding to lp */
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

  /* at this point:
      qp corresponds to the last vertex that has been deleted
      up corresponds to the last vertex that should be used to connect a new
      vertex to the newly created vertex above. In the normal case, qp and up
      are the same vertex, but qp and up can be different if the newly created
      vertex lies in the midplane
      qs contains the index of the edge of qp that is next in line to be tested:
      the edge that comes after the intersected edge that was deleted above
      us corresponds to the edge of up that was connected to the vertex that is
      now connected to the newly created vertex above
      q contains the projected distance between qp and the midplane, along dx
      cs contains the index of the last dangling edge of the last vertex that
      was created above; we still need to connect this edge to a vertex below */

  /* we have found one intersected edge (or at least an edge that lies inside
     the midplane) and created one new vertex that lies in the midplane, with
     dangling edges. We now try to find other intersected edges and create other
     new vertices that will be connected to the first new vertex. */

  int cp = -1;
  int iqs = -1;
  int new_double_edge = -1;

  /* cp and rp both contain the index of the last vertex that was created
     cp will be updated if we add more vertices, rp will be kept, as we need it
     to link the last new vertex to the first new vertex in the end */
  cp = vindex;
  rp = vindex;
  /* we traverse connections of the first removed vertex, until we arrive at an
     edge that links to this vertex (or its equivalent in the degenerate
     case) */
  safewhile(qp != up || qs != us) {
    /* test the next edge of qp */
    lp = voronoi_get_edge(c, qp, qs);
    lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                             &teststack_size);
    if (lw == 0) {

      /* degenerate case: next vertex lies inside the plane */

      k = 2;
      if (double_edge) {
        k = 1;
      }
      /* store the vertex and edge on the other side of the edge in qp and qs */
      qs = voronoi_get_edgeindex(c, qp, qs);
      qp = lp;

      /* move on to the next edge of qp and keep the original edge for
         reference */
      iqs = qs;
      ++qs;
      if (qs == c->orders[qp]) {
        qs = 0;
      }

      /* test the next edges, and try to find one that does NOT lie below the
         plane */
      lp = voronoi_get_edge(c, qp, qs);
      lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                               &teststack_size);
      safewhile(lw == -1) {
        ++k;
        ++qs;
        if (qs == c->orders[qp]) {
          qs = 0;
        }
        lp = voronoi_get_edge(c, qp, qs);
        lw = voronoi_test_vertex(&c->vertices[3 * lp], dx, r2, &l, teststack,
                                 &teststack_size);
      }

      /* qs now contains the next edge NOT below the plane
         k contains the order of the new vertex to create: the number of edges
         below the plane + 2 (+1 if we have a double edge) */

      /* if qp (the vertex in the plane) was already visited before, visitflags
         will contain the index of the newly created vertex that replaces it */
      j = visitflags[qp];

      /* we need to find out what the order of the new vertex will be, and if we
         are dealing with a new double edge or not */
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
                --k;
              } else {
                new_double_edge = 0;
              }
            } else {
              if (j == rp && lp == up && voronoi_get_edge(c, qp, qs) == us) {
                new_double_edge = 1;
                --k;
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
              --k;
            } else {
              new_double_edge = 0;
            }
          } else {
            new_double_edge = 0;
          }
        }
      }

      //      if (j > 0) {
      //        error("Case not handled!");
      //      }

      /* create new order k vertex */
      vindex = c->nvert;
      ++c->nvert;
      if (c->nvert == VORONOI3D_MAXNUMVERT) {
        error("Too many vertices!");
      }
      c->orders[vindex] = k;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      if (c->offsets[vindex] + k >= VORONOI3D_MAXNUMEDGE) {
        error("Too many edges!");
      }

      visitflags[vindex] = -vindex;
      c->vertices[3 * vindex + 0] = c->vertices[3 * qp + 0];
      c->vertices[3 * vindex + 1] = c->vertices[3 * qp + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * qp + 2];

      visitflags[qp] = vindex;
      dstack[dstack_size] = qp;
      ++dstack_size;
      j = vindex;
      i = 0;

      if (!double_edge) {
        voronoi_set_ngb(c, j, i, ngb);
        voronoi_set_edge(c, j, i, cp);
        voronoi_set_edgeindex(c, j, i, cs);
        voronoi_set_edge(c, cp, cs, j);
        voronoi_set_edgeindex(c, cp, cs, i);
        ++i;
      }

      qs = iqs;
      iqs = k - 1;
      if (new_double_edge) {
        iqs = k;
      }
      safewhile(i < iqs) {
        ++qs;
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
        ++i;
      }
      ++qs;
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

      /* normal case: next vertex lies below or above the plane */

      if (lw == 1) {

        /* vertex lies above the plane */

        /* we just delete the vertex and continue with the next edge of this
           vertex */

        qs = voronoi_get_edgeindex(c, qp, qs) + 1;
        if (qs == c->orders[lp]) {
          qs = 0;
        }
        qp = lp;
        q = l;
        dstack[dstack_size] = qp;
        ++dstack_size;
      } else {

        /* vertex lies below the plane */

        /* we have found our next intersected edge: create a new vertex and link
           it to the other vertices */

        if (q == l) {
          error("Upper and lower vertex are the same!");
        }

        r = q / (q - l);
        l = 1.0f - r;

        if (r > FLT_MAX || r < -FLT_MAX || l > FLT_MAX || l < -FLT_MAX) {
          error("Value out of bounds (r: %g, l: %g)!", r, l);
        }

        /* create new order 3 vertex */
        vindex = c->nvert;
        ++c->nvert;
        if (c->nvert == VORONOI3D_MAXNUMVERT) {
          error("Too many vertices!");
        }
        visitflags[vindex] = -vindex;
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

        /* link the edges:
           the first edge is connected to the last edge of the previous new
           vertex. The last edge will be connected to the next new vertex, and
           is left open for the moment */
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

        /* continue with the next edge of qp (the last vertex above the
           midplane */
        ++qs;
        if (qs == c->orders[qp]) {
          qs = 0;
        }
        /* store the last newly created vertex and its dangling edge for the
           next iteration */
        cp = vindex;
        cs = 2;
      } /* if(lw == 1) */

    } /* if(lw == 0) */

  } /* while() */

  /* we finished adding new vertices. Now connect the last dangling edge of the
     last newly created vertex to the first dangling edge of the first newly
     created vertex */
  voronoi_set_edge(c, cp, cs, rp);
  voronoi_set_edge(c, rp, 0, cp);
  voronoi_set_edgeindex(c, cp, cs, 0);
  voronoi_set_edgeindex(c, rp, 0, cs);

  /* now remove the vertices in the delete stack */

  /* the algorithm above did not necessarily visit all vertices above the plane.
     here we scan for vertices that are linked to vertices that are to be
     removed and add them to the delete stack if necessary
     this only works because we made sure that all deleted vertices no longer
     have edges that connect them to vertices that need to stay */
  for (i = 0; i < dstack_size; ++i) {
    for (j = 0; j < c->orders[dstack[i]]; ++j) {
      if (voronoi_get_edge(c, dstack[i], j) >= 0) {
        dstack[dstack_size] = voronoi_get_edge(c, dstack[i], j);
        ++dstack_size;
        voronoi_set_edge(c, dstack[i], j, -1);
        voronoi_set_edgeindex(c, dstack[i], j, -1);
      }
    }
  }

  /* collapse order 1 and 2 vertices: vertices with only 1 edge or 2 edges that
     can be created during the plane intersection routine */
  /* first flag them */
  int low_order_stack[VORONOI3D_MAXNUMVERT];
  int low_order_index = 0;
  for (i = 0; i < c->nvert; ++i) {
    if (voronoi_get_edge(c, i, 0) >= 0 && c->orders[i] < 3) {
      low_order_stack[low_order_index] = i;
      ++low_order_index;
    }
  }

  /* now remove them */
  safewhile(low_order_index) {
    int v = low_order_stack[low_order_index - 1];
    /* the vertex might already have been deleted by a previous operation */
    if (voronoi_get_edge(c, v, 0) < 0) {
      --low_order_index;
      continue;
    }
    if (c->orders[v] == 2) {
      int jj = voronoi_get_edge(c, v, 0);
      int kk = voronoi_get_edge(c, v, 1);
      int bb = voronoi_get_edgeindex(c, v, 1);
      int ll = 0;
      safewhile(ll < c->orders[jj] && voronoi_get_edge(c, jj, ll) != kk) {
        ++ll;
      }
      if (ll == c->orders[jj]) {
        int a = voronoi_get_edgeindex(c, v, 0);
        /* jj and kk are not joined together. Replace their edges pointing to v
           with a new edge pointing from jj to kk */
        voronoi_set_edge(c, jj, a, k);
        voronoi_set_edgeindex(c, jj, a, bb);
        voronoi_set_edge(c, kk, bb, jj);
        voronoi_set_edgeindex(c, kk, bb, a);
        /* no new elements added to the stack: decrease the counter */
        --low_order_index;
      } else {
        /* just remove the edges from jj to v and from kk to v: create two new
           vertices */
        /* vertex jj */
        vindex = c->nvert;
        ++c->nvert;
        c->vertices[3 * vindex] = c->vertices[3 * jj];
        c->vertices[3 * vindex + 1] = c->vertices[3 * jj + 1];
        c->vertices[3 * vindex + 2] = c->vertices[3 * jj + 2];
        c->orders[vindex] = c->orders[jj] - 1;
        c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
        int m = 0;
        for (int n = 0; n < c->orders[jj]; ++n) {
          int lll = voronoi_get_edge(c, jj, n);
          if (lll != v) {
            /* make a new edge */
            voronoi_set_edge(c, vindex, m, lll);
            voronoi_set_edgeindex(c, vindex, m,
                                  voronoi_get_edgeindex(c, jj, n));
            /* update the other vertex */
            voronoi_set_edge(c, lll, voronoi_get_edgeindex(c, jj, n), vindex);
            voronoi_set_edgeindex(c, lll, voronoi_get_edgeindex(c, jj, n), m);
            /* copy ngb information */
            voronoi_set_ngb(c, vindex, m, voronoi_get_ngb(c, jj, n));
            ++m;
          }
          /* remove the old vertex */
          voronoi_set_edge(c, jj, n, -1);
          voronoi_set_edgeindex(c, jj, n, -1);
        }
        /* vertex kk */
        vindex = c->nvert;
        ++c->nvert;
        c->vertices[3 * vindex] = c->vertices[3 * kk];
        c->vertices[3 * vindex + 1] = c->vertices[3 * kk + 1];
        c->vertices[3 * vindex + 2] = c->vertices[3 * kk + 2];
        c->orders[vindex] = c->orders[kk] - 1;
        c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
        m = 0;
        for (int n = 0; n < c->orders[kk]; ++n) {
          int lll = voronoi_get_edge(c, kk, n);
          if (lll != v) {
            /* make a new edge */
            voronoi_set_edge(c, vindex, m, lll);
            voronoi_set_edgeindex(c, vindex, m,
                                  voronoi_get_edgeindex(c, kk, n));
            /* update the other vertex */
            voronoi_set_edge(c, lll, voronoi_get_edgeindex(c, kk, n), vindex);
            voronoi_set_edgeindex(c, lll, voronoi_get_edgeindex(c, kk, n), m);
            /* copy ngb information */
            /* this one is special: we copy the ngb corresponding to the
               deleted edge and skip the one after that */
            if (n == bb + 1) {
              voronoi_set_ngb(c, vindex, m, voronoi_get_ngb(c, kk, bb));
            } else {
              voronoi_set_ngb(c, vindex, m, voronoi_get_ngb(c, kk, n));
            }
            ++m;
          }
          /* remove the old vertex */
          voronoi_set_edge(c, kk, n, -1);
          voronoi_set_edgeindex(c, kk, n, -1);
        }
        /* check if jj or kk has become an order 2 vertex */
        /* if they have become an order 1 vertex, they were already an order 2
           vertex, and they should already be in the list... */
        if (c->orders[vindex] == 2) {
          if (c->orders[vindex - 1] == 2) {
            low_order_stack[low_order_index] = vindex - 1;
            ++low_order_index;
            low_order_stack[low_order_index] = vindex;
            /* we do not increase the index here: we want this element to be the
               next element that is processed */
          } else {
            low_order_stack[low_order_index] = vindex;
          }
        } else {
          if (c->orders[vindex - 1] == 2) {
            low_order_stack[low_order_index] = vindex - 1;
          } else {
            /* no new vertices added to the stack: decrease the counter */
            --low_order_index;
          }
        }
      }
      /* Remove the vertex */
      voronoi_set_edge(c, v, 0, -1);
      voronoi_set_edgeindex(c, v, 0, -1);
      voronoi_set_edge(c, v, 1, -1);
      voronoi_set_edgeindex(c, v, 1, -1);
    } else if (c->orders[v] == 1) {
      int jj = voronoi_get_edge(c, v, 0);
      /* we have to remove the edge between j and v. We create a new vertex */
      vindex = c->nvert;
      ++c->nvert;
      c->vertices[3 * vindex] = c->vertices[3 * jj];
      c->vertices[3 * vindex + 1] = c->vertices[3 * jj + 1];
      c->vertices[3 * vindex + 2] = c->vertices[3 * jj + 2];
      c->orders[vindex] = c->orders[j] - 1;
      c->offsets[vindex] = c->offsets[vindex - 1] + c->orders[vindex - 1];
      int m = 0;
      for (int kk = 0; kk < c->orders[j]; ++kk) {
        int ll = voronoi_get_edge(c, jj, kk);
        if (ll != v) {
          /* make a new edge */
          voronoi_set_edge(c, vindex, m, ll);
          voronoi_set_edgeindex(c, vindex, m, voronoi_get_edgeindex(c, jj, kk));
          /* update the other vertex */
          voronoi_set_edge(c, ll, voronoi_get_edgeindex(c, jj, kk), vindex);
          voronoi_set_edgeindex(c, ll, voronoi_get_edgeindex(c, jj, kk), m);
          /* copy ngb information */
          voronoi_set_ngb(c, vindex, m, voronoi_get_ngb(c, jj, kk));
          ++m;
        }
        /* remove the old vertex */
        voronoi_set_edge(c, jj, kk, -1);
        voronoi_set_edgeindex(c, jj, kk, -1);
      }
      /* if the new vertex is a new order 2 vertex, add it to the stack */
      if (c->orders[vindex] == 2) {
        low_order_stack[low_order_index - 1] = vindex;
      } else {
        --low_order_index;
      }
      /* remove the order 1 vertex */
      voronoi_set_edge(c, v, 0, -1);
      voronoi_set_edgeindex(c, v, 0, -1);
    } else {
      error("Vertex with order %i. This should not happen!", c->orders[v]);
    }
  }

  /* remove deleted vertices from all arrays */
  struct voronoi_cell new_cell;
  /* make sure the contents of the new cell are the same as for the old cell */
  memcpy(&new_cell, c, sizeof(struct voronoi_cell));
  int m, n;
  for (vindex = 0; vindex < c->nvert; ++vindex) {
    j = vindex;
    /* find next edge that is not deleted */
    safewhile(j < c->nvert && voronoi_get_edge(c, j, 0) < 0) { ++j; }

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
    for (k = 0; k < c->orders[j]; ++k) {
      voronoi_set_edge(&new_cell, vindex, k, voronoi_get_edge(c, j, k));
      voronoi_set_edgeindex(&new_cell, vindex, k,
                            voronoi_get_edgeindex(c, j, k));
      voronoi_set_ngb(&new_cell, vindex, k, voronoi_get_ngb(c, j, k));
    }

    /* update other edges */
    for (k = 0; k < c->orders[j]; ++k) {
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

#ifdef VORONOI3D_EXPENSIVE_CHECKS
  voronoi_check_cell_consistency(c);
#endif
}

/**
 * @brief Get the volume of the tetrahedron made up by the four given vertices.
 *
 * The vertices are not expected to be oriented in a specific way. If the input
 * happens to be coplanar or colinear, the returned volume will just be zero.
 *
 * @param v1 First vertex.
 * @param v2 Second vertex.
 * @param v3 Third vertex.
 * @param v4 Fourth vertex.
 * @return Volume of the tetrahedron.
 */
__attribute__((always_inline)) INLINE float voronoi_volume_tetrahedron(
    const float *v1, const float *v2, const float *v3, const float *v4) {

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
 * @brief Get the centroid of the tetrahedron made up by the four given
 * vertices.
 *
 * The centroid is just the average of four vertex coordinates.
 *
 * @param centroid Array to store the centroid in.
 * @param v1 First vertex.
 * @param v2 Second vertex.
 * @param v3 Third vertex.
 * @param v4 Fourth vertex.
 */
__attribute__((always_inline)) INLINE void voronoi_centroid_tetrahedron(
    float *centroid, const float *v1, const float *v2, const float *v3,
    const float *v4) {

  centroid[0] = 0.25f * (v1[0] + v2[0] + v3[0] + v4[0]);
  centroid[1] = 0.25f * (v1[1] + v2[1] + v3[1] + v4[1]);
  centroid[2] = 0.25f * (v1[2] + v2[2] + v3[2] + v4[2]);
}

/**
 * @brief Calculate the volume and centroid of a 3D Voronoi cell.
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
  for (i = 1; i < cell->nvert; ++i) {

    v2[0] = cell->vertices[3 * i + 0];
    v2[1] = cell->vertices[3 * i + 1];
    v2[2] = cell->vertices[3 * i + 2];

    /*  loop over the edges of the vertex*/
    for (j = 0; j < cell->orders[i]; ++j) {

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

        int loopcount = 0;
        safewhile(m != i) {
          if (loopcount == 999) {
            voronoi_print_cell(cell);
            voronoi_print_gnuplot_c(cell);
          }
          ++loopcount;
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
 * Note that the face midpoints are calculated relative w.r.t. the cell
 * generator!
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

        safewhile(m != i) {
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
        ++cell->nface;

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
 * @brief Initialize a 3D Voronoi cell.
 *
 * @param cell 3D Voronoi cell to initialize.
 * @param x Position of the generator of the cell.
 * @param anchor Anchor of the simulation box.
 * @param side Side lengths of the simulation box.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, const double *x, const double *anchor,
    const double *side) {

  cell->x[0] = x[0];
  cell->x[1] = x[1];
  cell->x[2] = x[2];

  voronoi_initialize(cell, anchor, side);

  cell->volume = 0.0f;
  cell->centroid[0] = 0.0f;
  cell->centroid[1] = 0.0f;
  cell->centroid[2] = 0.0f;
  cell->nface = 0;
}

/**
 * @brief Interact a 3D Voronoi cell with a particle with given relative
 * position and ID.
 *
 * @param cell 3D Voronoi cell.
 * @param dx Relative position of the interacting generator w.r.t. the cell
 * generator (in fact: dx = generator - neighbour).
 * @param id ID of the interacting neighbour.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, const float *dx, unsigned long long id) {

  voronoi_intersect(cell, dx, id);
}

/**
 * @brief Finalize a 3D Voronoi cell.
 *
 * @param cell 3D Voronoi cell.
 * @return Maximal radius that could still change the structure of the cell.
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
 * cell and the given neighbour.
 *
 * @param cell 3D Voronoi cell.
 * @param ngb ID of a particle that is possibly a neighbour of this cell.
 * @param midpoint Array to store the relative position of the face in.
 * @return 0 if the given neighbour is not a neighbour, the surface area of
 * the face otherwise.
 */
__attribute__((always_inline)) INLINE float voronoi_get_face(
    const struct voronoi_cell *cell, unsigned long long ngb, float *midpoint) {

  int i = 0;
  while (i < cell->nface && cell->ngbs[i] != ngb) {
    ++i;
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
 * @brief Get the centroid of a 3D Voronoi cell.
 *
 * @param cell 3D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    const struct voronoi_cell *cell, float *centroid) {

  centroid[0] = cell->centroid[0];
  centroid[1] = cell->centroid[1];
  centroid[2] = cell->centroid[2];
}

#endif  // SWIFT_VORONOIXD_ALGORITHM_H
