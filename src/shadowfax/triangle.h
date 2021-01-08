/*******************************************************************************
 * This file is part of cVoronoi.
 * Copyright (c) 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

/**
 * @file triangle.h
 *
 * @brief Triangle object and functionality.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#ifndef SWIFT_TRIANGLE_H
#define SWIFT_TRIANGLE_H

/**
 * @brief Triangle.
 *
 * A triangle connects 3 points in 2D space, called the vertices of the
 * triangle. The connections between these vertices are line segments. Within
 * a Delaunay tessellation, each line segment forms the border between two
 * neighbouring triangles.
 *
 * Within this struct, we store the indices of the 3 vertices, as well as the
 * indices of the neighbouring triangles and the indices of this triangle within
 * the neighbour lists of its neighbours. The latter is done to speed up lookup
 * operations during incremental Delaunay tessellation construction.
 *
 * We use the following convention: neighbours[0] corresponds to the neighbour
 * across the line segment that connects vertices[1] and vertices[2].
 * index_in_neighbours[0] is then the position of the index of this triangle
 * within neighbours[0]->neighbours, so that
 * neighbours[0]->neighbours[index_in_neighbours[0]] is this triangle again.
 */
struct triangle {

  /*! @brief Indices of the particles that make up the triangle. */
  int vertices[3];

  /*! @brief Indices of the neighbour triangles. The neighbour at index i is the
   *  triangle on the other side of the segment opposite vertex i. */
  int neighbours[3];

  /*! @brief Index of this triangle in the neighbour list of its neighbours. */
  int index_in_neighbour[3];
};

/**
 * @brief Initialize a triangle with the given vertices.
 *
 * @param t Triangle.
 * @param v0, v1, v2 Vertices.
 */
inline static void triangle_init(struct triangle* restrict t, int v0, int v1,
                                 int v2) {

  t->vertices[0] = v0;
  t->vertices[1] = v1;
  t->vertices[2] = v2;

  t->neighbours[0] = -1;
  t->neighbours[1] = -1;
  t->neighbours[2] = -1;

  t->index_in_neighbour[0] = -1;
  t->index_in_neighbour[1] = -1;
  t->index_in_neighbour[2] = -1;
}

/**
 * @brief Replace the neighbour at the given index with the given new neighbour.
 *
 * @param index Index of the neighbour in the list of neighbours.
 * @param neighbour New neighbour.
 * @param index_in_neighbour Index of this triangle in the neighbour list of the
 * new neighbour.
 */
inline static void triangle_swap_neighbour(struct triangle* restrict t,
                                           int index, int neighbour,
                                           int index_in_neighbour) {
  t->neighbours[index] = neighbour;
  t->index_in_neighbour[index] = index_in_neighbour;
}

#endif /* SWIFT_TRIANGLE_H */
