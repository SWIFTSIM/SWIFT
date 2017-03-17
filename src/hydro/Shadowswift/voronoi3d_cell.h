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

#ifndef SWIFT_VORONOIXD_CELL_H
#define SWIFT_VORONOIXD_CELL_H

/* Maximal number of neighbours that can be stored in a voronoi_cell struct */
#define VORONOI3D_MAXNUMNGB 100
/* Maximal number of vertices that can be stored in a voronoi_cell struct */
#define VORONOI3D_MAXNUMVERT 500
/* Maximal number of edges that can be stored in a voronoi_cell struct */
#define VORONOI3D_MAXNUMEDGE 1500
/* Maximal number of faces that can be stored in a voronoi_cell struct */
#define VORONOI3D_MAXFACE 100

/* 3D Voronoi cell */
struct voronoi_cell {

  /* The position of the generator of the cell. */
  double x[3];

  /* The volume of the 3D cell. */
  float volume;

  /* The centroid of the cell. */
  float centroid[3];

  /* Number of cell vertices. */
  int nvert;

  /* Vertex coordinates. */
  float vertices[3 * VORONOI3D_MAXNUMVERT];

  /* Number of edges for every vertex. */
  char orders[VORONOI3D_MAXNUMVERT];

  /* Offsets of the edges, edgeindices and neighbours corresponding to a
     particular vertex in the internal arrays */
  int offsets[VORONOI3D_MAXNUMVERT];

  /* Edge information. Edges are ordered counterclockwise w.r.t. a vector
     pointing from the cell generator to the vertex. */
  int edges[VORONOI3D_MAXNUMEDGE];

  /* Additional edge information. */
  char edgeindices[VORONOI3D_MAXNUMEDGE];

  /* Neighbour information. This field is used differently depending on where we
     are in the algorithm. During cell construction, it contains, for every edge
     of every vertex, the index of the neighbour that generates the face
     counterclockwise of the edge w.r.t. a vector pointing from the vertex along
     the edge. After cell finalization, it contains a neighbour for every face,
     in the same order as the face_areas and face_midpoints arrays. */
  unsigned long long ngbs[VORONOI3D_MAXNUMEDGE];

  /* Number of faces of the cell. */
  unsigned char nface;

  /* Surface areas of the cell faces. */
  float face_areas[VORONOI3D_MAXFACE];

  /* Midpoints of the cell faces. */
  float face_midpoints[VORONOI3D_MAXFACE][3];
};

/**
 * @brief Copy the contents of the 3D Voronoi cell pointed to by source into the
 * 3D Voronoi cell pointed to by destination
 *
 * @param source Pointer to a 3D Voronoi cell to read from.
 * @param destination Pointer to a 3D Voronoi cell to write to.
 */
__attribute__((always_inline)) INLINE void voronoi3d_cell_copy(
    struct voronoi_cell *source, struct voronoi_cell *destination) {

  /* Copy the position of the generator of the cell. */
  destination->x[0] = source->x[0];
  destination->x[1] = source->x[1];
  destination->x[2] = source->x[2];

  /* Copy the volume of the 3D cell. */
  destination->volume = source->volume;

  /* Copy the centroid of the cell. */
  destination->centroid[0] = source->centroid[0];
  destination->centroid[1] = source->centroid[1];
  destination->centroid[2] = source->centroid[2];

  /* Copy the number of cell vertices. */
  destination->nvert = source->nvert;

  /* Copy the vertex coordinates. We only copy the 3*nvert first coordinates. */
  for (int i = 0; i < 3 * source->nvert; ++i) {
    destination->vertices[i] = source->vertices[i];
  }

  /* Copy the number of edges for every vertex. Again, we only copy the nvert
     first values. */
  for (int i = 0; i < source->nvert; ++i) {
    destination->orders[i] = source->orders[i];
  }

  /* Copy the nvert first values of the offsets. */
  for (int i = 0; i < source->nvert; ++i) {
    destination->offsets[i] = source->offsets[i];
  }

  /* Copy the edge information. No idea how many edges we have, so we copy
     everything. */
  for (int i = 0; i < VORONOI3D_MAXNUMEDGE; ++i) {
    destination->edges[i] = source->edges[i];
  }

  /* Copy all additional edge information. */
  for (int i = 0; i < VORONOI3D_MAXNUMEDGE; ++i) {
    destination->edgeindices[i] = source->edgeindices[i];
  }

  /* Copy neighbour information. Since neighbours are stored per edge, the total
     number of neighbours in this list is larger than numngb and we copy
     everything. */
  for (int i = 0; i < VORONOI3D_MAXNUMEDGE; ++i) {
    destination->ngbs[i] = source->ngbs[i];
  }
}

#endif  // SWIFT_VORONOIXD_CELL_H
