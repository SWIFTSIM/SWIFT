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

#ifndef SWIFT_VORONOIXD_CELL_H
#define SWIFT_VORONOIXD_CELL_H

/* Maximal number of vertices (and neighbours) that can be stored in a
   voronoi_cell struct. */
#define VORONOI2D_MAXNUMVERT 100

/* 2D Voronoi cell */
struct voronoi_cell {

  /* The position of the generator of the cell. */
  double x[2];

  /* The "volume" of the 2D cell. */
  float volume;

  /* The centroid of the cell. */
  float centroid[2];

  /* Number of cell vertices (and neighbours). */
  int nvert;

  /* We only need to store one of these at the same time. */
  union {
    /* The relative positions of the vertices of the cell. */
    float vertices[VORONOI2D_MAXNUMVERT][2];

    /* The midpoints of the faces. */
    float face_midpoints[VORONOI2D_MAXNUMVERT][2];
  };

  /* The ids of the neighbouring cells. */
  unsigned long long ngbs[VORONOI2D_MAXNUMVERT];

  /* The lengths of the faces. */
  float face_lengths[VORONOI2D_MAXNUMVERT];
};

#endif  // SWIFT_VORONOIXD_CELL_H
