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

/* 2D Voronoi cell */
struct voronoi_cell {

  /* The position of the generator of the cell. */
  double x[2];

  /* The "volume" of the 2D cell. */
  float volume;

  /* The centroid of the cell. */
  float centroid[2];
};

#endif  // SWIFT_VORONOIXD_CELL_H
