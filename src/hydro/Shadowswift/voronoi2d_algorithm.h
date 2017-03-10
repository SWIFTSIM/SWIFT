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
#include "voronoi2d_cell.h"

/**
 * @brief Initialize a 2D Voronoi cell
 *
 * @param cell 2D Voronoi cell to initialize.
 * @param x Position of the generator of the cell.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, double *x) {}

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

  return 1.0f;
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

  return 0.0f;
}

/**
 * @brief Get the centroid of a 2D Voronoi cell
 *
 * @param cell 2D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    struct voronoi_cell *cell, float *centroid) {}

#endif  // SWIFT_VORONOIXD_ALGORITHM_H
