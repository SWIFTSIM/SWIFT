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

#include "error.h"
#include "inline.h"
#include "voronoi1d_cell.h"

#include <math.h>
#include <stdlib.h>

/**
 * @brief Store the extents of the simulation box in the global variables.
 *
 * @param anchor Corner of the simulation box with the lowest coordinate values.
 * @param side Side lengths of the simulation box.
 */
__attribute__((always_inline)) INLINE static void voronoi_set_box(
    const float *anchor, const float *side) {}

/**
 * @brief Initialize a 1D Voronoi cell.
 *
 * Sets the positions of left and right neighbours to very large values, the
 * generator position to the given particle position, and all other quantities
 * to zero.
 *
 * @param cell 1D Voronoi cell to initialize.
 * @param x Position of the generator of the cell.
 * @param anchor Anchor of the simulation box.
 * @param side Side lengths of the simulation box.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_init(
    struct voronoi_cell *cell, const double *x, const double *anchor,
    const double *side) {
  cell->x = x[0];
  cell->xL = anchor[0] - cell->x;
  cell->xR = anchor[0] + side[0] - cell->x;
  cell->idL = 0;
  cell->idR = 0;
  cell->volume = 0.0f;
  cell->centroid = 0.0f;
}

/**
 * @brief Interact a 1D Voronoi cell with a particle with given relative
 * position and ID.
 *
 * This method checks if the given relative position is closer to the cell
 * generator than the current left or right neighbour and updates neighbours
 * accordingly.
 *
 * @param cell 1D Voronoi cell.
 * @param dx Relative position of the interacting generator w.r.t. the cell
 * generator (in fact: dx = generator - neighbour).
 * @param id ID of the interacting neighbour.
 */
__attribute__((always_inline)) INLINE void voronoi_cell_interact(
    struct voronoi_cell *cell, const float *dx, unsigned long long id) {

  /* Check for stupidity */
  if (dx[0] == 0.0f) {
    error("Cannot interact a Voronoi cell generator with itself!");
  }

  if (-dx[0] < 0.0f) {
    /* New left neighbour? */
    if (-dx[0] > cell->xL) {
      cell->xL = -dx[0];
      cell->idL = id;
    }
  } else {
    /* New right neighbour? */
    if (-dx[0] < cell->xR) {
      cell->xR = -dx[0];
      cell->idR = id;
    }
  }
}

/**
 * @brief Finalize a 1D Voronoi cell.
 *
 * Calculates the relative positions of the midpoints of the faces (which in
 * this case are just the midpoints of the segments connecting the generator
 * with the two neighbours) w.r.t. the generator, and the cell volume (length)
 * and centroid (midpoint of the segment connecting the midpoints of the faces).
 * This function returns the maximal radius at which a particle could still
 * change the structure of the cell, i.e. twice the largest distance between
 * the cell generator and one of its faces. If the cell has been interacted with
 * all neighbours within this radius, we know for sure that the cell is
 * complete.
 *
 * @param cell 1D Voronoi cell.
 * @return Maximal radius that could still change the structure of the cell.
 */
__attribute__((always_inline)) INLINE float voronoi_cell_finalize(
    struct voronoi_cell *cell) {

  float xL, xR;
  float max_radius;

  max_radius = fmax(-cell->xL, cell->xR);
  cell->xL = xL = 0.5f * cell->xL;
  cell->xR = xR = 0.5f * cell->xR;

  cell->volume = xR - xL;
  cell->centroid = cell->x + 0.5f * (xL + xR);

  return max_radius;
}

/**
 * @brief Get the oriented surface area and midpoint of the face between a
 * 1D Voronoi cell and the given neighbour.
 *
 * This function also checks if the given neighbour is in fact a neighbour of
 * this cell. Since we perform gradient and flux calculations for all neighbour
 * pairs within the smoothing length, which assumes the cell to be spherical,
 * it can happen that this is not the case. It is the responsibility of the
 * routine that calls this function to check for a zero return value and
 * deal with it appropriately.
 *
 * For this specific case, we simply check if the neighbour is the left or
 * right neighbour and set the surface area to 1. The midpoint is set to the
 * relative position vector of the appropriate face.
 *
 * @param cell 1D Voronoi cell.
 * @param ngb ID of a particle that is possibly a neighbour of this cell.
 * @param midpoint Array to store the relative position of the face in.
 * @return 0 if the given neighbour is not a neighbour, surface area 1.0f
 * otherwise.
 */
__attribute__((always_inline)) INLINE float voronoi_get_face(
    const struct voronoi_cell *cell, unsigned long long ngb, float *midpoint) {

  if (ngb != cell->idL && ngb != cell->idR) {
    /* this is perfectly possible: we interact with all particles within the
       smoothing length, and they do not need to be all neighbours.
       If this happens, we return 0, so that the flux method can return */
    return 0.0f;
  }

  if (ngb == cell->idL) {
    /* Left face */
    midpoint[0] = cell->xL;
  } else {
    /* Right face */
    midpoint[0] = cell->xR;
  }
  /* The other components of midpoint are just zero */
  midpoint[1] = 0.0f;
  midpoint[2] = 0.0f;

  return 1.0f;
}

/**
 * @brief Get the centroid of a 1D Voronoi cell.
 *
 * We store only the relevant coordinate of the centroid, but need to return
 * a 3D vector.
 *
 * @param cell 1D Voronoi cell.
 * @param centroid Array to store the centroid in.
 */
__attribute__((always_inline)) INLINE void voronoi_get_centroid(
    const struct voronoi_cell *cell, float *centroid) {

  centroid[0] = cell->centroid;
  centroid[1] = 0.0f;
  centroid[2] = 0.0f;
}

#endif  // SWIFT_VORONOIXD_ALGORITHM_H
