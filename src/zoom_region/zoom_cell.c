/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *               2024 Will J. Roper (w.roper@sussex.ac.uk)
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

/* Includes */
#include <float.h>

/* Local includes */
#include "../cell.h"
#include "../space.h"
#include "zoom_cell.h"

/**
 * @brief For a given particle location, what TL cell does it belong in nested
 *        grid?
 *
 * This is a simple helper function to reduce repetition.
 *
 * @param cdim The cell grid dimensions.
 * @param bounds The edges of this nested region.
 * @param x, y, z Location of particle.
 * @param iwidth The width of a cell in this grid.
 * @param offset The offset of this cell type in cells_top.
 */
int cell_getid_with_bounds(const int *cdim, const double *bounds,
                           const double x, const double y, const double z,
                           const double *iwidth, const int offset) {

  /* Get the cell ijk coordinates in this grid. */
  const int i = (x - bounds[0]) * iwidth[0];
  const int j = (y - bounds[1]) * iwidth[1];
  const int k = (z - bounds[2]) * iwidth[2];

  /* Which zoom TL cell are we in? */
  const int cell_id = cell_getid_offset(cdim, offset, i, j, k);

  return cell_id;
}

/**
 * @brief For a given particle location, what TL cell does it belong to?
 *
 * Slightly more complicated in the zoom case, as there are now 1/2 embedded TL
 * grids.
 *
 * First see if the particle is in the background grid, if it is an empty or
 * void cell check the nested cell types.
 *
 * @param s The space.
 * @param x, y, z Location of particle.
 */
int zoom_cell_getid(const struct space *s, const double x, const double y,
                    const double z) {

  /* Initilaise the cell id to return. */
  int cell_id;

  /* Lets get some properties of the zoom region. */
  const struct zoom_region_properties *zoom_props = s->zoom_props;
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_lower_bounds[3] = {zoom_props->region_lower_bounds[0],
                                       zoom_props->region_lower_bounds[1],
                                       zoom_props->region_lower_bounds[2]};
  const int buffer_cell_offset = zoom_props->buffer_cell_offset;
  const double buffer_lower_bounds[3] = {zoom_props->buffer_lower_bounds[0],
                                         zoom_props->buffer_lower_bounds[1],
                                         zoom_props->buffer_lower_bounds[2]};

  /* Here we go down the heirarchy to get the cell_id, it's marginally slower
   * but guarantees that void cells are handled properly. */

  /* Get the background cell ijk coordinates. */
  const int bkg_i = x * s->iwidth[0];
  const int bkg_j = y * s->iwidth[1];
  const int bkg_k = z * s->iwidth[2];

  /* Which background cell is this? */
  cell_id = cell_getid(s->cdim, bkg_i, bkg_j, bkg_k) + bkg_cell_offset;

  /* If this is a void cell we are in the zoom region. */
  if (s->cells_top[cell_id].subtype == void_cell) {

    /* Which zoom TL cell are we in? */
    cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_lower_bounds, x,
                                     y, z, s->zoom_props->iwidth,
                                     /*offset*/ 0);

  }

  /* If this is an empty cell we are in the buffer cells.
   * Otherwise, It's a legitimate background cell, and we'll return it. */
  else if (s->cells_top[cell_id].subtype == empty) {

    /* Which buffer TL cell are we in? */
    cell_id = cell_getid_with_bounds(
        s->zoom_props->buffer_cdim, buffer_lower_bounds, x, y, z,
        s->zoom_props->buffer_iwidth, buffer_cell_offset);

    /* Here we need to check if this is the void buffer cell.
     * Otherwise, It's a legitimate buffer cell, and we'll return it. */
    if (s->cells_top[cell_id].subtype == void_cell) {

      /* Which zoom TL cell are we in? */
      cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_lower_bounds,
                                       x, y, z, s->zoom_props->iwidth,
                                       /*offset*/ 0);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (cell_id < 0 || cell_id >= s->nr_cells)
    error("cell_id out of range: %i (%f %f %f)", cell_id, x, y, z);
#endif

  return cell_id;
}
