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
#ifndef SWIFT_ZOOM_CELL_H
#define SWIFT_ZOOM_CELL_H

/* Local includes */
#include "space.h"
#include "timeline.h"

/* Zoom specific cell_getid to handle different levels. */
void zoom_construct_tl_cells(struct space *s, const integertime_t ti_current,
                             int verbose);

/**
 * @brief Is this cell within the buffer region?
 *
 * @param c The #cell.
 * @param s The #space.
 */
__attribute__((always_inline)) INLINE static int cell_inside_buffer_region(
    const struct cell *c, const struct space *s) {

  /* Get the middle of the cell (since the cell grids align this eliminates
   * any issues from rounding). */
  const double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                         c->loc[1] + 0.5 * c->width[1],
                         c->loc[2] + 0.5 * c->width[2]};

  return ((mid[0] > s->zoom_props->buffer_lower_bounds[0]) &&
          (mid[0] < s->zoom_props->buffer_upper_bounds[0]) &&
          (mid[1] > s->zoom_props->buffer_lower_bounds[1]) &&
          (mid[1] < s->zoom_props->buffer_upper_bounds[1]) &&
          (mid[2] > s->zoom_props->buffer_lower_bounds[2]) &&
          (mid[2] < s->zoom_props->buffer_upper_bounds[2]));
}

/**
 * @brief Is this cell within the zoom region?
 *
 * @param c The #cell.
 * @param s The #space.
 */
__attribute__((always_inline)) INLINE static int cell_inside_zoom_region(
    const struct cell *c, const struct space *s) {

  /* Get the middle of the cell (since the cell grids align this eliminates
   * any issues from rounding). */
  const double mid[3] = {c->loc[0] + 0.5 * c->width[0],
                         c->loc[1] + 0.5 * c->width[1],
                         c->loc[2] + 0.5 * c->width[2]};

  return ((mid[0] > s->zoom_props->region_lower_bounds[0]) &&
          (mid[0] < s->zoom_props->region_upper_bounds[0]) &&
          (mid[1] > s->zoom_props->region_lower_bounds[1]) &&
          (mid[1] < s->zoom_props->region_upper_bounds[1]) &&
          (mid[2] > s->zoom_props->region_lower_bounds[2]) &&
          (mid[2] < s->zoom_props->region_upper_bounds[2]));
}

#endif /* SWIFT_ZOOM_CELL_H */
