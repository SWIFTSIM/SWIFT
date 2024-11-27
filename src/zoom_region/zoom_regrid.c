/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "scheduler.h"
#include "space.h"
#include "zoom.h"

/**
 * @breif Check whether we need a regrid based on the zoom region.
 *
 * This function is called in space_regrid in space_regrid.c and provides
 * the zoom specific regrid check.
 *
 * TODO: In the future this check needs to take into account the particle
 * distribution to ensure that the zoom region is sufficiently placed given
 * any bulk motion.
 *
 * @param s The #space.
 * @param new_cdim The new top-level cell dimensions (based on current hmax).
 *
 * @return 1 if a zoom regrid is needed, 0 otherwise.
 */
int zoom_need_regrid(const struct space *s, const int new_cdim[3]) {
  /* If we are running a zoom do we need to regrid based on the new cdim? */
  return (new_cdim[0] < s->zoom_props->cdim[0] ||
          new_cdim[1] < s->zoom_props->cdim[1] ||
          new_cdim[2] < s->zoom_props->cdim[2]);
}

/**
 * @brief Find an acceptable geometry given the required zoom cdim.
 *
 * Unlike a uniform box, we have options with a zoom region. We can either
 * increase the background cdim to better resolve the zoom region or we can
 * increase the depth of the zoom region.
 *
 * We'll first try to increase the background cdim a reasonable amount since
 * this will carry less of a performance penalty for small increases than
 * doubling the number of zoom cells. This is also most likely to produce a
 * valid set up in all but the most extreme cases.
 *
 * @param s The #space.
 * @param new_cdim The new top-level cell dimensions (based on current hmax).
 */
void zoom_regrid_find_acceptable_geometry(struct space *s,
                                          const int new_cdim[3]) {

  /* Loop until we've found an acceptable geometry. */
  while (zoom_need_regrid(s, new_cdim)) {

    /* First try increasing the background cdim to a maximum increase of 50% its
     * current value. */
    int old_bkg_cdim = s->cdim[0];
    while (zoom_need_regrid(s, new_cdim) && s->cdim[0] < 1.5 * old_bkg_cdim) {

      /* Increment the background cdim. */
      s->cdim[0]++;
      s->cdim[1]++;
      s->cdim[2]++;

      /* Recalculate the zoom region geometry. (silently) */
      zoom_region_init(s, 0);
    }

    /* If this worked we can stop here. */
    if (!zoom_need_regrid(s, new_cdim)) {
      break;
    }

    /* If it didn't work we'll try increasing the depth of the zoom region (with
     * the original background cdim). */
    s->zoom_props->zoom_cell_depth++;
    zoom_region_init(s, 0);
  }
}

/**
 * @brief Prepare the cells for the zoom region.
 *
 * This function is called in space_regrid in space_regrid.c and provides
 * the zoom specific cell preparation, freeing only the zoom specific cell
 * pointer arrays and computing the new zoom region geometry if necessary.
 *
 * The cells are counted here and stored on the space. In the zoom case this
 * includes counts from each top-level cell grid (zoom, bkg, and buffer if
 * used).
 *
 * @param s The #space.
 * @param zoom_cdim The new top-level cell dimensions (based on current hmax).
 * @param verbose Whether to print verbose output.
 */
void zoom_prepare_cells(struct space *s, const int zoom_cdim[3], int verbose) {
  /* Free the old zoom specific cells, if they were allocated. */
  if (s->cells_top != NULL) {
    swift_free("local_zoom_cells_top", s->zoom_props->local_zoom_cells_top);
    swift_free("local_bkg_cells_top", s->zoom_props->local_bkg_cells_top);
    swift_free("local_zoom_cells_with_particles_top",
               s->zoom_props->local_zoom_cells_with_particles_top);
    swift_free("local_bkg_cell_with_particless_top",
               s->zoom_props->local_bkg_cells_with_particles_top);
    swift_free("local_buffer_cell_with_particless_top",
               s->zoom_props->local_buffer_cells_with_particles_top);
    swift_free("void_cell_indices", s->zoom_props->void_cell_indices);
    swift_free("neighbour_cells_top", s->zoom_props->neighbour_cells_top);

    /* Find an acceptable geometry given the required zoom cdim. */
    zoom_regrid_find_acceptable_geometry(s, zoom_cdim);

    /* The above function found the geometry silently, if we're running in
     * verbose mode print the cell properties report. */
    if (verbose) {
      zoom_report_cell_properties(s);
    }
  }

  /* Count the number of top level cells. */
  s->tot_cells = s->nr_cells =
      (s->cdim[0] * s->cdim[1] * s->cdim[2]) +
      (s->zoom_props->cdim[0] * s->zoom_props->cdim[1] *
       s->zoom_props->cdim[2]) +
      (s->zoom_props->buffer_cdim[0] * s->zoom_props->buffer_cdim[1] *
       s->zoom_props->buffer_cdim[2]);
}

/**
 * @brief Allocate the cell indices arrays used for the zoom region.
 *
 * @param s The #space.
 */
void zoom_allocate_cells(struct space *s) {

  /* Allocate the indices of local zoom cells */
  if (swift_memalign("local_zoom_cells_top",
                     (void **)&s->zoom_props->local_zoom_cells_top,
                     SWIFT_STRUCT_ALIGNMENT,
                     s->zoom_props->nr_zoom_cells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level zoom cells.");
  bzero(s->zoom_props->local_zoom_cells_top,
        s->zoom_props->nr_zoom_cells * sizeof(int));

  /* Allocate the indices of local bkg cells */
  if (swift_memalign("local_bkg_cells_top",
                     (void **)&s->zoom_props->local_bkg_cells_top,
                     SWIFT_STRUCT_ALIGNMENT,
                     s->zoom_props->nr_bkg_cells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->local_bkg_cells_top,
        s->zoom_props->nr_bkg_cells * sizeof(int));

  /* Allocate the indices of local buffer cells */
  if (s->zoom_props->with_buffer_cells) {
    if (swift_memalign("local_buffer_cells_top",
                       (void **)&s->zoom_props->local_buffer_cells_top,
                       SWIFT_STRUCT_ALIGNMENT,
                       s->zoom_props->nr_buffer_cells * sizeof(int)) != 0)
      error("Failed to allocate indices of local top-level background cells.");
    bzero(s->zoom_props->local_buffer_cells_top,
          s->zoom_props->nr_buffer_cells * sizeof(int));
  }

  /* Allocate the indices of local zoom cells with particles */
  if (swift_memalign(
          "local_zoom_cells_with_particles_top",
          (void **)&s->zoom_props->local_zoom_cells_with_particles_top,
          SWIFT_STRUCT_ALIGNMENT,
          s->zoom_props->nr_zoom_cells * sizeof(int)) != 0)
    error(
        "Failed to allocate indices of local top-level zoom cells with "
        "particles.");
  bzero(s->zoom_props->local_zoom_cells_with_particles_top,
        s->zoom_props->nr_zoom_cells * sizeof(int));

  /* Allocate the indices of local bkg cells with particles */
  if (swift_memalign(
          "local_bkg_cells_with_particles_top",
          (void **)&s->zoom_props->local_bkg_cells_with_particles_top,
          SWIFT_STRUCT_ALIGNMENT,
          s->zoom_props->nr_bkg_cells * sizeof(int)) != 0)
    error(
        "Failed to allocate indices of local top-level background cells with "
        "particles.");
  bzero(s->zoom_props->local_bkg_cells_with_particles_top,
        s->zoom_props->nr_bkg_cells * sizeof(int));

  /* Allocate the indices of local buffer cells with particles */
  if (s->zoom_props->with_buffer_cells) {
    if (swift_memalign(
            "local_buffer_cells_with_particles_top",
            (void **)&s->zoom_props->local_buffer_cells_with_particles_top,
            SWIFT_STRUCT_ALIGNMENT,
            s->zoom_props->nr_buffer_cells * sizeof(int)) != 0)
      error(
          "Failed to allocate indices of local top-level buffer cells with "
          "particles.");
    bzero(s->zoom_props->local_buffer_cells_with_particles_top,
          s->zoom_props->nr_buffer_cells * sizeof(int));
  }
}
