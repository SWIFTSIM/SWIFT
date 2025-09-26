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
 * @brief Check whether we need a regrid based on the zoom region motion.
 *
 * If the zoom region has moved more than a certain fraction of its size
 * in any dimension then we need to regrid. This is to ensure that the
 * zoom region remains well centred on the particle distribution.
 *
 * @param s The #space.
 * @return 1 if a zoom regrid is needed, 0 otherwise.
 */
static int zoom_need_regrid_motion(const struct space *s) {

  const double dx[3] = {s->zoom_props->zoom_shift[0],
                        s->zoom_props->zoom_shift[1],
                        s->zoom_props->zoom_shift[2]};
  const double max_shift = s->zoom_props->max_com_dx;
  if (dx[0] > max_shift * s->zoom_props->dim[0] ||
      dx[1] > max_shift * s->zoom_props->dim[1] ||
      dx[2] > max_shift * s->zoom_props->dim[2]) {
    message(
        "Zoom region shift exceeds %.2f%% of the zoom region in at least one "
        "dimension (shift=(%.3e,%.3e,%.3e), zoom_dim=(%.3e,%.3e,%.3e)).",
        max_shift * 100.0, dx[0], dx[1], dx[2], s->zoom_props->dim[0],
        s->zoom_props->dim[1], s->zoom_props->dim[2]);
    return 1;
  }

  return 0;
}

/**
 * @brief Check whether we need a regrid based on the particle extent.
 *
 * If the particle distribution exceeds a certain fraction of the zoom region
 * size in any dimension then we need to regrid. This is to ensure that the
 * zoom region remains well centred on the particle distribution.
 *
 * @param s The #space.
 * @return 1 if a zoom regrid is needed, 0 otherwise.
 */
static int zoom_need_regrid_extent(const struct space *s) {

  /* Derive the maximum allowed particle extent from the user specified
   * target padding fraction, if the particle distribution is more than
   * this fraction then we are no longer doing what the user asked for. */
  const double max_part_dim_frac = 1.0 / s->zoom_props->user_region_pad_factor;

  /* Get the maximum allowed extent of the particle distribution. */
  const double max_dim[3] = {s->zoom_props->dim[0] * max_part_dim_frac,
                             s->zoom_props->dim[1] * max_part_dim_frac,
                             s->zoom_props->dim[2] * max_part_dim_frac};

  /* If the particle distribution is more than 90% of the zoom region width
   * (based on the zoom cells themselves) in any dimension we need to regrid. */
  if (s->zoom_props->part_dim[0] > max_dim[0] ||
      s->zoom_props->part_dim[1] > max_dim[1] ||
      s->zoom_props->part_dim[2] > max_dim[2]) {
    message(
        "Particle distribution exceeds %.2f%% of the zoom region "
        "in at least one dimension (part_dim=(%.3e,%.3e,%.3e), "
        "zoom_dim=(%.3e,%.3e,%.3e)).",
        max_part_dim_frac * 100.0, s->zoom_props->part_dim[0],
        s->zoom_props->part_dim[1], s->zoom_props->part_dim[2],
        s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
    return 1;
  }
  return 0;
}

/**
 * @brief Check whether we need a regrid based on hmax requiring larger cells.
 *
 * If the current hmax requires larger cells in the zoom region than we
 * currently have then we need to regrid.
 *
 * @param s The #space.
 * @param new_cdim The new top-level cell dimensions (based on current hmax).
 * @return 1 if a zoom regrid is needed, 0 otherwise.
 */
static int zoom_need_regrid_hmax(const struct space *s, const int new_cdim[3]) {
  return (new_cdim[0] < s->zoom_props->cdim[0] ||
          new_cdim[1] < s->zoom_props->cdim[1] ||
          new_cdim[2] < s->zoom_props->cdim[2]);
}

/**
 * @brief Check whether we need a regrid based on the zoom region.
 *
 * This function is called in space_regrid in space_regrid.c and provides
 * the zoom specific regrid check.
 *
 * Unlike a uniform box, in a zoom simulation there are three reasons we might
 * need to regrid:
 * 1) The high resolution particle distribution has moved too far from the
 *   centre of the zoom region.
 * 2) The high resolution particle distribution has expanded too far and is
 *  likely to be no longer well contained within the zoom region soon.
 * 3) [Same as a uniform box] The maximum smoothing length (hmax) has increased
 * such that the current zoom cells are too small to properly resolve the
 * hydrodynamics.
 *
 * @param s The #space.
 * @param new_cdim The new top-level cell dimensions (based on current hmax).
 *
 * @return 1 if a zoom regrid is needed, 0 otherwise.
 */
int zoom_need_regrid(const struct space *s, const int new_cdim[3]) {

  /* Have we exceeded the allowed shift of the zoom region? */
  if (zoom_need_regrid_motion(s)) {
    return 1;
  }

  /* Has the particle distribution exceeded the allowed extent? */
  if (zoom_need_regrid_extent(s)) {
    return 1;
  }

  /* Has hmax increased such that we need larger zoom cells? */
  if (zoom_need_regrid_hmax(s, new_cdim)) {
    return 1;
  }
  return 0;
}

/**
 * @brief Find an acceptable geometry given the required zoom cdim.
 *
 * This function has multiple use cases:
 * 1) We are regridding due to the zoom region motion. In this case the
 *  particle distribution itself needs shifting and the cells need to be
 *  rebuilt around it.
 * 2) We are regridding due to the zoom region extent. In this case the
 *  particle distribution itself needs shifting and the cells need to be
 *  rebuilt around expanding the width of the zoom region.
 * 3) We are regridding due to hmax requiring larger cells in the zoom region.
 *  In this case we need to find a new geometry that can accommodate the new
 *  cdim (note that the particles will be shifted in this case too).
 *
 * If we are doing the first case we simply call zoom_region_init to shift and
 * recalculate the geometry and we're done. The second and third cases are more
 * complex and require changes to the cell dimensions.
 *
 * For the second and third cases we'll first try to decrease the background
 * cdim a reasonable amount since this will carry less of a performance penalty
 * than doubling the size of zoom cells. This is also most likely to produce a
 * valid set up in all but the most extreme cases.
 *
 * @param s The #space.
 * @param new_cdim The new top-level cell dimensions (based on current hmax).
 */
void zoom_regrid_find_acceptable_geometry(struct space *s,
                                          const int new_cdim[3]) {

  /* If we do not need a regrid just exit. */
  if (!zoom_need_regrid(s, new_cdim)) {
    return;
  }

  /* If we are regridding due to motion then just try to recentre the zoom
   * region. */
  if (zoom_need_regrid_motion(s) && !zoom_need_regrid_extent(s) &&
      !zoom_need_regrid_hmax(s, new_cdim)) {

    /* Recalculate the zoom region geometry. (silently) */
    zoom_region_init(s, /*regridding=*/1, /*verbose=*/0);
    return;
  }

  /* Otherwise, we need to also adjust input cell geometry. */

  /* Loop until we've found an acceptable geometry. */
  int old_bkg_cdim = s->cdim[0];
  while (zoom_need_regrid(s, new_cdim)) {

    /* First try decreasing the background cdim to a minimum of 50% its
     * current value. */
    while (zoom_need_regrid(s, new_cdim) && s->cdim[0] > 0.5 * old_bkg_cdim) {

      /* Decrement the background cdim. */
      s->zoom_props->bkg_cdim[0]--;
      s->zoom_props->bkg_cdim[1]--;
      s->zoom_props->bkg_cdim[2]--;

      /* Recalculate the zoom region geometry. (silently) */
      zoom_region_init(s, /*regridding=*/1, /*verbose=*/0);
    }

    /* If this worked we can stop here. */
    if (!zoom_need_regrid(s, new_cdim)) {
      break;
    }

    /* If we aren't doing a hmax based regrid we have now failed, adjusting
     * the zoom region depth won't help us. */
    if (!zoom_need_regrid_hmax(s, new_cdim)) {
      error(
          "Failed to find an acceptable zoom region before reaching 50%% of "
          "the initial background cell size (i.e. not a hmax based regrid). "
          "Try to decrease the initial background cell size.");
    }

    /* Reset the background cdim to its original value, we'll try decreasing
     * the zoom region depth instead, and then loop again. */
    s->zoom_props->bkg_cdim[0] = old_bkg_cdim;
    s->zoom_props->bkg_cdim[1] = old_bkg_cdim;
    s->zoom_props->bkg_cdim[2] = old_bkg_cdim;

    /* If adjusting the background didn't work then decrease the zoom region
     * depth by one (doubling the zoom cell size). */
    s->zoom_props->zoom_cell_depth--;

    /* Ensure we don't go below depth 1. */
    if (s->zoom_props->zoom_cell_depth < 1) {
      error(
          "Failed to find an acceptable zoom region before reaching depth=0 "
          "(i.e. not a zoom). Try to decrease the initial background cell "
          "size.");
    }

    /* Recalculate the zoom region geometry. (silently) */
    zoom_region_init(s, /*regridding=*/1, /*verbose=*/0);
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
