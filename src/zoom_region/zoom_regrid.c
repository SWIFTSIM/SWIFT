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
  if (fabs(dx[0]) > max_shift * s->zoom_props->dim[0] ||
      fabs(dx[1]) > max_shift * s->zoom_props->dim[1] ||
      fabs(dx[2]) > max_shift * s->zoom_props->dim[2]) {
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
    message("Regridding due to zoom region motion.");
    return 1;
  }

  /* Has the particle distribution exceeded the allowed extent? */
  if (zoom_need_regrid_extent(s)) {
    message("Regridding due to zoom region extent.");
    return 1;
  }

  /* Has hmax increased such that we need larger zoom cells? */
  if (zoom_need_regrid_hmax(s, new_cdim)) {
    message("Regridding due to hmax requiring larger zoom cells.");
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

    /* Recalculate the zoom region geometry. */
    zoom_region_init(s, /*regridding=*/1, /*verbose=*/1);
    return;
  }

  /* Otherwise, we need to also adjust input cell geometry. */

  /* Loop until we've found an acceptable geometry. */
  int old_bkg_cdim = s->cdim[0];
  while (zoom_need_regrid(s, new_cdim)) {

    /* First try decreasing the background cdim to a minimum of 50% its
     * current value. */
    while (zoom_need_regrid(s, new_cdim) && s->cdim[0] > 0.5 * old_bkg_cdim) {

      message(
          "Adjusting background cell size to find acceptable zoom geometry ("
          "bkg_cdim=(%d,%d,%d)).",
          s->zoom_props->bkg_cdim[0] - 2, s->zoom_props->bkg_cdim[1] - 2,
          s->zoom_props->bkg_cdim[2] - 2);

      /* Preserve grid parity: changing it moves the centred zoom region from a
       * background-cell centre to a boundary, or vice versa, discontinuously
       * changing the void geometry. */
      s->zoom_props->bkg_cdim[0] -= 2;
      s->zoom_props->bkg_cdim[1] -= 2;
      s->zoom_props->bkg_cdim[2] -= 2;

      /* Recalculate the zoom region geometry. */
      zoom_region_init(s, /*regridding=*/1, /*verbose=*/1);
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

    message(
        "Adjusting zoom cell size to find acceptable zoom geometry "
        "(depth=%d).",
        s->zoom_props->zoom_cell_depth - 1);

    /* If adjusting the background didn't work then decrease the zoom region
     * depth by one (doubling the zoom cell size). */
    const int old_zoom_cell_depth = s->zoom_props->zoom_cell_depth;
    s->zoom_props->zoom_cell_depth--;

    /* Ensure we don't go below depth 1. */
    if (s->zoom_props->zoom_cell_depth < 1) {
      error(
          "Failed to find an acceptable zoom region before reaching depth=0 "
          "(i.e. not a zoom). Try to decrease the initial background cell "
          "size.");
    }

    /* If the geometry was originally driven by a requested cdim, update that
     * target too, otherwise zoom_region_init() will derive the old depth again.
     */
    if (s->zoom_props->zoom_cell_min_cdim > 0) {
      const int parent_cdim =
          s->zoom_props->cdim[0] / (1 << old_zoom_cell_depth);
      s->zoom_props->zoom_cell_min_cdim =
          parent_cdim * (1 << s->zoom_props->zoom_cell_depth);
    }

    /* Recalculate the zoom region geometry. */
    zoom_region_init(s, /*regridding=*/1, /*verbose=*/1);
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
 * includes counts from each top-level cell grid (zoom and bkg).
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
    swift_free("void_cell_indices", s->zoom_props->void_cell_indices);
    swift_free("neighbour_cells_top", s->zoom_props->neighbour_cells_top);

    /* Find an acceptable geometry given the required zoom cdim. */
    zoom_regrid_find_acceptable_geometry(s, zoom_cdim);
  }

  /* Count the number of top level cells. */
  s->tot_cells = s->nr_cells =
      (s->cdim[0] * s->cdim[1] * s->cdim[2]) +
      (s->zoom_props->cdim[0] * s->zoom_props->cdim[1] *
       s->zoom_props->cdim[2]);
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
}

/**
 * @brief Put existing "extra" particles at the centres of their zoom cells.
 *
 * Zoom regridding can move cells while leaving the number of spares unchanged.
 * In that case space_allocate_extras will not run its placement loops. When
 * the count changes, those loops already place the old and new spares.
 *
 * @param s The #space with newly constructed zoom cells and old spares.
 */
void zoom_space_shift_extras(struct space *s) {

  const size_t nr_zoom_cells = s->zoom_props->nr_zoom_cells;
  const size_t nr_gparts = s->nr_gparts;
  const size_t nr_parts = s->nr_parts;
  const size_t nr_sinks = s->nr_sinks;
  const size_t nr_sparts = s->nr_sparts;
  const size_t nr_bparts = s->nr_bparts;
  const size_t nr_extra_gparts = s->nr_extra_gparts;
  const size_t nr_extra_parts = s->nr_extra_parts;
  const size_t nr_extra_sinks = s->nr_extra_sinks;
  const size_t nr_extra_sparts = s->nr_extra_sparts;
  const size_t nr_extra_bparts = s->nr_extra_bparts;
  const size_t extra_gparts_per_cell = space_extra_gparts;
  const size_t extra_parts_per_cell = space_extra_parts;
  const size_t extra_sinks_per_cell = space_extra_sinks;
  const size_t extra_sparts_per_cell = space_extra_sparts;
  const size_t extra_bparts_per_cell = space_extra_bparts;

  /* Nothing to do if there are no extras. */
  if (nr_extra_parts == 0 && nr_extra_gparts == 0 && nr_extra_sparts == 0 &&
      nr_extra_bparts == 0 && nr_extra_sinks == 0)
    return;

  /* Count local zoom cells to see whether the number of extras changed. */
  size_t nr_local_cells = 0;
  for (size_t i = 0; i < nr_zoom_cells; ++i) {
    if (s->cells_top[i].nodeID == s->e->nodeID) ++nr_local_cells;
  }

  /* Shift extra gparts to the centre of their new cells. */
  if (nr_local_cells > 0 && nr_extra_gparts > 0 &&
      nr_extra_gparts == nr_local_cells * extra_gparts_per_cell) {
    size_t next = 0;

    /* Loop over the new zoom cells owned by this rank. */
    for (size_t i = 0; i < nr_zoom_cells; ++i) {
      const struct cell *c = &s->cells_top[i];
      if (c->nodeID != s->e->nodeID) continue;

      /* Skip particles which aren't extras. */
      for (; next < nr_gparts; ++next) {
        if (s->gparts[next].time_bin == time_bin_not_created) break;
      }

      /* Each cell must still have its full block of extras. */
      if (nr_gparts - next < extra_gparts_per_cell) {
        error("Missing extra gparts");
      }

      /* Move this cell's extras to its new centre. */
      for (size_t j = 0; j < extra_gparts_per_cell; ++j) {
        struct gpart *gp = &s->gparts[next + j];
        if (gp->time_bin != time_bin_not_created) {
          error("Extra gparts are not contiguous");
        }
        for (int k = 0; k < 3; ++k) {
          gp->x[k] = c->loc[k] + c->width[k] / 2.;
        }
      }

      /* Continue with the next block of extras. */
      next += extra_gparts_per_cell;
    }
  }

  /* Shift extra parts to the centre of their new cells. */
  if (nr_local_cells > 0 && nr_extra_parts > 0 &&
      nr_extra_parts == nr_local_cells * extra_parts_per_cell) {
    size_t next = 0;

    /* Loop over the new zoom cells owned by this rank. */
    for (size_t i = 0; i < nr_zoom_cells; ++i) {
      const struct cell *c = &s->cells_top[i];
      if (c->nodeID != s->e->nodeID) continue;

      /* Skip particles which aren't extras. */
      for (; next < nr_parts; ++next) {
        if (s->parts[next].time_bin == time_bin_not_created) break;
      }

      /* Each cell must still have its full block of extras. */
      if (nr_parts - next < extra_parts_per_cell) {
        error("Missing extra parts");
      }

      /* Move this cell's extras to its new centre. */
      for (size_t j = 0; j < extra_parts_per_cell; ++j) {
        struct part *p = &s->parts[next + j];
        if (p->time_bin != time_bin_not_created) {
          error("Extra parts are not contiguous");
        }
        for (int k = 0; k < 3; ++k) {
          p->x[k] = c->loc[k] + c->width[k] / 2.;
        }
      }

      /* Continue with the next block of extras. */
      next += extra_parts_per_cell;
    }
  }

  /* Shift extra sinks to the centre of their new cells. */
  if (nr_local_cells > 0 && nr_extra_sinks > 0 &&
      nr_extra_sinks == nr_local_cells * extra_sinks_per_cell) {
    size_t next = 0;

    /* Loop over the new zoom cells owned by this rank. */
    for (size_t i = 0; i < nr_zoom_cells; ++i) {
      const struct cell *c = &s->cells_top[i];
      if (c->nodeID != s->e->nodeID) continue;

      /* Skip particles which aren't extras. */
      for (; next < nr_sinks; ++next) {
        if (s->sinks[next].time_bin == time_bin_not_created) break;
      }

      /* Each cell must still have its full block of extras. */
      if (nr_sinks - next < extra_sinks_per_cell) {
        error("Missing extra sinks");
      }

      /* Move this cell's extras to its new centre. */
      for (size_t j = 0; j < extra_sinks_per_cell; ++j) {
        struct sink *sink = &s->sinks[next + j];
        if (sink->time_bin != time_bin_not_created) {
          error("Extra sinks are not contiguous");
        }
        for (int k = 0; k < 3; ++k) {
          sink->x[k] = c->loc[k] + c->width[k] / 2.;
        }
      }

      /* Continue with the next block of extras. */
      next += extra_sinks_per_cell;
    }
  }

  /* Shift extra sparts to the centre of their new cells. */
  if (nr_local_cells > 0 && nr_extra_sparts > 0 &&
      nr_extra_sparts == nr_local_cells * extra_sparts_per_cell) {
    size_t next = 0;

    /* Loop over the new zoom cells owned by this rank. */
    for (size_t i = 0; i < nr_zoom_cells; ++i) {
      const struct cell *c = &s->cells_top[i];
      if (c->nodeID != s->e->nodeID) continue;

      /* Skip particles which aren't extras. */
      for (; next < nr_sparts; ++next) {
        if (s->sparts[next].time_bin == time_bin_not_created) break;
      }

      /* Each cell must still have its full block of extras. */
      if (nr_sparts - next < extra_sparts_per_cell) {
        error("Missing extra sparts");
      }

      /* Move this cell's extras to its new centre. */
      for (size_t j = 0; j < extra_sparts_per_cell; ++j) {
        struct spart *sp = &s->sparts[next + j];
        if (sp->time_bin != time_bin_not_created) {
          error("Extra sparts are not contiguous");
        }
        for (int k = 0; k < 3; ++k) {
          sp->x[k] = c->loc[k] + c->width[k] / 2.;
        }
      }

      /* Continue with the next block of extras. */
      next += extra_sparts_per_cell;
    }
  }

  /* Shift extra bparts to the centre of their new cells. */
  if (nr_local_cells > 0 && nr_extra_bparts > 0 &&
      nr_extra_bparts == nr_local_cells * extra_bparts_per_cell) {
    size_t next = 0;

    /* Loop over the new zoom cells owned by this rank. */
    for (size_t i = 0; i < nr_zoom_cells; ++i) {
      const struct cell *c = &s->cells_top[i];
      if (c->nodeID != s->e->nodeID) continue;

      /* Skip particles which aren't extras. */
      for (; next < nr_bparts; ++next) {
        if (s->bparts[next].time_bin == time_bin_not_created) break;
      }

      /* Each cell must still have its full block of extras. */
      if (nr_bparts - next < extra_bparts_per_cell) {
        error("Missing extra bparts");
      }

      /* Move this cell's extras to its new centre. */
      for (size_t j = 0; j < extra_bparts_per_cell; ++j) {
        struct bpart *bp = &s->bparts[next + j];
        if (bp->time_bin != time_bin_not_created) {
          error("Extra bparts are not contiguous");
        }
        for (int k = 0; k < 3; ++k) {
          bp->x[k] = c->loc[k] + c->width[k] / 2.;
        }
      }

      /* Continue with the next block of extras. */
      next += extra_bparts_per_cell;
    }
  }
}
