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

/* Config header */
#include <config.h>

/* Standard includes */
#include <float.h>

/* Local includes */
#include "cell.h"
#include "engine.h"
#include "multipole.h"
#include "space.h"
#include "zoom.h"

/**
 * @brief Find which background or buffer cells contain the zoom region.
 *
 * A void cell is a low resolution cell above the zoom region (or part
 * of it).
 *
 * A void cell is always in the cell grid directly above the zoom cells, i.e. if
 * there are buffer cells, the void cells are the buffer cells, if there are no
 * buffer cells, the void cells are the background cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void zoom_find_void_cells(struct space *s, const int verbose) {

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;

  /* Zero the void cell count. */
  zoom_props->nr_void_cells = 0;

  /* Get the offset and the number of cells we're dealing with. */
  int offset = zoom_props->bkg_cell_offset;
  int ncells = zoom_props->nr_bkg_cells;

  /* Allocate the indices of void cells */
  if (swift_memalign(
          "void_cell_indices", (void **)&s->zoom_props->void_cell_indices,
          SWIFT_STRUCT_ALIGNMENT, zoom_props->nr_bkg_cells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->void_cell_indices,
        zoom_props->nr_bkg_cells * sizeof(int));

  /* Loop over the background cells and find cells containing
   * the void region. */
  for (int cid = offset; cid < offset + ncells; cid++) {

    /* Get the cell */
    struct cell *c = &cells[cid];

    /* Label this cell if it contains the zoom region. */
    if (zoom_cell_overlaps_zoom_region(c, s)) {
      c->subtype = cell_subtype_void;
      zoom_props->void_cell_indices[zoom_props->nr_void_cells++] = cid;
    }
  }

  if (verbose) message("Found %d void cells", zoom_props->nr_void_cells);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check the void cells are in the right place. */
  for (int i = 0; i < zoom_props->nr_void_cells; i++) {
    const int cid = zoom_props->void_cell_indices[i];
    if (cid < offset || cid >= offset + ncells)
      error("Void cell index is out of range (cid=%d, offset=%d, ncells=%d)",
            cid, offset, ncells);
    if (zoom_cell_overlaps_zoom_region(&cells[cid], s) == 0)
      error("Void cell is not inside the zoom region (cid=%d)", cid);
  }
#endif

  /* We also need to label the buffer cells as void cells if they
   * are within the zoom region. */
  int nr_buffer_void = 0;
  if (zoom_props->with_buffer_cells) {
    offset = zoom_props->buffer_cell_offset;
    ncells = zoom_props->nr_buffer_cells;

    /* Loop over the buffer cells and find cells containing
     * the zoom region. */
    for (int cid = offset; cid < offset + ncells; cid++) {

      /* Get the cell */
      struct cell *c = &cells[cid];

      /* Label this cell if it contains the zoom region. */
      if (zoom_cell_overlaps_zoom_region(c, s)) {
        c->subtype = cell_subtype_void;
        nr_buffer_void++;
      }
    }

    if (verbose)
      message("Found %d void cells in the buffer region", nr_buffer_void);
  }
}

/**
 * @brief Find what background or buffer cells surround the zoom region.
 *
 * When interacting background cells and zoom TL cells, it helps to know
 * which background cells are within the mesh distance of the zoom region.
 * These cells are the cells which will interact via pair/grav or long_range
 * gravity tasks, and get tagged as "neighbour" cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void zoom_find_neighbouring_cells(struct space *s, const int verbose) {

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get the right cell cdim. */
  int cdim[3];
  double iwidth[3];
  if (zoom_props->with_buffer_cells) {
    cdim[0] = zoom_props->buffer_cdim[0];
    cdim[1] = zoom_props->buffer_cdim[1];
    cdim[2] = zoom_props->buffer_cdim[2];
    iwidth[0] = zoom_props->buffer_iwidth[0];
    iwidth[1] = zoom_props->buffer_iwidth[1];
    iwidth[2] = zoom_props->buffer_iwidth[2];
  } else {
    cdim[0] = s->cdim[0];
    cdim[1] = s->cdim[1];
    cdim[2] = s->cdim[2];
    iwidth[0] = s->iwidth[0];
    iwidth[1] = s->iwidth[1];
    iwidth[2] = s->iwidth[2];
  }

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;

  /* Get the right offset and the number of cells we're dealing with. */
  int offset;
  int ncells;
  if (zoom_props->with_buffer_cells) {
    offset = zoom_props->buffer_cell_offset;
    ncells = zoom_props->nr_buffer_cells;
  } else {
    offset = zoom_props->bkg_cell_offset;
    ncells = zoom_props->nr_bkg_cells;
  }

  /* At this point we can only define neighbour cells by cell properties,
   * leaving the fancy gravity distance criterion for task creation later.
   * Here we just make sure all possible neighbour cells are flagged
   * as such. */

  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
  const int delta_cells = ceil(zoom_props->neighbour_distance *
                               max3(iwidth[0], iwidth[1], iwidth[2])) +
                          1;

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Let's be verbose about this choice */
  if (verbose)
    message(
        "Looking for neighbouring background cells up to %d background "
        "top-level "
        "cells away from the zoom region (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Allocate the indices of neighbour background cells. */
  /* We don't know how many we will need at this point so have
   * to allocate assuming we will have all cells in range. */
  if (swift_memalign("neighbour_cells_top",
                     (void **)&s->zoom_props->neighbour_cells_top,
                     SWIFT_STRUCT_ALIGNMENT, ncells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->neighbour_cells_top, ncells * sizeof(int));

  /* Get a pointer to the neighbour cells index array. */
  int *neighbour_cells_top = s->zoom_props->neighbour_cells_top;

  /* Loop over background/buffer cells. We'll talk out delta_cells cells from
   * any void cells. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get this cell. */
        const size_t cid = cell_getid_offset(cdim, offset, i, j, k);
        struct cell *ci = &cells[cid];

        /* Skip non-void cells, we only want to walk out. */
        if (ci->subtype != cell_subtype_void) continue;

        /* Loop over every other cell within (Manhattan) range delta */
        for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

          /* Escape beyond range */
          if (ii < 0 || ii >= cdim[0]) continue;

          for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

            /* Escape beyond range */
            if (jj < 0 || jj >= cdim[1]) continue;

            for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

              /* Escape beyond range */
              if (kk < 0 || kk >= cdim[2]) continue;

              /* Get this cell. */
              const int cjd = cell_getid_offset(cdim, offset, ii, jj, kk);
              struct cell *cj = &cells[cjd];

              /* Ensure this neighbour isn't a void cell or an already
               * counted neighbour. */
              if (cj->subtype != cell_subtype_void &&
                  cj->subtype != cell_subtype_neighbour) {

                /* Record that we've found a neighbour. */
                cj->subtype = cell_subtype_neighbour;
                neighbour_cells_top[zoom_props->nr_neighbour_cells++] = cjd;
              }
            } /* neighbour k loop */
          } /* neighbour j loop */
        } /* neighbour i loop */
      } /* k loop */
    } /* j loop */
  } /* i loop */

  if (verbose)
    message("%i cells neighbour the zoom region",
            zoom_props->nr_neighbour_cells);
}

/**
 * @brief Run through all cells and ensure they have the correct cell type and
 * width for their position in s->cells_top.
 *
 * @param s The space.
 */
void zoom_verify_cell_type(struct space *s) {
#ifdef SWIFT_DEBUG_CHECKS

  /* Get the cells array and cell properties */
  struct cell *cells = s->cells_top;
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;
  const int buffer_offset = s->zoom_props->buffer_cell_offset;
  const double *zoom_width = s->zoom_props->width;
  const double *width = s->width;

  /* Loop over all cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Check cell type */
    if (cid < bkg_cell_offset && cells[cid].type != cell_type_zoom)
      error(
          "Cell has the wrong cell type for it's array position (cid=%d, "
          "c->type=%s, "
          "s->zoom_props->bkg_cell_offset=%d)",
          cid, cellID_names[cells[cid].type], bkg_cell_offset);
    if (cid >= bkg_cell_offset && cells[cid].type == cell_type_zoom)
      error(
          "Cell has the wrong cell type for it's array position (cid=%d, "
          "c->type=%s, "
          "s->zoom_props->bkg_cell_offset=%d)",
          cid, cellID_names[cells[cid].type], bkg_cell_offset);

    /* Check cell widths */
    for (int i = 0; i < 3; i++) {
      if (cid < bkg_cell_offset && cells[cid].width[i] != zoom_width[i])
        error(
            "Cell has the wrong cell width for it's array position (cid=%d, "
            "c->type=%s, "
            "s->zoom_props->bkg_cell_offset=%d, c->width=[%f %f %f], "
            "s->zoom_props->width=[%f %f %f])",
            cid, cellID_names[cells[cid].type], bkg_cell_offset,
            cells[cid].width[0], cells[cid].width[1], cells[cid].width[2],
            s->zoom_props->width[0], s->zoom_props->width[1],
            s->zoom_props->width[2]);
      if ((cid >= bkg_cell_offset && cid < buffer_offset) &&
          cells[cid].width[i] != width[i])
        error(
            "Cell has the wrong cell width for it's array position (cid=%d, "
            "c->type=%s, "
            "s->zoom_props->bkg_cell_offset=%d, c->width=[%f %f %f], "
            "s->zoom_props->width=[%f %f %f])",
            cid, cellID_names[cells[cid].type], bkg_cell_offset,
            cells[cid].width[0], cells[cid].width[1], cells[cid].width[2],
            s->zoom_props->width[0], s->zoom_props->width[1],
            s->zoom_props->width[2]);
    }
  }

  /* Ensure the cell and region boundaries align. */
  if (s->zoom_props->with_buffer_cells) {
    int found_bkg_bufferi_low = 0;
    int found_bkg_bufferj_low = 0;
    int found_bkg_bufferk_low = 0;
    int found_bkg_bufferi_up = 0;
    int found_bkg_bufferj_up = 0;
    int found_bkg_bufferk_up = 0;
    double tol = 0.001 * s->width[0]; /* Tolerence on edge matching */
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {

          /* Define the cell position. */
          const double pos[3] = {s->width[0] * i, s->width[1] * j,
                                 s->width[2] * k};

          /* Test the lower boundary. */
          if (fabs(pos[0] - s->zoom_props->buffer_lower_bounds[0]) < tol)
            found_bkg_bufferi_low = 1;
          if (fabs(pos[1] - s->zoom_props->buffer_lower_bounds[1]) < tol)
            found_bkg_bufferj_low = 1;
          if (fabs(pos[2] - s->zoom_props->buffer_lower_bounds[2]) < tol)
            found_bkg_bufferk_low = 1;

          /* Test the upper boundary. */
          if (fabs(pos[0] - s->zoom_props->buffer_upper_bounds[0]) < tol)
            found_bkg_bufferi_up = 1;
          if (fabs(pos[1] - s->zoom_props->buffer_upper_bounds[1]) < tol)
            found_bkg_bufferj_up = 1;
          if (fabs(pos[2] - s->zoom_props->buffer_upper_bounds[2]) < tol)
            found_bkg_bufferk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    if (!found_bkg_bufferi_low || !found_bkg_bufferj_low ||
        !found_bkg_bufferk_low || !found_bkg_bufferi_up ||
        !found_bkg_bufferj_up || !found_bkg_bufferk_up)
      error("The background cell and buffer region edges don't match!");

    /* And for the zoom cells. */
    int found_buffer_zoomi_low = 0;
    int found_buffer_zoomj_low = 0;
    int found_buffer_zoomk_low = 0;
    int found_buffer_zoomi_up = 0;
    int found_buffer_zoomj_up = 0;
    int found_buffer_zoomk_up = 0;
    tol = 0.001 * s->zoom_props->buffer_width[0];
    for (int i = 0; i < s->zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < s->zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < s->zoom_props->buffer_cdim[2]; k++) {

          /* Define the cell position. */
          const double pos[3] = {s->zoom_props->buffer_lower_bounds[0] +
                                     s->zoom_props->buffer_width[0] * i,
                                 s->zoom_props->buffer_lower_bounds[1] +
                                     s->zoom_props->buffer_width[1] * j,
                                 s->zoom_props->buffer_lower_bounds[2] +
                                     s->zoom_props->buffer_width[2] * k};

          /* Test the lower boundary. */
          if (fabs(pos[0] - s->zoom_props->region_lower_bounds[0]) < tol)
            found_buffer_zoomi_low = 1;
          if (fabs(pos[1] - s->zoom_props->region_lower_bounds[1]) < tol)
            found_buffer_zoomj_low = 1;
          if (fabs(pos[2] - s->zoom_props->region_lower_bounds[2]) < tol)
            found_buffer_zoomk_low = 1;

          /* Test the upper boundary. */
          if (fabs(pos[0] - s->zoom_props->region_upper_bounds[0]) < tol)
            found_buffer_zoomi_up = 1;
          if (fabs(pos[1] - s->zoom_props->region_upper_bounds[1]) < tol)
            found_buffer_zoomj_up = 1;
          if (fabs(pos[2] - s->zoom_props->region_upper_bounds[2]) < tol)
            found_buffer_zoomk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    if (!found_buffer_zoomi_low || !found_buffer_zoomj_low ||
        !found_buffer_zoomk_low || !found_buffer_zoomi_up ||
        !found_buffer_zoomj_up || !found_buffer_zoomk_up)
      error("The buffer cell and zoom region edges don't match!");
  } else {

    /* Check the background and zoom cells align. */
    int found_bkg_zoomi_low = 0;
    int found_bkg_zoomj_low = 0;
    int found_bkg_zoomk_low = 0;
    int found_bkg_zoomi_up = 0;
    int found_bkg_zoomj_up = 0;
    int found_bkg_zoomk_up = 0;
    const double tol = 0.001 * s->width[0];
    for (int i = 0; i < s->cdim[0]; i++) {
      for (int j = 0; j < s->cdim[1]; j++) {
        for (int k = 0; k < s->cdim[2]; k++) {

          /* Define the cell position. */
          const double pos[3] = {s->width[0] * i, s->width[1] * j,
                                 s->width[2] * k};

          /* Test the lower boundary. */
          if (fabs(pos[0] - s->zoom_props->region_lower_bounds[0]) < tol)
            found_bkg_zoomi_low = 1;
          if (fabs(pos[1] - s->zoom_props->region_lower_bounds[1]) < tol)
            found_bkg_zoomj_low = 1;
          if (fabs(pos[2] - s->zoom_props->region_lower_bounds[2]) < tol)
            found_bkg_zoomk_low = 1;

          /* Test the upper boundary. */
          if (fabs(pos[0] - s->zoom_props->region_upper_bounds[0]) < tol)
            found_bkg_zoomi_up = 1;
          if (fabs(pos[1] - s->zoom_props->region_upper_bounds[1]) < tol)
            found_bkg_zoomj_up = 1;
          if (fabs(pos[2] - s->zoom_props->region_upper_bounds[2]) < tol)
            found_bkg_zoomk_up = 1;
        }
      }
    }
    /* Did we find the boundaries?. */
    if (!found_bkg_zoomi_low || !found_bkg_zoomj_low || !found_bkg_zoomk_low ||
        !found_bkg_zoomi_up || !found_bkg_zoomj_up || !found_bkg_zoomk_up)
      error("The background cell and zoom region edges don't match!");
  }
#else
  error("zoom_verify_cell_type() called without SWIFT_DEBUG_CHECKS defined");
#endif
}

/**
 * @brief Build the TL cells, with a zoom region.
 *
 * This replaces the loop in space_regrid when running with a zoom region. It
 * constructs all zoom, background and buffer cells (if required) and sets their
 * initial values.
 *
 * Zoom cells occupy the first cells in the cells array, followed by background
 * cells and then buffer cells (if required).
 *
 * @param s The space.
 * @param ti_current The current time.
 * @param verbose Are we talking?
 */
void zoom_construct_tl_cells(struct space *s, const integertime_t ti_current,
                             int verbose) {

  const ticks tic = getticks();

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get some zoom region properties */
  const float dmin_zoom =
      min3(zoom_props->width[0], zoom_props->width[1], zoom_props->width[2]);
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_bounds[3] = {zoom_props->region_lower_bounds[0],
                                 zoom_props->region_lower_bounds[1],
                                 zoom_props->region_lower_bounds[2]};

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < zoom_props->cdim[0]; i++) {
    for (int j = 0; j < zoom_props->cdim[1]; j++) {
      for (int k = 0; k < zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(zoom_props->cdim, i, j, k);

        /* Create the zoom cell and it's multipoles */
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->width[0]) + zoom_bounds[0];
        c->loc[1] = (j * zoom_props->width[1]) + zoom_bounds[1];
        c->loc[2] = (k * zoom_props->width[2]) + zoom_bounds[2];
        c->width[0] = zoom_props->width[0];
        c->width[1] = zoom_props->width[1];
        c->width[2] = zoom_props->width[2];
        if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
        c->type = cell_type_zoom;
        c->subtype = cell_subtype_regular;
        c->dmin = dmin_zoom;
        c->h_min_allowed =
            c->dmin * 0.5 * (1. / kernel_gamma); /* Only zoom cells do hydro */
        c->h_max_allowed =
            c->dmin * (1. / kernel_gamma); /* Only zoom cells do hydro */
        c->depth = 0;
        c->split = 0;
        c->hydro.count = 0;
        c->grav.count = 0;
        c->stars.count = 0;
        c->sinks.count = 0;
        c->top = c;
        c->super = c;
        c->hydro.super = c;
        c->grav.super = c;
        c->hydro.ti_old_part = ti_current;
        c->grav.ti_old_part = ti_current;
        c->stars.ti_old_part = ti_current;
        c->sinks.ti_old_part = ti_current;
        c->black_holes.ti_old_part = ti_current;
        c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
        c->mpi.tag = -1;
        c->mpi.recv = NULL;
        c->mpi.send = NULL;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif
      }
    }
  }

  /* Attach zoom cells to their cell array (they are the first cells in
   * cells_top so nothing fancy needed here). */
  zoom_props->zoom_cells_top = s->cells_top;

  if (verbose)
    message("Set zoom cell dimensions to [ %i %i %i ].", zoom_props->cdim[0],
            zoom_props->cdim[1], zoom_props->cdim[2]);

  /* Get the minimum cell size. */
  const double dmin = min3(s->width[0], s->width[1], s->width[2]);

  /* Loop over background cells and set locations and initial values */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        const size_t cid = cell_getid_offset(s->cdim, bkg_cell_offset, i, j, k);

        /* Background top level cells. */
        struct cell *c = &s->cells_top[cid];
        c->loc[0] = i * s->width[0];
        c->loc[1] = j * s->width[1];
        c->loc[2] = k * s->width[2];
        c->width[0] = s->width[0];
        c->width[1] = s->width[1];
        c->width[2] = s->width[2];
        c->dmin = dmin;
        if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
        c->type = cell_type_bkg;
        c->subtype = cell_subtype_regular;
        c->depth = 0;
        c->split = 0;
        c->hydro.count = 0;
        c->grav.count = 0;
        c->stars.count = 0;
        c->sinks.count = 0;
        c->top = c;
        c->super = c;
        c->hydro.super = c;
        c->grav.super = c;
        c->hydro.ti_old_part = ti_current;
        c->grav.ti_old_part = ti_current;
        c->stars.ti_old_part = ti_current;
        c->sinks.ti_old_part = ti_current;
        c->black_holes.ti_old_part = ti_current;
        c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
        c->mpi.tag = -1;
        c->mpi.recv = NULL;
        c->mpi.send = NULL;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif
      }
    }
  }

  /* Attach background cells to their cell array. */
  zoom_props->bkg_cells_top = &s->cells_top[bkg_cell_offset];

  if (verbose)
    message("Set background cell dimensions to [ %i %i %i ].", s->cdim[0],
            s->cdim[1], s->cdim[2]);

  /* If we have a buffer region create buffer cells. */
  if (zoom_props->with_buffer_cells) {

    /* Get relevant information. */
    const float dmin_buffer =
        min3(zoom_props->buffer_width[0], zoom_props->buffer_width[1],
             zoom_props->buffer_width[2]);
    const int buffer_offset = zoom_props->buffer_cell_offset;
    const double buffer_bounds[3] = {zoom_props->buffer_lower_bounds[0],
                                     zoom_props->buffer_lower_bounds[1],
                                     zoom_props->buffer_lower_bounds[2]};

    /* Loop over buffer cells and set locations and initial values */
    for (int i = 0; i < zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < zoom_props->buffer_cdim[2]; k++) {
          const size_t cid = cell_getid_offset(zoom_props->buffer_cdim,
                                               buffer_offset, i, j, k);

          /* Buffer top level cells. */

          struct cell *c = &s->cells_top[cid];
          c->loc[0] = (i * zoom_props->buffer_width[0]) + buffer_bounds[0];
          c->loc[1] = (j * zoom_props->buffer_width[1]) + buffer_bounds[1];
          c->loc[2] = (k * zoom_props->buffer_width[2]) + buffer_bounds[2];
          c->width[0] = zoom_props->buffer_width[0];
          c->width[1] = zoom_props->buffer_width[1];
          c->width[2] = zoom_props->buffer_width[2];
          c->dmin = dmin_buffer;
          if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
          c->type = cell_type_buffer;
          c->subtype = cell_subtype_regular;
          c->depth = 0;
          c->split = 0;
          c->hydro.count = 0;
          c->grav.count = 0;
          c->stars.count = 0;
          c->sinks.count = 0;
          c->top = c;
          c->super = c;
          c->hydro.super = c;
          c->grav.super = c;
          c->hydro.ti_old_part = ti_current;
          c->grav.ti_old_part = ti_current;
          c->stars.ti_old_part = ti_current;
          c->sinks.ti_old_part = ti_current;
          c->black_holes.ti_old_part = ti_current;
          c->grav.ti_old_multipole = ti_current;
#ifdef WITH_MPI
          c->mpi.tag = -1;
          c->mpi.recv = NULL;
          c->mpi.send = NULL;
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
          cell_assign_top_level_cell_index(c, s);
#endif
        }
      }
    }

    /* Attach background cells to their cell array. */
    zoom_props->buffer_cells_top = &s->cells_top[buffer_offset];

    if (verbose)
      message("Set buffer cell dimensions to [ %i %i %i ].",
              zoom_props->buffer_cdim[0], zoom_props->buffer_cdim[1],
              zoom_props->buffer_cdim[2]);
  }

  /* Now find and label what cells contain the zoom region. */
  zoom_find_void_cells(s, verbose);

  /* Find neighbours of the zoom region (cells within the distance past which
   * the gravity mesh takes over). */
  zoom_find_neighbouring_cells(s, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Lets check all the cells are in the right place with the correct widths */
  zoom_verify_cell_type(s);
#endif

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Collect the timestep from the zoom region to the void level.
 *
 * This will recurse to the void leaves and grab the timestep from the zoom
 * cells, populating the void cell tree as it goes.
 *
 * @param c The #cell.
 */
static void zoom_void_timestep_collect_recursive(struct cell *c) {

  /* Define the timestep info we'll be collecting. */
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_grav_end_min = max_nr_timesteps, ti_grav_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we have the right kind of cell. */
  if (c->subtype != cell_subtype_void) {
    error(
        "Trying to update timesteps on cell which isn't a void cell! "
        "(c->type=%s, c->subtype=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);
  }
#endif

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL) {

      /* Recurse, if we're still in the void cell tree. */
      if (cp->subtype == cell_subtype_void) {
        zoom_void_timestep_collect_recursive(cp);
      }

      /* And update */
      ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
      ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
      ti_rt_end_min = min(cp->rt.ti_rt_end_min, ti_rt_end_min);
      ti_rt_beg_max = max(cp->rt.ti_rt_beg_max, ti_rt_beg_max);
      ti_grav_end_min = min(ti_grav_end_min, cp->grav.ti_end_min);
      ti_grav_beg_max = max(ti_grav_beg_max, cp->grav.ti_beg_max);
      ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);
      ti_stars_beg_max = max(ti_stars_beg_max, cp->stars.ti_beg_max);
      ti_black_holes_end_min =
          min(ti_black_holes_end_min, cp->black_holes.ti_end_min);
      ti_black_holes_beg_max =
          max(ti_black_holes_beg_max, cp->black_holes.ti_beg_max);
      ti_sinks_end_min = min(ti_sinks_end_min, cp->sinks.ti_end_min);
      ti_sinks_beg_max = max(ti_sinks_beg_max, cp->sinks.ti_beg_max);
    }
  }

  /* Store the collected values in the cell. */
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->rt.ti_rt_end_min = ti_rt_end_min;
  c->rt.ti_rt_beg_max = ti_rt_beg_max;
  c->grav.ti_end_min = ti_grav_end_min;
  c->grav.ti_beg_max = ti_grav_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
}

/**
 * @brief Collect the timestep from the zoom region to the void level.
 *
 * @param e The engine.
 */
void zoom_void_timestep_collect(struct engine *e) {

  ticks tic = getticks();

  /* Get the space and void cells. */
  struct space *s = e->s;
  struct zoom_region_properties *zoom_props = s->zoom_props;
  const int nr_void_cells = zoom_props->nr_void_cells;
  const int *void_cell_indices = zoom_props->void_cell_indices;

  /* Loop over all void cells and collect the timesteps. */
  for (int i = 0; i < nr_void_cells; i++) {
    struct cell *c = &s->cells_top[void_cell_indices[i]];
    zoom_void_timestep_collect_recursive(c);
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
