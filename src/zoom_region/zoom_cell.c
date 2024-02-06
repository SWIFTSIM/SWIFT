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

/**
 * @brief Find what background or buffer cells contain the zoom region.
 *
 * A void cell is a cell that contains the zoom region. If there is no buffer
 * cell grid, then the void cells are the background cells. If there is a
 * buffer cell grid, then the void cells are the buffer cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_void_cells(struct space *s, const int verbose) {

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

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

  /* Work out how many void cells we should have. */
  int target_void_count;
  if (zoom_props->with_buffer_cells) {
    /* With buffer cells this is simple since it's requested by the user. */
    target_void_count = zoom_props->region_buffer_ratio *
                        zoom_props->region_buffer_ratio *
                        zoom_props->region_buffer_ratio;
  } else {
    /* With only background cells we need to work this out from the
     * background cell width */
    const int void_cdim = zoom_props->dim[0] * s->iwidth[0];
    target_void_count = void_cdim * void_cdim * void_cdim;
  }

  /* Allocate the indices of void cells */
  if (swift_memalign("void_cells_top", (void **)&s->zoom_props->void_cells_top,
                     SWIFT_STRUCT_ALIGNMENT,
                     target_void_count * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->void_cells_top, target_void_count * sizeof(int));

  /* Loop over the background/buffer cells and find cells containing
   * the zoom region. */
  for (int cid = offset; cid < ncells; cid++) {

    /* Get the cell */
    struct cell *c = &cells[cid];

    /* Label this cell if it contains the zoom region. */
    if (cell_inside_zoom_region(c, s)) {
      c->subtype = void_cell;
      zoom_props->void_cells_top[zoom_props->nr_void_cells++] = cid;
    }
  }

  if (target_void_count != zoom_props->nr_void_cells)
    error(
        "Not all void cells were found and labelled! (target_void_count=%d, "
        "found_void_count=%d)",
        target_void_count, zoom_props->nr_void_cells);

  if (verbose)
    message("%i void cells contain the zoom region", zoom_props->nr_void_cells);
}

/**
 * @brief Find what background cells contain the buffer region.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_empty_cells(struct space *s, const int verbose) {

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;

  /* Get the offset and the number of cells we're dealing with. */
  int offset = zoom_props->bkg_cell_offset;
  int ncells = zoom_props->nr_bkg_cells;

  /* Loop over natural cells and find cells containing the zoom region. */
  for (int cid = offset; cid < offset + ncells; cid++) {

    /* Get this cell. */
    struct cell *c = &cells[cid];

    /* Assign the cell type. */
    if (cell_inside_buffer_region(c, s)) {
      c->subtype = empty;
      zoom_props->nr_empty_cells++;
    }
  }

  if (zoom_props->nr_empty_cells == 0)
    error("No empty cells were found! (nr_bkg_cells=%d)", ncells);

  if (verbose)
    message("%i background cells contain the buffer region",
            zoom_props->nr_empty_cells);
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
void find_neighbouring_cells(struct space *s,
                             struct gravity_props *gravity_properties,
                             const int verbose) {

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

  /* Get gravity mesh distance. */
  const double max_distance =
      gravity_properties->r_s * gravity_properties->r_cut_max_ratio;

  /* At this point we can only define neighbour cells by cell properties,
   * leaving the fancy gravity distance criterion for task creation later.
   * Here we just make sure all possible neighbour cells are flagged
   * as such. */

  /* Maximal distance any interaction can take place
   * before the mesh kicks in, rounded up to the next integer */
  const int delta_cells =
      ceil(max_distance * max3(iwidth[0], iwidth[1], iwidth[2])) + 1;

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
        "Looking for neighbouring natural cells up to %d natural top-level "
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
        if (ci->subtype != void_cell) continue;

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
              if (cj->subtype != void_cell && cj->subtype != neighbour) {

                /* Record that we've found a neighbour. */
                cj->subtype = neighbour;
                neighbour_cells_top[zoom_props->nr_neighbour_cells++] = cjd;
              }
            } /* neighbour k loop */
          }   /* neighbour j loop */
        }     /* neighbour i loop */
      }       /* k loop */
    }         /* j loop */
  }           /* i loop */

  if (verbose)
    message("%i cells neighbour the zoom region",
            zoom_props->nr_neighbour_cells);
}

/**
 * @brief Build the TL cells, with a zoom region.
 *
 * This replaces the loop in space_regrid when running with a zoom region.
 *
 * Construct an additional set of TL "zoom" cells embedded within the TL cell
 * structure with the dimensions of each cell structure being the same (with
 * differing widths).
 *
 * Therefore the new TL cell structure is 2*cdim**3, with the "natural" TL cells
 * occupying the first half of the TL cell list, and the "zoom" TL cells
 * ocupying the second half.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void zoom_construct_tl_cells(struct space *s, const int *cdim, const float dmin,
                             const integertime_t ti_current,
                             struct gravity_props *gravity_properties,
                             int verbose) {

  /* Get the zoom properties */
  struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get some zoom region properties */
  const float dmin_zoom =
      min3(zoom_props->width[0], zoom_props->width[1], zoom_props->width[2]);
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_bounds[3] = {zoom_props->region_lower_bounds[0],
                                 zoom_props->region_lower_bounds[1],
                                 zoom_props->region_lower_bounds[2]};

  struct cell *restrict c;

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < zoom_props->cdim[0]; i++) {
    for (int j = 0; j < zoom_props->cdim[1]; j++) {
      for (int k = 0; k < zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(zoom_props->cdim, i, j, k);

        /* Create the zoom cell and it's multipoles */
        c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->width[0]) + zoom_bounds[0];
        c->loc[1] = (j * zoom_props->width[1]) + zoom_bounds[1];
        c->loc[2] = (k * zoom_props->width[2]) + zoom_bounds[2];
        c->width[0] = zoom_props->width[0];
        c->width[1] = zoom_props->width[1];
        c->width[2] = zoom_props->width[2];
        const size_t parent_cid =
            cell_getid(cdim, (int)(c->loc[0] + (c->width[0] / 2)) / s->width[0],
                       (int)(c->loc[1] + (c->width[1] / 2)) / s->width[1],
                       (int)(c->loc[2] + (c->width[2] / 2)) / s->width[2]) +
            bkg_cell_offset;
        if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
        c->type = zoom;
        c->subtype = none_sub;
        c->dmin = dmin_zoom;
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
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
        c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif
      }
    }
  }

  /* Loop over natural cells and set locations and initial values */
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        const size_t cid = cell_getid_offset(s->cdim, bkg_cell_offset, i, j, k);

        /* Natural top level cells. */
        c = &s->cells_top[cid];
        c->loc[0] = i * s->width[0];
        c->loc[1] = j * s->width[1];
        c->loc[2] = k * s->width[2];
        c->width[0] = s->width[0];
        c->width[1] = s->width[1];
        c->width[2] = s->width[2];
        c->dmin = dmin;
        if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
        c->type = bkg;
        c->subtype = none_sub;
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
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
        c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
        cell_assign_top_level_cell_index(c, s);
#endif
      }
    }
  }

  /* If we have a buffer region create buffer cells. */
  if (zoom_props->with_buffer_cells) {

    if (verbose)
      message("%i background cells are in the buffer region",
              zoom_props->nr_empty_cells);

    /* Get relevant information. */
    const float dmin_buffer =
        min3(zoom_props->buffer_width[0], zoom_props->buffer_width[1],
             zoom_props->buffer_width[2]);
    const int buffer_offset = zoom_props->buffer_cell_offset;
    const double buffer_bounds[6] = {zoom_props->buffer_lower_bounds[0],
                                     zoom_props->buffer_lower_bounds[1],
                                     zoom_props->buffer_lower_bounds[2]};

    /* Loop over buffer cells and set locations and initial values */
    for (int i = 0; i < zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < zoom_props->buffer_cdim[2]; k++) {
          const size_t cid =
              cell_getid(zoom_props->buffer_cdim, i, j, k) + buffer_offset;

          /* Buffer top level cells. */
          c = &s->cells_top[cid];
          c->loc[0] = (i * zoom_props->buffer_width[0]) + buffer_bounds[0];
          c->loc[1] = (j * zoom_props->buffer_width[1]) + buffer_bounds[1];
          c->loc[2] = (k * zoom_props->buffer_width[2]) + buffer_bounds[2];
          c->width[0] = zoom_props->buffer_width[0];
          c->width[1] = zoom_props->buffer_width[1];
          c->width[2] = zoom_props->buffer_width[2];
          c->dmin = dmin_buffer;
          const size_t parent_cid =
              cell_getid(
                  zoom_props->buffer_cdim,
                  (int)((c->loc[0] + (c->width[0] / 2)) - buffer_bounds[0]) *
                      zoom_props->buffer_iwidth[0],
                  (int)((c->loc[1] + (c->width[1] / 2)) - buffer_bounds[2]) *
                      zoom_props->buffer_iwidth[1],
                  (int)((c->loc[2] + (c->width[2] / 2)) - buffer_bounds[4]) *
                      zoom_props->buffer_iwidth[2]) +
              buffer_offset;
          if (s->with_self_gravity) c->grav.multipole = &s->multipoles_top[cid];
          c->type = buffer;
          c->subtype = none_sub;
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
#if (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
          c->nr_vertex_edges = 0;
#endif
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
          cell_assign_top_level_cell_index(c, s);
#endif
        }
      }
    }
  }

  /* Now find and label what cells contain the zoom region. */
  find_void_cells(s, verbose);

  /* If we have a buffer region then find and label the empty cells. */
  if (zoom_props->with_buffer_cells) {
    find_empty_cells(s, verbose);
  }

  /* Find neighbours of the zoom region (cells within the distance past which
   * the gravity mesh takes over). */
  find_neighbouring_cells(s, gravity_properties, verbose);

  /* NOTE: The below code relies on functions that will be implemented down
   * the line in later partition focuses PRs but keeping this here as a
   * reminder */
  /* #if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS)) */

  /*   /\* Find the number of edges we will need for the domain decomp. *\/ */
  /*   find_vertex_edges(s, verbose); */

  /*   if (verbose) */
  /*     message("%i vertex 'edges' found in total", zoom_props->nr_edges); */

  /* #endif */

#ifdef SWIFT_DEBUG_CHECKS
  /* Lets check all the cells are in the right place with the correct widths */
  debug_cell_type(s);
#endif
}
