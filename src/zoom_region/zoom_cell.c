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
  const int j = (y - bounds[2]) * iwidth[1];
  const int k = (z - bounds[4]) * iwidth[2];

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
  int cell_id;

  /* Lets get some properties of the zoom region. */
  const struct zoom_region_properties *zoom_props = s->zoom_props;
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_region_bounds[6] = {
      zoom_props->region_bounds[0], zoom_props->region_bounds[1],
      zoom_props->region_bounds[2], zoom_props->region_bounds[3],
      zoom_props->region_bounds[4], zoom_props->region_bounds[5]};
  const int buffer_cell_offset = zoom_props->buffer_cell_offset;
  const double buffer_bounds[6] = {
      zoom_props->buffer_bounds[0], zoom_props->buffer_bounds[1],
      zoom_props->buffer_bounds[2], zoom_props->buffer_bounds[3],
      zoom_props->buffer_bounds[4], zoom_props->buffer_bounds[5]};

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
    cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_region_bounds, x,
                                     y, z, s->zoom_props->iwidth,
                                     /*offset*/ 0);

  }

  /* If this is an empty cell we are in the buffer cells.
   * Otherwise, It's a legitimate background cell, and we'll return it. */
  else if (s->cells_top[cell_id].subtype == empty) {

    /* Which buffer TL cell are we in? */
    cell_id = cell_getid_with_bounds(s->zoom_props->buffer_cdim, buffer_bounds,
                                     x, y, z, s->zoom_props->buffer_iwidth,
                                     buffer_cell_offset);

    /* Here we need to check if this is the void buffer cell.
     * Otherwise, It's a legitimate buffer cell, and we'll return it. */
    if (s->cells_top[cell_id].subtype == void_cell) {

      /* Which zoom TL cell are we in? */
      cell_id = cell_getid_with_bounds(s->zoom_props->cdim, zoom_region_bounds,
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
  const struct zoom_region_properties *zoom_props = s->zoom_props;

  /* Get some zoom region properties */
  const float dmin_zoom =
      min3(zoom_props->width[0], zoom_props->width[1], zoom_props->width[2]);
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_region_bounds[6] = {
      zoom_props->region_bounds[0], zoom_props->region_bounds[1],
      zoom_props->region_bounds[2], zoom_props->region_bounds[3],
      zoom_props->region_bounds[4], zoom_props->region_bounds[5]};

  struct cell *restrict c;

  /* Loop over zoom cells and set locations and initial values */
  for (int i = 0; i < zoom_props->cdim[0]; i++) {
    for (int j = 0; j < zoom_props->cdim[1]; j++) {
      for (int k = 0; k < zoom_props->cdim[2]; k++) {
        const size_t cid = cell_getid(zoom_props->cdim, i, j, k);

        /* Create the zoom cell and it's multipoles */
        c = &s->cells_top[cid];
        c->loc[0] = (i * zoom_props->width[0]) + zoom_region_bounds[0];
        c->loc[1] = (j * zoom_props->width[1]) + zoom_region_bounds[2];
        c->loc[2] = (k * zoom_props->width[2]) + zoom_region_bounds[4];
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
        c->subtype = none;
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
        const size_t cid = cell_getid(s->cdim, i, j, k) + bkg_cell_offset;

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
        c->subtype = none;
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

        /* Assign the cell type. */
        if (zoom_props->with_buffer_cells && cell_inside_buffer_region(c, s)) {
          c->subtype = empty;
          zoom_props->nr_empty_cells++;
        }
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
    const double buffer_bounds[6] = {
        zoom_props->buffer_bounds[0], zoom_props->buffer_bounds[1],
        zoom_props->buffer_bounds[2], zoom_props->buffer_bounds[3],
        zoom_props->buffer_bounds[4], zoom_props->buffer_bounds[5]};

    /* Loop over buffer cells and set locations and initial values */
    for (int i = 0; i < zoom_props->buffer_cdim[0]; i++) {
      for (int j = 0; j < zoom_props->buffer_cdim[1]; j++) {
        for (int k = 0; k < zoom_props->buffer_cdim[2]; k++) {
          const size_t cid =
              cell_getid(zoom_props->buffer_cdim, i, j, k) + buffer_offset;

          /* Buffer top level cells. */
          c = &s->cells_top[cid];
          c->loc[0] = (i * zoom_props->buffer_width[0]) + buffer_bounds[0];
          c->loc[1] = (j * zoom_props->buffer_width[1]) + buffer_bounds[2];
          c->loc[2] = (k * zoom_props->buffer_width[2]) + buffer_bounds[4];
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
          c->subtype = none;
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

  /* Now find what cells contain the zoom region. */
  find_void_cells(s, verbose);

  /* Now find what cells neighbour the zoom region. */
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

/**
 * @brief Find what TL cells contain the zoom region.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_void_cells(struct space *s, const int verbose) {
#ifdef WITH_ZOOM_REGION

  /* Get the right cell cdim. */
  int cdim[3];
  if (zoom_props->with_buffer_cells) {
    cdim[0] = s->zoom_props->buffer_cdim[0];
    cdim[1] = s->zoom_props->buffer_cdim[1];
    cdim[2] = s->zoom_props->buffer_cdim[2];
  } else {
    cdim[0] = s->cdim[0];
    cdim[1] = s->cdim[1];
    cdim[2] = s->cdim[2];
  }

  /* Get the cell pointers. */
  struct cell *cells = s->cells_top;

  /* Get the right offset and the number of cells we're dealing with. */
  int offset;
  int ncells;
  if (s->zoom_props->with_buffer_cells) {
    offset = s->zoom_props->buffer_cell_offset;
    ncells = s->zoom_props->nr_buffer_cells;
  } else {
    offset = s->zoom_props->bkg_cell_offset;
    ncells = s->zoom_props->nr_bkg_cells;
  }

  /* Allocate the indices of void cells */
  /* TODO: We can know ahead of time how many void cells there are. We don't
   * need to allocate this many.  */
  int void_count = 0;
  if (swift_memalign("void_cells_top", (void **)&s->zoom_props->void_cells_top,
                     SWIFT_STRUCT_ALIGNMENT, ncells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->void_cells_top, ncells * sizeof(int));

  /* Loop over natural cells and find cells containing the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {
        const size_t cid = cell_getid(cdim, i, j, k) + offset;

        /* Label this background cell. */
        if (cell_contains_zoom_region(&cells[cid], s)) {
          cells[cid].subtype = void_cell;
          s->zoom_props->void_cells_top[void_count++] = cid;
        }
      }
    }
  }

  /* Store the number of neighbour cells */
  s->zoom_props->nr_void_cells = void_count;

  if (void_count == 0)
    error("No void cells were found! (nr_buffer_cells=%d)",
          s->zoom_props->nr_buffer_cells);

  if (verbose) message("%i cells contain the zoom region", void_count);

#endif
}

/**
 * @brief Find what TL cells surround the zoom region.
 *
 * When interacting "natural" TL cells and "zoom" TL cells, it helps to know
 * what natural TL cells surround the zoom region. These cells then get tagged
 * as "neighbour" cells.
 *
 * @param s The space.
 * @param verbose Are we talking?
 */
void find_neighbouring_cells(struct space *s,
                             struct gravity_props *gravity_properties,
                             const int verbose) {
#ifdef WITH_ZOOM_REGION

  /* Get the right cell cdim. */
  int cdim[3];
  double iwidth[3];
  if (s->zoom_props->with_buffer_cells) {
    cdim[0] = s->zoom_props->buffer_cdim[0];
    cdim[1] = s->zoom_props->buffer_cdim[1];
    cdim[2] = s->zoom_props->buffer_cdim[2];
    iwidth[0] = s->zoom_props->buffer_iwidth[0];
    iwidth[1] = s->zoom_props->buffer_iwidth[1];
    iwidth[2] = s->zoom_props->buffer_iwidth[2];
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
  if (s->zoom_props->with_buffer_cells) {
    offset = s->zoom_props->buffer_cell_offset;
    ncells = s->zoom_props->nr_buffer_cells;
  } else {
    offset = s->zoom_props->bkg_cell_offset;
    ncells = s->zoom_props->nr_bkg_cells;
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

  /* Allocate the indices of neighbour background cells */
  if (swift_memalign("neighbour_cells_top",
                     (void **)&s->zoom_props->neighbour_cells_top,
                     SWIFT_STRUCT_ALIGNMENT, ncells * sizeof(int)) != 0)
    error("Failed to allocate indices of local top-level background cells.");
  bzero(s->zoom_props->neighbour_cells_top, ncells * sizeof(int));

  int neighbour_count = 0;

  /* Loop over natural cells and find cells neighbouring the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get this cell. */
        const size_t cid = cell_getid(cdim, i, j, k) + offset;
        struct cell *ci = &cells[cid];

        /* Skip non-void cells. */
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
              const int cjd = cell_getid(cdim, ii, jj, kk) + offset;

              /* Ensure this neighbour isn't a void cell or an already
               * counted neighbour. */
              if (cells[cjd].subtype != void_cell &&
                  cells[cjd].subtype != neighbour) {

                /* Record that we've found a neighbour. */
                cells[cjd].subtype = neighbour;
                s->zoom_props->neighbour_cells_top[neighbour_count++] = cjd;
              }
            } /* neighbour k loop */
          }   /* neighbour j loop */
        }     /* neighbour i loop */
      }       /* natural k loop */
    }         /* natural j loop */
  }           /* natural i loop */

  /* Store the number of neighbour cells */
  s->zoom_props->nr_neighbour_cells = neighbour_count;

  if (verbose)
    message("%i cells neighbouring the zoom region", neighbour_count);
#endif
}
