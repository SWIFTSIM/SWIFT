/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will J. Roper (w.roper@sussex.ac.uk)
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

/* Config */
#include <config.h>

/* Local includes */
#include "cell.h"
#include "engine.h"
#include "space.h"
#include "threadpool.h"
#include "zoom.h"

/**
 * @brief Link the top level cells in zoom grid to a void parent cell.
 *
 * The leaves of the void cell hierarchy are top level cells in the nested grid.
 * This function sets the progeny of the highest res void cell to be the top
 * level cell it contains (either zoom or buffer cells).
 *
 * NOTE: The void cells with top level progeny are not treated as split cells
 * since they are linked into the top level "progeny". We don't want to
 * accidentally treat them as split cells and recurse from void cells straight
 * through to the zoom cells unless explictly desired.
 *
 * @param s The space.
 * @param c The void cell progeny to link
 */
static void zoom_link_void_zoom_leaves(struct space *s, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we have the right kind of cell. */
  if (c->subtype != cell_subtype_void) {
    error(
        "Trying to split cell which isn't a void cell! (c->type=%s, "
        "c->subtype=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);
  }

  /* Check that the widths are right. */
  if (fabs((c->width[0] / 2) - s->zoom_props->width[0]) >
      (0.001 * s->zoom_props->width[0]))
    error(
        "The width of the zoom cell is not half the width of the void "
        "cell were about to link to! (c->width[0]=%f, "
        "s->zoom_props->width[0]=%f)",
        c->width[0] / 2, s->zoom_props->width[0]);

#endif

  /* We need to ensure this bottom level isn't treated like a
   * normal split cell since it's linked into top level "progeny". */
  c->split = 0;

  /* Loop over the 8 progeny cells which are now the nested top level cells. */
  for (int k = 0; k < 8; k++) {

    /* Establish the location of the fake progeny cell. */
    double loc[3] = {c->loc[0] + (c->width[0] / 4),
                     c->loc[1] + (c->width[1] / 4),
                     c->loc[2] + (c->width[2] / 4)};
    if (k & 4) loc[0] += c->width[0] / 2;
    if (k & 2) loc[1] += c->width[1] / 2;
    if (k & 1) loc[2] += c->width[2] / 2;

    /* Which cell are we in? */
    int cid = cell_getid_from_pos(s, loc[0], loc[1], loc[2]);

    /* Get the zoom cell. */
    struct cell *zoom_cell = &s->cells_top[cid];

    // /* Is this zoom progeny empty? Empty means no particles (and if running
    //  * with self-gravity no multipole mass). When running with self-gravity
    //  * over MPI we can have 0 particles but we will have multipole mass from
    //  * the other ranks so we need to keep the cell in the tree. */
    // int empty = 0;
    // if (!s->with_self_gravity && zoom_cell->grav.count == 0 &&
    //     zoom_cell->hydro.count == 0 && zoom_cell->stars.count == 0 &&
    //     zoom_cell->sinks.count == 0 && zoom_cell->black_holes.count == 0) {
    //   empty = 1;
    // } else if (s->with_self_gravity &&
    //            zoom_cell->grav.multipole->m_pole.M_000 == 0.f) {
    //   empty = 1;
    // }
    //
    // /* If this top level cell is empty, don't link it in. */
    // if (empty) {
    //   c->progeny[k] = NULL;
    //   continue;
    // }

    /* Link this nested cell into the void cell hierarchy. */
    c->progeny[k] = zoom_cell;

    /* Flag this void cell "progeny" as the cell's void cell parent. */
    zoom_cell->void_parent = c;
  }
}

/**
 * @brief Link the top level cells in buffer grid to a void parent cell.
 *
 * The leaves of the void cell hierarchy are top level cells in the nested grid.
 * This function sets the progeny of the highest res void cell to be the top
 * level cell it contains (either zoom or buffer cells).
 *
 * @param s The space.
 * @param c The void cell progeny to link
 */
void zoom_link_void_buffer_leaves(struct space *s, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we have the right kind of cell. */
  if (c->subtype != cell_subtype_void) {
    error(
        "Trying to split cell which isn't a void cell! (c->type=%s, "
        "c->subtype=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);
  }

  /* Check that the widths are right. */
  if (fabs((c->width[0] / 2) - s->zoom_props->buffer_width[0]) >
      (0.001 * s->zoom_props->buffer_width[0]))
    error(
        "The width of the buffer cell is not half the width of the void "
        "cell were about to link to! (c->width[0]=%f, "
        "s->zoom_props->buffer_width[0]=%f)",
        c->width[0] / 2, s->zoom_props->buffer_width[0]);

#endif

  /* Loop over the 8 progeny cells which are now the nested top level cells. */
  for (int k = 0; k < 8; k++) {

    /* Establish the location of the fake progeny cell. */
    double loc[3] = {c->loc[0] + (c->width[0] / 4),
                     c->loc[1] + (c->width[1] / 4),
                     c->loc[2] + (c->width[2] / 4)};
    if (k & 4) loc[0] += c->width[0] / 2;
    if (k & 2) loc[1] += c->width[1] / 2;
    if (k & 1) loc[2] += c->width[2] / 2;

    /* Which cell are we in? */
    int cid = cell_getid_below_bkg(s->zoom_props->buffer_cdim,
                                   s->zoom_props->buffer_lower_bounds, loc[0],
                                   loc[1], loc[2], s->zoom_props->buffer_iwidth,
                                   s->zoom_props->buffer_cell_offset);

    /* Get the zoom cell. */
    struct cell *buffer_cell = &s->cells_top[cid];

    /* Link this nested cell into the void cell hierarchy. */
    c->progeny[k] = buffer_cell;

    /* Flag this void cell "progeny" as the cell's void cell parent. */
    buffer_cell->void_parent = c;

    /* Set the parent of the buffer cell to be the void cell. */
    buffer_cell->parent = c;

    /* If we're in a void cell continue the void tree depth. */
    if (buffer_cell->subtype == cell_subtype_void) {
      buffer_cell->depth = c->depth + 1;
    }
  }
}

/**
 * @brief Recursively split a cell.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param tpid The thread id.
 */
void zoom_void_split_recursive(struct space *s, struct cell *c,
                               const short int tpid) {

  const int depth = c->depth;
  int maxdepth = 0;

  /* Initialise the timestep information we need to collect. */
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_rt_min_step_size = max_nr_timesteps;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;

  /* Set the top level void cell tpid. Doing it here ensures top level void
   * cells have the same tpid as their progeny. */
  if (depth == 0) c->tpid = tpid;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we haven't found a void cell with particles. */
  if (c->subtype == cell_subtype_void && c->grav.count > 0)
    error(
        "Trying to split a Void with particles! "
        "(c->type=%s, c->subtype=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);
#endif

  /* If the depth is too large, we have a problem and should stop. */
  if (depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting the void cells, aborting. "
        "This is most likely due to having the zoom region too deep within the "
        "background cells. Try using buffer cells or increasing "
        "region_buffer_cell_ratio (depth=%d)",
        space_cell_maxdepth, depth);
  }

  /* Construct or attach the progeny ready to populate the multipoles (if
   * doing gravity). If we are above one of the nested top level cell grids
   * we will attach those existing cells rather than grab new ones. */

  /* If we're above the zoom level we need to link in the zoom cells. */
  if (c->depth == s->zoom_props->zoom_cell_depth - 1) {
    zoom_link_void_zoom_leaves(s, c);
  }

  /* If we're above the buffer level we need to link in the buffer cells. */
  else if (c->depth == s->zoom_props->buffer_cell_depth - 1) {
    zoom_link_void_buffer_leaves(s, c);
  }

  /* Otherwise, we're in the void tree and need to construct new progeny. */
  else {
    space_construct_progeny(s, c, tpid);
  }

  for (int k = 0; k < 8; k++) {

    /* Get the progenitor */
    struct cell *cp = c->progeny[k];

    /* Skip NULL progeny. */
    if (cp == NULL) continue;

    /* If the progeny is a void cell, we need to recurse. */
    if (cp->subtype == cell_subtype_void) {

      /* Recurse */
      zoom_void_split_recursive(s, cp, tpid);

      /* Increase the depth */
      maxdepth = max(maxdepth, cp->maxdepth);
    }

    /* If the progeny is a non-local non-void cell we're done. */
    else if (cp->nodeID != engine_rank) {
      continue;
    }

    /* Update the timestep information. */
    ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
    ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
    ti_rt_end_min = min(ti_rt_end_min, cp->rt.ti_rt_end_min);
    ti_rt_beg_max = max(ti_rt_beg_max, cp->rt.ti_rt_beg_max);
    ti_rt_min_step_size = min(ti_rt_min_step_size, cp->rt.ti_rt_min_step_size);
    ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
    ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);
    ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);
    ti_stars_beg_max = max(ti_stars_beg_max, cp->stars.ti_beg_max);
    ti_sinks_end_min = min(ti_sinks_end_min, cp->sinks.ti_end_min);
    ti_sinks_beg_max = max(ti_sinks_beg_max, cp->sinks.ti_beg_max);
    ti_black_holes_end_min =
        min(ti_black_holes_end_min, cp->black_holes.ti_end_min);
    ti_black_holes_beg_max =
        max(ti_black_holes_beg_max, cp->black_holes.ti_beg_max);
  }

  /* Update the properties of the void cell. */
  c->maxdepth = maxdepth;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->rt.ti_rt_end_min = ti_rt_end_min;
  c->rt.ti_rt_beg_max = ti_rt_beg_max;
  c->rt.ti_rt_min_step_size = ti_rt_min_step_size;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;

  /* Deal with the multipole */
  if (s->with_self_gravity) {
    space_populate_multipole(c);
  }
}

/**
 * @brief Split particles between cells of the void cell hierarchy.
 *
 * This has to be done after proxy exchange to ensure we have all the
 * foreign multipoles, so is separated from all other splitting done in the
 * usual space_split.
 *
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void zoom_void_space_split(struct space *s, int verbose) {

#ifdef SWIFT_DEBUG_CHECKS
  /* We should never get here when not running with a zoom region. */
  if (!s->with_zoom_region) {
    error("zoom_void_space_split called when running without a zoom region.");
  }
#endif

  const ticks tic = getticks();

  /* Unpack some useful information. */
  struct cell *cells_top = s->cells_top;
  int *void_cell_indices = s->zoom_props->void_cell_indices;
  int nr_void_cells = s->zoom_props->nr_void_cells;

  /* Create the void cell trees and populate their multipoles. This is only
   * a handful of cells so no threadpool. */

  /* Loop over the void cells */
  for (int ind = 0; ind < nr_void_cells; ind++) {
    struct cell *c = &cells_top[void_cell_indices[ind]];
    zoom_void_split_recursive(s, c, /*tpid*/ 0);

    /* Reset the contains zoom cells flag, we're about to recompute it. */
    c->contains_zoom_cells = 0;
  }

  /* Flag which void cells are useful (i.e. contain zoom cells) by looping
   * over zoom cells. */
  int nr_useful_voids = 0;
  for (int cid = 0; cid < s->zoom_props->nr_zoom_cells; cid++) {
    struct cell *zoom_cell = &cells_top[cid];

    /* Get the void parent. */
    struct cell *void_parent = zoom_cell->void_parent;

    /* Skip if we don't have a void parent (empty zoom cell). */
    if (void_parent == NULL) continue;

    /* Skip empty zoom cells. */
    if (zoom_cell->grav.count == 0) continue;

    /* Get the top level void cell. */
    struct cell *top = void_parent->top;

    /* Flag that this void cell contains zoom cells. */
    if (!top->contains_zoom_cells) {
      nr_useful_voids++;
      top->contains_zoom_cells = 1;
    }
  }

  if (verbose) {
    message("Found %d void cells with zoom cells out of %d total void cells.",
            nr_useful_voids, nr_void_cells);
    message("Void cell tree and multipole construction took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

#ifdef SWIFT_DEBUG_CHECKS

  /* Ensure all buffer cells are linked into the tree. */
  int notlinked = 0;
  if (s->zoom_props->with_buffer_cells) {
    for (int k = s->zoom_props->buffer_cell_offset;
         k < s->zoom_props->buffer_cell_offset + s->zoom_props->nr_buffer_cells;
         k++) {
      if (cells_top[k].void_parent == NULL) notlinked++;
    }
    if (notlinked > 0)
      error("%d buffer cells are not linked into a void cell tree!", notlinked);
  }

  /* Ensure all zoom cells are linked into the tree. */
  notlinked = 0;
  for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {
    if (cells_top[k].void_parent == NULL && cells_top[k].grav.count > 0)
      notlinked++;
  }
  if (notlinked > 0)
    error("%d zoom cells are not linked into a void cell tree!", notlinked);

  /* When running with gravity we need to make sure our multipole counts make
   * sense between cell levels. */
  if (s->with_self_gravity) {
    /* Collect the number of particles in the void multipoles. */
    int nr_gparts_in_void = 0;
    for (int i = 0; i < nr_void_cells; i++) {
      nr_gparts_in_void +=
          s->multipoles_top[s->zoom_props->void_cell_indices[i]]
              .m_pole.num_gpart;
    }

    /* Collect the number of particles in the buffer multipoles. */
    int nr_gparts = 0;
    if (s->zoom_props->with_buffer_cells) {
      for (int k = s->zoom_props->buffer_cell_offset;
           k <
           s->zoom_props->buffer_cell_offset + s->zoom_props->nr_buffer_cells;
           k++) {
        nr_gparts += s->multipoles_top[k].m_pole.num_gpart;
      }
    } else {
      /* With no buffer cells we instead need to count up the zoom cell
       * gpart counts. */
      /* Collect the number of particles in the zoom multipoles. */
      for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {
        nr_gparts += s->multipoles_top[k].m_pole.num_gpart;
      }
    }

    /* Check the number of gparts is consistent. */
    if (s->zoom_props->with_buffer_cells && nr_gparts_in_void != nr_gparts) {
      error(
          "Number of gparts is inconsistent between buffer cells and "
          "void multipole (nr_gparts_in_void=%d, nr_gparts=%d)",
          nr_gparts_in_void, nr_gparts);
    }

    /* Check the number of particles in the void cells. */
    else if (!s->zoom_props->with_buffer_cells &&
             nr_gparts_in_void != nr_gparts) {
      error(
          "Number of gparts is inconsistent between zoom cells and "
          "void multipole (nr_gparts_in_void=%d, nr_gparts=%d)",
          nr_gparts_in_void, nr_gparts);
    }

    /* Ensure the zoom count agrees between the buffer void cells and zoom cells
     * (if we have buffer cells). */
    if (s->zoom_props->with_buffer_cells) {
      int zoom_count = 0;
      int buffer_count = 0;
      for (int k = 0; k < s->nr_cells; k++) {
        struct cell *c = &s->cells_top[k];
        if (c->type == cell_type_zoom) {
          zoom_count += c->grav.multipole->m_pole.num_gpart;
        } else if (c->type == cell_type_buffer &&
                   c->subtype == cell_subtype_void) {
          buffer_count += c->grav.multipole->m_pole.num_gpart;
        }
      }

      if (zoom_count != buffer_count) {
        error(
            "Buffer and zoom multipoles disagree (buffer_count=%d, "
            "zoom_count=%d)",
            buffer_count, zoom_count);
      }
    }

#ifdef WITH_MPI
    /* Ensure all void cells agree on their time zone across ranks and the
     * timesteps have been initialised correctly. */
    integertime_t *global_ti_beg_max =
        (integertime_t *)malloc(nr_void_cells * sizeof(integertime_t));
    integertime_t *global_ti_end_min =
        (integertime_t *)malloc(nr_void_cells * sizeof(integertime_t));
    for (int ind = 0; ind < nr_void_cells; ind++) {
      struct cell *c = &cells_top[void_cell_indices[ind]];
      integertime_t local_ti_beg_max = c->grav.ti_beg_max;
      integertime_t local_ti_end_min = c->grav.ti_end_min;
      MPI_Allreduce(&local_ti_beg_max, &global_ti_beg_max[ind], 1,
                    MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&local_ti_end_min, &global_ti_end_min[ind], 1,
                    MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
    }

    /* Count how many void cells have bad timesteps. */
    int bad_min_timesteps = 0;
    int bad_max_timesteps = 0;
    for (int ind = 0; ind < nr_void_cells; ind++) {
      struct cell *c = &cells_top[void_cell_indices[ind]];
      /* We only care if the void cell matters (i.e. contains zoom cells). */
      if (!c->contains_zoom_cells) continue;
      if (c->grav.ti_end_min == -1) bad_min_timesteps++;
      if (c->grav.ti_beg_max == -1) bad_max_timesteps++;
    }

    /* If we have any bad timesteps we have a problem. */
    if (bad_min_timesteps > 0) {
      message(
          "%d void cells have incorrect minimum end timesteps after MPI "
          "reduction!",
          bad_min_timesteps);
    }
    if (bad_max_timesteps > 0) {
      message(
          "%d void cells have incorrect maximum begin timesteps after MPI "
          "reduction!",
          bad_max_timesteps);
    }
#endif
  }
#endif
}
