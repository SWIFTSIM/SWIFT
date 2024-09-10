/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2024 Will J. Roper (w.roper@sussex.ac.uk).
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

/* Standard includes. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "multipole.h"
#include "proxy.h"
#include "space.h"
#include "zoom_region/zoom.h"

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for the background cells.
 *
 * This mapper only considers bkg->bkg interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 *
 * This will create pair tasks between void cells and background cells (if
 * running without buffer cells). These pair tasks will be split into smaller
 * tasks during task splitting to make the most of any possible void mm
 * interactions above the zoom level in the cell tree.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_bkg_cells(void *map_data,
                                                     int num_elements,
                                                     void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->zoom_props->bkg_cells_top;

  /* We always use the mesh if the volume is periodic. */
  const int use_mesh = s->periodic;

  /* Unlike buffer or zoom cells, background cells are periodic at the box
   * boundaries if the space is periodic. */
  const int periodic = s->periodic;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[0].width[0], s->max_softening,
      s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  /* NOTE: The 2 in the max below may not be necessary but does insure some
   * safety buffer. */
  const int delta = max((int)(sqrt(3) * distance / cells[0].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (periodic) {
    if (delta >= cdim[0] / 2) {
      if (cdim[0] % 2 == 0) {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2 - 1;
      } else {
        delta_m = cdim[0] / 2;
        delta_p = cdim[0] / 2;
      }
    }
  } else {
    if (delta > cdim[0]) {
      delta_m = cdim[0];
      delta_p = cdim[0];
    }
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Create a self task, and loop over neighbouring cells making pair tasks
     * where appropriate. */
    engine_gravity_make_task_loop(e, (size_t)(map_data) + ind, cdim, cells,
                                  periodic, use_mesh, delta_m, delta_p);
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for the buffer cells.
 *
 * This mapper only considers buffer->buffer interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 *
 * This will create pair tasks between void cells and buffer cells.
 * These pair tasks will be split into smallerrtasks during task splitting to
 * make the most of any possible void mm interactions above the zoom level in
 * the cell tree.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_buffer_cells(void *map_data,
                                                        int num_elements,
                                                        void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  const int cdim[3] = {s->zoom_props->buffer_cdim[0],
                       s->zoom_props->buffer_cdim[1],
                       s->zoom_props->buffer_cdim[2]};
  struct cell *cells = s->zoom_props->buffer_cells_top;

  /* We always use the mesh if the volume is periodic. */
  const int use_mesh = s->periodic;

  /* The buffer region is never periodic at it's boundaries. */
  const int periodic = 0;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[0].width[0], s->max_softening,
      s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  /* NOTE: The 2 in the max below may not be necessary but does insure some
   * safety buffer. */
  const int delta = max((int)(sqrt(3) * distance / cells[0].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Create a self task, and loop over neighbouring cells making pair tasks
     * where appropriate. */
    engine_gravity_make_task_loop(e, (size_t)(map_data) + ind, cdim, cells,
                                  periodic, use_mesh, delta_m, delta_p);
  }
}

/**
 * @brief Constructs the top-level self gravity tasks zoom cells.
 *
 * All local top-level zoom cells get a self task if they contain particles.
 *
 * zoom pair tasks are created during task splitting as the pair tasks
 * including void cells are split into smaller tasks. This is done so we can
 * make the most of any possible void mm interactions above the zoom level in
 * the cell tree.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  struct space *s = e->s;
  struct cell *cells = s->zoom_props->zoom_cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *c = &cells[cid];

    /* Skip cells without gravity particles. */
    if (c->grav.count == 0) return;

    /* If the cell is local build a self-interaction */
    if (c->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, c,
                        NULL);
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between natural level cells
 * and zoom level cells.
 *
 * This mapper only consider buffer->bkg interactions. It will only be called
 * if buffer cells are used. This is also the only place we create tasks
 * between different cell types since these will not be covered by splitting
 * a cell hierarchy (zoom cells lie within the void cell tree, buffer cells
 * do not lie within a background cell tree).
 *
 * - All top-cells get a self task.
 * - All pairs of differing sized cells within range according to
 *   the multipole acceptance criterion get a pair task.
 *
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_buffer_bkg(void *map_data,
                                                      int num_elements,
                                                      void *extra_data) {

  /* Useful local information */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;

  /* Some info about the zoom domain */
  const int buffer_offset = s->zoom_props->buffer_cell_offset;
  const int bkg_offset = s->zoom_props->bkg_cell_offset;

  /* We always use the mesh if the volume is periodic. */
  const int use_mesh = s->periodic;

  /* We need to account for periodic boundary conditions when background cells
   * are involved. */
  const int periodic = s->periodic;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind + buffer_offset;

    /* Get the cell */
    struct cell *ci = &cells[cid];

#ifdef SWIFT_DEBUG_CHECKS
    if (ci->type != cell_type_buffer) {
      error("Buffer cell is not a buffer cell. (ci->type = %s)",
            cellID_names[ci->type]);
    }
#endif

    /* Skip cells without gravity particles unless they're voids. */
    if (ci->grav.count == 0 && ci->subtype != cell_subtype_void) continue;

    /* Loop over every neighbouring background cells */
    for (int cjd = bkg_offset; cjd < buffer_offset; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

#ifdef SWIFT_DEBUG_CHECKS
      if (cj->type != cell_type_bkg) {
        error("Background cell is not a background cell. (cj->type = %s)",
              cellID_names[cj->type]);
      }
#endif

      /* Avoid empty cells and completely foreign pairs */
      if (cj->grav.count == 0 || (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      /* Do we need a pair interaction for these cells? */
      if (engine_gravity_need_cell_pair_task(e, ci, cj, periodic, use_mesh)) {

        /* Ok, we need to add a direct pair calculation */
        engine_make_pair_gravity_task(e, sched, ci, cj, nodeID, cid, cjd);
      }
    }
  }
}

/**
 * @brief Construct the hierarchical tasks for the void cell tree recursively.
 *
 * This will construct:
 * - The init for preparing void cell multipoles.
 * - The init implicit task for the void cells.
 * - The long-range gravity task for the void cells.
 * - The down-pass gravity task for the void cells.
 * - The down-pass implicit task for the void cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
static void zoom_engine_make_hierarchical_void_tasks_recursive(struct engine *e,
                                                               struct cell *c) {

  struct scheduler *s = &e->sched;
  const int is_self_gravity = (e->policy & engine_policy_self_gravity);

  /* At the top level we have a few different tasks to make. */
  if (c->subtype == cell_subtype_void && c->grav.super == c) {

    if (is_self_gravity) {

      /* Initialisation of the multipoles */
      c->grav.init = scheduler_addtask(s, task_type_init_grav,
                                       task_subtype_none, 0, 0, c, NULL);

      /* Gravity recursive down-pass */
      c->grav.down = scheduler_addtask(s, task_type_grav_down,
                                       task_subtype_none, 0, 0, c, NULL);

      /* Implicit tasks for the up and down passes */
      c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                           task_subtype_none, 0, 1, c, NULL);
      c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                          task_subtype_none, 0, 1, c, NULL);

      /* Link in the implicit tasks */
      scheduler_addunlock(s, c->grav.init, c->grav.init_out);
      scheduler_addunlock(s, c->grav.down_in, c->grav.down);
    }
  } else if ((c->grav.super != NULL) && c->subtype == cell_subtype_void) {

    /* Below the top level we just need to link in the implicit tasks. */

    if (is_self_gravity) {

      c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                           task_subtype_none, 0, 1, c, NULL);

      c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                          task_subtype_none, 0, 1, c, NULL);

      scheduler_addunlock(s, c->parent->grav.init_out, c->grav.init_out);
      scheduler_addunlock(s, c->grav.down_in, c->parent->grav.down_in);
    }
  } else if (c->type == cell_type_zoom) {

    /* Below the zoom top level we just need to link in the down implicit
     * tasks and handle the down unlocks between levels. */

    if (is_self_gravity) {

      /* Only create the down in if it doesn't exist. */
      if (c->grav.down_in == NULL) {
        c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                            task_subtype_none, 0, 1, c, NULL);
      }

      /* If we are at the zoom top level we need to link in the void parent
       * down_in. */
      if (c->top == c) {
        scheduler_addunlock(s, c->grav.down_in, c->void_parent->grav.down_in);
      } else {
        scheduler_addunlock(s, c->grav.down_in, c->parent->grav.down_in);
      }

      /* If we're at the super level, unlock the zoom super cell's grav down
       * after the void cell grav down and long range tasks. Afterwards we're
       * done recursing. */
      if (c->grav.super == c) {
        scheduler_addunlock(s, c->top->void_parent->grav.super->grav.down,
                            c->grav.down);
        scheduler_addunlock(s, c->grav.long_range,
                            c->top->void_parent->grav.super->grav.down);

        /* We also need to unlock the zoom's end_force after the void cells
         * down. */
        scheduler_addunlock(s, c->top->void_parent->grav.super->grav.down,
                            c->grav.super->grav.end_force);
      }
    }
  }

  /* Recurse but don't go deeper then the zoom super level. */
  if (c->subtype == cell_subtype_void || c->grav.super == NULL) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        zoom_engine_make_hierarchical_void_tasks_recursive(e, c->progeny[k]);
      }
    }
  }
}

/**
 * @brief Construct the hierarchical tasks for the void cell tree.
 *
 * This will construct:
 * - The long-range gravity task for the void cells.
 * - The down-pass gravity task for the void cells.
 * - The down-pass implicit task for the void cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void zoom_engine_make_hierarchical_void_tasks(struct engine *e) {

  ticks tic = getticks();

  /* Get a handle on the zoom properties. */
  struct space *s = e->s;
  struct zoom_region_properties *zoom_props = s->zoom_props;
  const int nr_void_cells = zoom_props->nr_void_cells;
  const int *void_cells = zoom_props->void_cells_top;
  struct cell *cells = s->cells_top;

  /* Loop through the void cells and make the hierarchical tasks. */
  for (int i = 0; i < nr_void_cells; i++) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Ensure we have a void cell. */
    if (cells[void_cells[i]].subtype != cell_subtype_void) {
      error("Cell is not a void cell.");
    }
#endif

    zoom_engine_make_hierarchical_void_tasks_recursive(e,
                                                       &cells[void_cells[i]]);
  }

  if (e->verbose)
    message("Making void cell tree tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for all combinations of cell types.
 *
 * - All top level cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 *
 * This is a wrapper around the various mappers defined above for all the
 * possible combinations of cell types including:
 * - zoom->zoom
 * - bkg->bkg
 * - zoom->bkg (if buffer cells are not used)
 * - buffer->buffer (if buffer cells are used)
 * - zoom->buffer (if buffer cells are used)
 * - buffer->bkg (if buffer cells are used)
 *
 * This replaces the function in engine_maketasks when running with a zoom
 * region.
 *
 * @param s The #space.
 * @param e The #engine.
 */
void zoom_engine_make_self_gravity_tasks(struct space *s, struct engine *e) {

  ticks tic = getticks();

  /* Zoom -> Zoom */
  threadpool_map(
      &e->threadpool, engine_make_self_gravity_tasks_mapper_zoom_cells, NULL,
      s->zoom_props->nr_zoom_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making zoom->zoom gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Background -> Background */
  threadpool_map(&e->threadpool,
                 engine_make_self_gravity_tasks_mapper_bkg_cells, NULL,
                 s->zoom_props->nr_bkg_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making bkg->bkg gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (s->zoom_props->with_buffer_cells) {

    tic = getticks();

    /* Buffer -> Buffer (only if we have a buffer region). */
    threadpool_map(
        &e->threadpool, engine_make_self_gravity_tasks_mapper_buffer_cells,
        NULL, s->zoom_props->nr_buffer_cells, 1, threadpool_auto_chunk_size, e);

    if (e->verbose)
      message("Making buffer->buffer gravity tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());

    tic = getticks();

    /* Buffer -> Background (only if we have a buffer region). */
    threadpool_map(
        &e->threadpool, engine_make_self_gravity_tasks_mapper_buffer_bkg, NULL,
        s->zoom_props->nr_buffer_cells, 1, threadpool_auto_chunk_size, e);

    if (e->verbose)
      message("Making buffer->bkg gravity tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }
}
