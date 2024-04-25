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
 * @brief Recurses through the void cells and constructs the top-level tasks
 * at the zoom level but only if an mm interaction couldn't be used higher in
 * the void hierarchy.
 *
 * @param e The #engine.
 * @param ci The neighbour #cell.
 * @param cj The void or zoom #cell to test for a pair task.
 * @param periodic The periodicity of the space.
 * @param use_mesh Are we using the mesh?
 */
static void engine_make_self_gravity_tasks_void_cell_recursive(struct engine *e,
                                                               struct cell *ci,
                                                               struct cell *cj,
                                                               int periodic,
                                                               int use_mesh) {

  /* At least one of these must be on this rank. For the void cell here we need
   * the nodeID to be a valid rank, < 0 indicates it's zoom leaves are on
   * different ranks and we can't know for sure that we can exit. */
  if (ci->nodeID != e->nodeID && (cj->nodeID >= 0 && cj->nodeID != e->nodeID))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that cj is either a void cell or a zoom cell. */
  if (cj->subtype != cell_subtype_void && cj->type != cell_type_zoom) {
    error(
        "Invalid cell type in "
        "engine_make_self_gravity_tasks_void_cell_recursive. (cj->type=%s, "
        "cj->subtype=%s)",
        cellID_names[cj->type], subcellID_names[cj->subtype]);
  }
#endif

  /* Recover the multipole information */
  const struct gravity_tensors *multi_j = cj->grav.multipole;

  /* Skip empty cells */
  if (multi_j->m_pole.M_000 == 0.f) return;

  /* Do we need a pair interaction for these cells? */
  if (engine_gravity_need_cell_pair_task(e, ci, cj, periodic, use_mesh)) {

    /* If so we can only make a direct interaction at the zoom level. */
    if (cj->type == cell_type_zoom) {

      /* Exit if there are no particles. */
      if (cj->grav.count == 0) return;

      /* Get the indices of ci and cj. */
      const int cid = ci - e->s->cells_top;
      const int cjd = cj - e->s->cells_top;

      /* Ok, we need to add a direct pair calculation */
      engine_make_pair_gravity_task(e, &e->sched, ci, cj, e->nodeID, cid, cjd);

    } else {
      /* Since we aren't at the zoom level we must recurse. */
      for (int k = 0; k < 8; k++) {
        engine_make_self_gravity_tasks_void_cell_recursive(
            e, ci, cj->progeny[k], periodic, use_mesh);
      }
    }
  }

  /* If not then this interaction will be handled by an mm interction in the
   * long range gravity task or the mesh. */
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions for zoom cells.
 *
 * This mapper only considers zoom->zoom interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 *
 * Unlike the other mappers considered with pairs of cells with the same type,
 * this mapper recurses through the void heirarchy to capture any interactions
 * which can be handled by a mm interaction in the long-range gravity task above
 * the zoom level.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  const int cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                       s->zoom_props->cdim[2]};
  struct cell *cells = s->zoom_props->zoom_cells_top;

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
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between zoom cells
 * and background cells.
 *
 * This mapper only considers interactions between zoom cells and their
 * neighbours (i.e. the cells within the gravity criterion distance
 * above them in hierarchy). When buffer cells are present, interactions
 between
 * zoom cells and buffer cells are considered, otherwise only interactions
 * between zoom cells and background cells are considered.
 *
 * - All top-cells get a self task.
 * - All pairs of differing sized cells within range according to
 *   the multipole acceptance criterion get a pair task.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.

 */
void engine_make_self_gravity_tasks_mapper_zoom_neighbour(void *map_data,
                                                          int num_elements,
                                                          void *extra_data) {

  /* Useful local information */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;

  /* We always use the mesh if the volume is periodic. */
  const int use_mesh = s->periodic;

  /* Zoom->Bkg interactions will never cross the box boundary. */
  const int periodic = 0;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Get the neighbouring background cells. */
  const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

  /* Get the void cells. */
  const int *void_cells = s->zoom_props->void_cells_top;
  const int nr_void_cells = s->zoom_props->nr_void_cells;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *ci = &cells[neighbour_cells[cid]];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Loop over every void cell. */
    for (int k = 0; k < nr_void_cells; k++) {

      /* Get the cell */
      int cjd = void_cells[k];
      struct cell *cj = &cells[cjd];

      /* Recurse through the void cell tree to see where we need pair tasks. */
      engine_make_self_gravity_tasks_void_cell_recursive(e, ci, cj, periodic,
                                                         use_mesh);
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between natural level cells
 * and zoom level cells.
 *
 * This mapper only consider buffer->bkg interactions. It will only be called
 * if buffer cells are used.
 *
 * - All top-cells get a self task.
 * - All pairs of differing sized cells within range according to
 *   the multipole acceptance criterion get a pair task.
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

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Loop over every neighbouring background cells */
    for (int cjd = bkg_offset; cjd < buffer_offset; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

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

  /* Report the distance at which we'll search for zoom->zoom pair tasks. */
  if (e->verbose) {
    const float distance = gravity_M2L_min_accept_distance(
        e->gravity_properties,
        sqrtf(3) * s->zoom_props->zoom_cells_top[0].width[0], s->max_softening,
        s->min_a_grav, s->max_mpole_power, /*periodic*/ 0);
    message("Searching for zoom->zoom pair tasks within %.3f U_L", distance);
  }

  /* Zoom -> Zoom */
  threadpool_map(
      &e->threadpool, engine_make_self_gravity_tasks_mapper_zoom_cells, NULL,
      s->zoom_props->nr_zoom_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making zoom->zoom gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Report the distance at which we'll search for bkg->bkg pair tasks. */
  if (e->verbose) {
    const float distance = gravity_M2L_min_accept_distance(
        e->gravity_properties,
        sqrtf(3) * s->zoom_props->bkg_cells_top[0].width[0], s->max_softening,
        s->min_a_grav, s->max_mpole_power, s->periodic);
    message("Searching for bkg->bkg pair tasks within %.3f U_L", distance);
  }

  /* Background -> Background */
  threadpool_map(&e->threadpool,
                 engine_make_self_gravity_tasks_mapper_bkg_cells, NULL,
                 s->zoom_props->nr_bkg_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making bkg->bkg gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (s->zoom_props->with_buffer_cells) {

    tic = getticks();

    /* Report the distance at which we'll search for buffer->buffer pair tasks.
     */
    if (e->verbose) {
      const float distance = gravity_M2L_min_accept_distance(
          e->gravity_properties,
          sqrtf(3) * s->zoom_props->buffer_cells_top[0].width[0],
          s->max_softening, s->min_a_grav, s->max_mpole_power, /*periodic*/ 0);
      message("Searching for buffer->buffer pair tasks within %.3f U_L",
              distance);
    }

    /* Buffer -> Buffer (only if we have a buffer region). */
    threadpool_map(
        &e->threadpool, engine_make_self_gravity_tasks_mapper_buffer_cells,
        NULL, s->zoom_props->nr_buffer_cells, 1, threadpool_auto_chunk_size, e);

    if (e->verbose)
      message("Making buffer->buffer gravity tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

  tic = getticks();

  /* Zoom -> Neighbour (A neighbour is a Buffer cell if we have a buffer region,
   * otherwise it's a Background cell). */
  threadpool_map(&e->threadpool,
                 engine_make_self_gravity_tasks_mapper_zoom_neighbour, NULL,
                 s->zoom_props->nr_neighbour_cells, 1,
                 threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making zoom->buffer gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (s->zoom_props->with_buffer_cells) {

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
