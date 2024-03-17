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
#include "zoom_region/zoom_maketasks.h"

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
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->cells_top;

  /* Some info about the zoom domain */
  const int bkg_cell_offset = s->zoom_props->bkg_cell_offset;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[bkg_cell_offset].width[0],
      s->max_softening, s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  const int delta =
      max((int)(sqrt(3) * distance / cells[bkg_cell_offset].width[0]) + 1, 2);
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

    /* Get the cell index, including background cell offset. */
    const int cid = (size_t)(map_data) + ind + bkg_cell_offset;

    /* Integer indices of the cell in the top-level grid */
    const int i = (cid - bkg_cell_offset) / (cdim[1] * cdim[2]);
    const int j = ((cid - bkg_cell_offset) / cdim[2]) % cdim[1];
    const int k = (cid - bkg_cell_offset) % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Ensure we haven't found a void cell with particles */
    if (ci->subtype == cell_subtype_void || ci->subtype == cell_subtype_empty)
      error("A void/empty cell (cid=%d) has got particles!", cid);

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Escape if non-periodic and beyond range */
      if (!periodic && (ii < 0 || ii >= cdim[0])) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

        /* Escape if non-periodic and beyond range */
        if (!periodic && (jj < 0 || jj >= cdim[1])) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Escape if non-periodic and beyond range */
          if (!periodic && (kk < 0 || kk >= cdim[2])) continue;

          /* Apply periodic BC (not harmful if not using periodic BC) */
          const int iii = (ii + cdim[0]) % cdim[0];
          const int jjj = (jj + cdim[1]) % cdim[1];
          const int kkk = (kk + cdim[2]) % cdim[2];

          /* Get the second cell */
          const int cjd =
              cell_getid_offset(cdim, bkg_cell_offset, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

#ifdef SWIFT_DEBUG_CHECKS
          /* Ensure both cells are background cells */
          if (ci->type != cell_type_bkg || cj->type != cell_type_bkg) {
            error(
                "Cell %d and cell %d are not both background cells! "
                "(ci->type=%d, cj->type=%d)",
                cid, cjd, cellID_names[ci->type], cellID_names[cj->type]);
          }
#endif

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Do we need a pair interaction for these cells? */
          if (engine_gravity_test_cell_pair(e, ci, cj)) {

            /* Ok, we need to add a direct pair calculation */
            engine_make_pair_gravity_task(e, sched, ci, cj, nodeID);
          }
        }
      }
    }
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
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const int cdim[3] = {s->zoom_props->buffer_cdim[0],
                       s->zoom_props->buffer_cdim[1],
                       s->zoom_props->buffer_cdim[2]};
  struct cell *cells = s->cells_top;

  /* Some info about the zoom domain */
  const int buffer_offset = s->zoom_props->buffer_cell_offset;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[buffer_offset].width[0],
      s->max_softening, s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
  const int delta =
      max((int)(sqrt(3) * distance / cells[buffer_offset].width[0]) + 1, 2);
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index, including background cell offset. */
    const int cid = (size_t)(map_data) + ind + buffer_offset;

    /* Integer indices of the cell in the top-level grid */
    const int i = (cid - buffer_offset) / (cdim[1] * cdim[2]);
    const int j = ((cid - buffer_offset) / cdim[2]) % cdim[1];
    const int k = (cid - buffer_offset) % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Ensure we haven't found a void cell with particles */
    if (ci->subtype == cell_subtype_void)
      error("A void cell (cid=%d) has got particles!", cid);

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta
     * (buffer cells don't care about periodicity) */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Buffer cells are never periodic. */
      if (ii < 0 || ii >= cdim[0]) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

        /* Buffer cells are never periodic. */
        if (jj < 0 || jj >= cdim[1]) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Buffer cells are never periodic. */
          if (kk < 0 || kk >= cdim[2]) continue;

          /* Get the second cell */
          const int cjd = cell_getid_offset(cdim, buffer_offset, ii, jj, kk);
          struct cell *cj = &cells[cjd];

#ifdef SWIFT_DEBUG_CHECKS
          /* Ensure both cells are buffer cells */
          if (ci->type != cell_type_buffer || cj->type != cell_type_buffer) {
            error(
                "Cell %d and cell %d are not both buffer cells! "
                "(ci->type=%d, cj->type=%d)",
                cid, cjd, cellID_names[ci->type], cellID_names[cj->type]);
          }
#endif

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Do we need a pair interaction for these cells? */
          if (engine_gravity_test_cell_pair(e, ci, cj)) {

            /* Ok, we need to add a direct pair calculation */
            engine_make_pair_gravity_task(e, sched, ci, cj, nodeID);
          }
        }
      }
    }
  }
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
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_self_gravity_tasks_mapper_zoom_cells(void *map_data,
                                                      int num_elements,
                                                      void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const int cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                       s->zoom_props->cdim[2]};
  struct cell *cells = s->cells_top;

  /* Compute maximal distance where we can expect a direct interaction */
  const float distance = gravity_M2L_min_accept_distance(
      e->gravity_properties, sqrtf(3) * cells[0].width[0], s->max_softening,
      s->min_a_grav, s->max_mpole_power, periodic);

  /* Convert the maximal search distance to a number of cells
   * Define a lower and upper delta in case things are not symmetric */
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

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the first cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta
     * (zoom cells don't care about periodicity) */
    for (int ii = i - delta_m; ii <= i + delta_p; ii++) {

      /* Zoom cells are never periodic, exit if beyond zoom region */
      if (ii < 0 || ii >= cdim[0]) continue;

      for (int jj = j - delta_m; jj <= j + delta_p; jj++) {

        /* Zoom cells are never periodic, exit if beyond zoom region */
        if (jj < 0 || jj >= cdim[1]) continue;

        for (int kk = k - delta_m; kk <= k + delta_p; kk++) {

          /* Zoom cells are never periodic, exit if beyond zoom region */
          if (kk < 0 || kk >= cdim[2]) continue;

          /* Get the second cell */
          const int cjd = cell_getid(cdim, ii, jj, kk);
          struct cell *cj = &cells[cjd];

#ifdef SWIFT_DEBUG_CHECKS
          /* Ensure both cells are zoom cells */
          if (ci->type != cell_type_zoom || cj->type != cell_type_zoom) {
            error(
                "Cell %d and cell %d are not both zoom cells! "
                "(ci->type=%d, cj->type=%d)",
                cid, cjd, cellID_names[ci->type], cellID_names[cj->type]);
          }
#endif

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Do we need a pair interaction for these cells? */
          if (engine_gravity_test_cell_pair(e, ci, cj)) {

            /* Ok, we need to add a direct pair calculation */
            engine_make_pair_gravity_task(e, sched, ci, cj, nodeID);
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions between zoom cells
 * and background cells.
 *
 * This mapper only considers interactions between zoom cells and their
 * neighbours (i.e. the cells within the gravity criterion distance
 * above them in hierarchy). When buffer cells are present, interactions between
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
void engine_make_self_gravity_tasks_mapper_zoom_bkg(void *map_data,
                                                    int num_elements,
                                                    void *extra_data) {

  /* Useful local information */
  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Get the neighbouring background cells. */
  const int nr_neighbours = s->zoom_props->nr_neighbour_cells;
  const int *neighbour_cells = s->zoom_props->neighbour_cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Loop over every neighbouring background cells */
    for (int k = 0; k < nr_neighbours; k++) {

      /* Get the cell */
      int cjd = neighbour_cells[k];
      struct cell *cj = &cells[cjd];

      /* Avoid duplicates, empty cells and completely foreign pairs */
      if (cid >= cjd || cj->grav.count == 0 ||
          (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

#ifdef SWIFT_DEBUG_CHECKS
      /* Ensure both cells are not in the same level */
      if (ci->type == cj->type) {
        error(
            "Cell %d and cell %d are the same cell type "
            "(One should be a neighbour)! "
            "(ci->type=%d, cj->type=%d)",
            cid, cjd, ci->type, cj->type);
      }
#endif

      /* Do we need a pair interaction for these cells? */
      if (engine_gravity_test_cell_pair(e, ci, cj)) {

        /* Ok, we need to add a direct pair calculation */
        engine_make_pair_gravity_task(e, sched, ci, cj, nodeID);
      }
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
      if (engine_gravity_test_cell_pair(e, ci, cj)) {

        /* Ok, we need to add a direct pair calculation */
        engine_make_pair_gravity_task(e, sched, ci, cj, nodeID);
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

  threadpool_map(
      &e->threadpool, engine_make_self_gravity_tasks_mapper_zoom_cells, NULL,
      s->zoom_props->nr_zoom_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making zoom->zoom gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  threadpool_map(&e->threadpool,
                 engine_make_self_gravity_tasks_mapper_bkg_cells, NULL,
                 s->zoom_props->nr_bkg_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making bkg->bkg gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (s->zoom_props->with_buffer_cells) {

    tic = getticks();

    threadpool_map(
        &e->threadpool, engine_make_self_gravity_tasks_mapper_buffer_cells,
        NULL, s->zoom_props->nr_buffer_cells, 1, threadpool_auto_chunk_size, e);

    if (e->verbose)
      message("Making buffer->buffer gravity tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }

  tic = getticks();

  threadpool_map(&e->threadpool, engine_make_self_gravity_tasks_mapper_zoom_bkg,
                 NULL, s->zoom_props->nr_zoom_cells, 1,
                 threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making zoom->buffer gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  if (s->zoom_props->with_buffer_cells) {

    tic = getticks();

    threadpool_map(
        &e->threadpool, engine_make_self_gravity_tasks_mapper_buffer_bkg, NULL,
        s->zoom_props->nr_buffer_cells, 1, threadpool_auto_chunk_size, e);

    if (e->verbose)
      message("Making buffer->bkg gravity tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());
  }
}

/**
 * @brief Constructs the top-level pair tasks for the first hydro loop over
 * neighbours
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * This replaces the function in engine_maketasks but simply removes
 * periodicity (the zoom level never has periodicity), otherwise the function
 * is identical.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_hydroloop_tasks_mapper_with_zoom(void *map_data,
                                                  int num_elements,
                                                  void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_sinks = (e->policy & engine_policy_sinks);
  const int with_black_holes = (e->policy & engine_policy_black_holes);

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->zoom_props->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without hydro or star particles */
    if ((ci->hydro.count == 0) && (!with_stars || ci->stars.count == 0) &&
        (!with_sinks || ci->sinks.count == 0) &&
        (!with_black_holes || ci->black_holes.count == 0))
      continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_density, 0, 0, ci,
                        NULL);
    }

    /* Now loop over all the neighbours of this cell */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (iii < 0 || iii >= cdim[0]) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (jjj < 0 || jjj >= cdim[1]) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (kkk < 0 || kkk >= cdim[2]) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Is that neighbour local and does it have gas or star particles ?
           */
          if ((cid >= cjd) ||
              ((cj->hydro.count == 0) &&
               (!with_feedback || cj->stars.count == 0) &&
               (!with_sinks || cj->sinks.count == 0) &&
               (!with_black_holes || cj->black_holes.count == 0)) ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Construct the pair task */
          const int sid = sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
          scheduler_addtask(sched, task_type_pair, task_subtype_density, sid, 0,
                            ci, cj);

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI

          /* Let's cross-check that we had a proxy for that cell */
          if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

            /* Find the proxy for this node */
            const int proxy_id = e->proxy_ind[cj->nodeID];
            if (proxy_id < 0)
              error("No proxy exists for that foreign node %d!", cj->nodeID);

            const struct proxy *p = &e->proxies[proxy_id];

            /* Check whether the cell exists in the proxy */
            int n = 0;
            for (n = 0; n < p->nr_cells_in; n++)
              if (p->cells_in[n] == cj) break;
            if (n == p->nr_cells_in)
              error(
                  "Cell %d not found in the proxy but trying to construct "
                  "hydro task!",
                  cjd);
          } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

            /* Find the proxy for this node */
            const int proxy_id = e->proxy_ind[ci->nodeID];
            if (proxy_id < 0)
              error("No proxy exists for that foreign node %d!", ci->nodeID);

            const struct proxy *p = &e->proxies[proxy_id];

            /* Check whether the cell exists in the proxy */
            int n = 0;
            for (n = 0; n < p->nr_cells_in; n++)
              if (p->cells_in[n] == ci) break;
            if (n == p->nr_cells_in)
              error(
                  "Cell %d not found in the proxy but trying to construct "
                  "hydro task!",
                  cid);
          }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level self + pair tasks for the FOF loop over
 * neighbours.
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * This replaces the function in engine_maketasks but simply removes
 * periodicity (the zoom level never has periodicity), otherwise the function
 * is identical.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_fofloop_tasks_mapper_with_zoom(void *map_data,
                                                int num_elements,
                                                void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->zoom_props->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cells is local build a self-interaction */
    if (ci->nodeID == nodeID)
      scheduler_addtask(sched, task_type_fof_self, task_subtype_none, 0, 0, ci,
                        NULL);
    else
      continue;

    /* Now loop over all the neighbours of this cell */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (iii < 0 || iii >= cdim[0]) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (jjj < 0 || jjj >= cdim[1]) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (kkk < 0 || kkk >= cdim[2]) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Is that neighbour local and does it have particles ? */
          if (cid >= cjd || cj->grav.count == 0 || (ci->nodeID != cj->nodeID))
            continue;

          /* Construct the pair task */
          scheduler_addtask(sched, task_type_fof_pair, task_subtype_none, 0, 0,
                            ci, cj);
        }
      }
    }
  }
}
