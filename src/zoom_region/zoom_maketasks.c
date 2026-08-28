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
#include "multipole_accept.h"
#include "proxy.h"
#include "space.h"
#include "zoom_region/zoom.h"

/**
 * @brief Compute the maximal distance at which a pair of background cells can
 *        still require a direct interaction.
 *
 * The multipole acceptance criterion combines the maximal multipole powers,
 * the maximal softening and the minimal particle acceleration with the size of
 * the cells being paired. In a zoom simulation the zoom and background grids
 * have different widths and their particle populations differ by orders of
 * magnitude in mass, softening and acceleration. A single reduction over every
 * top-level cell therefore pairs a quantity drawn from one grid with the cell
 * size of the other, which is meaningless. The background search uses the
 * background reduction only.
 *
 * @param s The #space.
 * @param props The properties of the gravity scheme.
 *
 * @return The maximal distance for a direct interaction.
 */
float zoom_bkg_M2L_min_accept_distance(const struct space *s,
                                       const struct gravity_props *props) {

  const struct cell *bkg_cells = s->zoom_props->bkg_cells_top;

  return gravity_M2L_min_accept_distance(
      props, sqrtf(3.f) * bkg_cells[0].width[0], s->bkg_max_softening,
      s->bkg_min_a_grav, s->bkg_max_mpole_power, s->periodic);
}

/**
 * @brief Compute the maximal distance at which a pair of zoom cells can still
 *        require a direct interaction.
 *
 * The zoom counterpart of zoom_bkg_M2L_min_accept_distance, built entirely
 * from the zoom cell width and the zoom particle population.
 *
 * @param s The #space.
 * @param props The properties of the gravity scheme.
 *
 * @return The maximal distance for a direct interaction.
 */
float zoom_zoom_M2L_min_accept_distance(const struct space *s,
                                        const struct gravity_props *props) {

  const struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

  return gravity_M2L_min_accept_distance(
      props, sqrtf(3.f) * zoom_cells[0].width[0], s->zoom_max_softening,
      s->zoom_min_a_grav, s->zoom_max_mpole_power, s->periodic);
}

/**
 * @brief Convert the background search distance into a number of cells.
 *
 * Defines a lower and upper delta in case things are not symmetric. If every
 * cell is in range of every other one the deltas are clamped so that the pair
 * loops cover the grid exactly once.
 *
 * When the mesh is in use the search is additionally bounded by the mesh
 * cut-off. Beyond r_cut_max the truncated forces are zero and cell_can_use_mesh
 * hands the pair to the mesh unconditionally, before any multipole acceptance
 * test is reached, so no pair further apart than that can ever need a task.
 * This is the same ball the long-range task searches over this same grid (see
 * runner_do_grav_long_range_zoom_periodic), and the two are complements: the
 * long-range task performs the M-M interactions that the pair tasks do not.
 * Searching further here creates work the long-range task will never ask for.
 *
 * NOTE: The 2 in the max below may not be necessary but does insure some
 * safety buffer.
 *
 * @param e The #engine.
 * @param delta_m (return) The number of cells to search in the -ve direction.
 * @param delta_p (return) The number of cells to search in the +ve direction.
 */
void zoom_bkg_gravity_search_delta(const struct engine *e, int *delta_m,
                                   int *delta_p) {

  const struct space *s = e->s;
  const struct cell *bkg_cells = s->zoom_props->bkg_cells_top;
  const int cdim = s->cdim[0];

  const float distance =
      zoom_bkg_M2L_min_accept_distance(s, e->gravity_properties);

  int delta = max((int)(sqrt(3) * distance / bkg_cells[0].width[0]) + 1, 2);

  /* Bound the search by the mesh cut-off. */
  if (s->periodic) {
    const int delta_mesh =
        ceil(e->mesh->r_cut_max *
             max3(s->iwidth[0], s->iwidth[1], s->iwidth[2])) +
        1;
    if (delta_mesh < delta) delta = delta_mesh;
  }

  *delta_m = delta;
  *delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (s->periodic) {
    if (delta >= cdim / 2) {
      if (cdim % 2 == 0) {
        *delta_m = cdim / 2;
        *delta_p = cdim / 2 - 1;
      } else {
        *delta_m = cdim / 2;
        *delta_p = cdim / 2;
      }
    }
  } else {
    if (delta > cdim) {
      *delta_m = cdim;
      *delta_p = cdim;
    }
  }
}

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
 * This will create pair tasks between void cells and background cells. These
 * pair tasks will be split into smaller
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

  /* Unlike zoom cells, background cells are periodic at the box boundaries if
   * the space is periodic. */
  const int periodic = s->periodic;

  /* Compute the range of background cells we need to search, using the
   * background cell width and the background particle population. */
  int delta_m, delta_p;
  zoom_bkg_gravity_search_delta(e, &delta_m, &delta_p);

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    const int cid = (size_t)(map_data) + ind;
    struct cell *c = &cells[cid];

    /* Skip void cells that do not contain any zoom cells. We do not want
     * to generate self/pair gravity tasks for these top-level cells. */
    if (c->subtype == cell_subtype_void && !c->contains_zoom_cells) continue;

    /* Create a self task, and loop over neighbouring cells making pair tasks
     * where appropriate. */
    engine_gravity_make_task_loop(e, cid, cdim, cells, periodic, use_mesh,
                                  delta_m, delta_p);
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
 * - bkg->bkg
 *
 * This replaces the function in engine_maketasks when running with a zoom
 * region.
 *
 * @param s The #space.
 * @param e The #engine.
 */
void zoom_engine_make_self_gravity_tasks(struct space *s, struct engine *e) {

  ticks tic = getticks();

  /* Report the search range the background mappers will use, alongside the
   * equivalent zoom quantity. Each is built from its own cell width and its
   * own particle population. When delta saturates at cdim/2 the pair loops
   * degenerate to all-pairs. */
  if (e->verbose) {
    const struct cell *bkg_cells = s->zoom_props->bkg_cells_top;
    const struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

    const float bkg_distance =
        zoom_bkg_M2L_min_accept_distance(s, e->gravity_properties);
    const float zoom_distance =
        zoom_zoom_M2L_min_accept_distance(s, e->gravity_properties);

    int delta_m, delta_p;
    zoom_bkg_gravity_search_delta(e, &delta_m, &delta_p);
    const int bkg_delta =
        max((int)(sqrt(3) * bkg_distance / bkg_cells[0].width[0]) + 1, 2);
    const int zoom_delta =
        max((int)(sqrt(3) * zoom_distance / zoom_cells[0].width[0]) + 1, 2);

    message(
        "Background pair search: distance=%.6e (%.2f bkg cell widths) "
        "delta=%d (clamped to -%d/+%d), bkg_cdim=%d, saturates at delta>=%d%s",
        bkg_distance, bkg_distance / bkg_cells[0].width[0], bkg_delta, delta_m,
        delta_p, s->cdim[0], s->cdim[0] / 2,
        bkg_delta >= s->cdim[0] / 2 ? " <== SATURATED" : "");
    message(
        "Zoom pair search:       distance=%.6e (%.2f zoom cell widths) "
        "delta=%d, zoom_cdim=%d",
        zoom_distance, zoom_distance / zoom_cells[0].width[0], zoom_delta,
        s->zoom_props->cdim[0]);
  }

  /* Background -> Background */
  threadpool_map(&e->threadpool,
                 engine_make_self_gravity_tasks_mapper_bkg_cells, NULL,
                 s->zoom_props->nr_bkg_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making bkg->bkg gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
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
void zoom_engine_make_hierarchical_void_tasks_recursive(struct engine *e,
                                                        struct cell *c) {

  struct scheduler *s = &e->sched;

  /* Nothing to do if there's no gravity. */
  if (!(e->policy & engine_policy_self_gravity)) return;

  /* At the super level we have a few different tasks to make. (We don't need
   * any tasks above the super level) */
  if (c->grav.super == c) {

    /* Initialisation of the multipoles */
    c->grav.init = scheduler_addtask(s, task_type_init_grav, task_subtype_none,
                                     0, 0, c, NULL);

    /* Gravity recursive down-pass */
    c->grav.down = scheduler_addtask(s, task_type_grav_down, task_subtype_none,
                                     0, 0, c, NULL);

    /* Implicit tasks for the up and down passes */
    c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                         task_subtype_none, 0, 1, c, NULL);
    c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                        task_subtype_none, 0, 1, c, NULL);

    /* Link in the implicit tasks */
    scheduler_addunlock(s, c->grav.init, c->grav.init_out);
    scheduler_addunlock(s, c->grav.down_in, c->grav.down);

  } else if (c->grav.super != NULL) {

    /* Below the super level we just need to link in the implicit tasks. */
    c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                         task_subtype_none, 0, 1, c, NULL);
    c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                        task_subtype_none, 0, 1, c, NULL);

    scheduler_addunlock(s, c->parent->grav.init_out, c->grav.init_out);
    scheduler_addunlock(s, c->grav.down_in, c->parent->grav.down_in);
  }

  /* Recurse but only in void cells. */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL && c->progeny[k]->subtype == cell_subtype_void &&
        c->progeny[k]->contains_zoom_cells) {
      zoom_engine_make_hierarchical_void_tasks_recursive(e, c->progeny[k]);
    }
  }
}

/**
 * @brief Construct the hierarchical tasks for the void cell tree.
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
void zoom_engine_make_hierarchical_void_tasks(struct engine *e) {

  ticks tic = getticks();

  /* Get a handle on the zoom properties. */
  struct space *s = e->s;
  struct zoom_region_properties *zoom_props = s->zoom_props;
  const int nr_void_cells = zoom_props->nr_void_cells;
  const int *void_cells = zoom_props->void_cell_indices;
  struct cell *cells = s->cells_top;

  /* Loop through the void cells and make the hierarchical tasks. */
  for (int i = 0; i < nr_void_cells; i++) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Ensure we have a void cell. */
    if (cells[void_cells[i]].subtype != cell_subtype_void) {
      error("Cell is not a void cell.");
    }
#endif

    /* Get the void cell. */
    struct cell *c = &cells[void_cells[i]];

    /* Skip "non-useful" void cells that do not contain any zoom cells.
     * We don't want, nor need, tasks for these. */
    if (c->contains_zoom_cells == 0) continue;

    zoom_engine_make_hierarchical_void_tasks_recursive(e, c);
  }

  if (e->verbose)
    message("Making void cell tree tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}
