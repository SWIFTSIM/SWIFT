/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/**
 * @brief Recursively dismantle a cell tree.
 *
 * @param s The #space.
 * @param c The #cell to recycle.
 * @param cell_rec_begin Pointer to the start of the list of cells to recycle.
 * @param cell_rec_end Pointer to the end of the list of cells to recycle.
 * @param multipole_rec_begin Pointer to the start of the list of multipoles to
 * recycle.
 * @param multipole_rec_end Pointer to the end of the list of multipoles to
 * recycle.
 */
void space_rebuild_recycle_rec(struct space *s, struct cell *c,
                               struct cell **cell_rec_begin,
                               struct cell **cell_rec_end,
                               struct gravity_tensors **multipole_rec_begin,
                               struct gravity_tensors **multipole_rec_end) {
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        space_rebuild_recycle_rec(s, c->progeny[k], cell_rec_begin,
                                  cell_rec_end, multipole_rec_begin,
                                  multipole_rec_end);

        c->progeny[k]->next = *cell_rec_begin;
        *cell_rec_begin = c->progeny[k];

        if (s->with_self_gravity) {
          c->progeny[k]->grav.multipole->next = *multipole_rec_begin;
          *multipole_rec_begin = c->progeny[k]->grav.multipole;
        }

        if (*cell_rec_end == NULL) *cell_rec_end = *cell_rec_begin;
        if (s->with_self_gravity && *multipole_rec_end == NULL)
          *multipole_rec_end = *multipole_rec_begin;

        c->progeny[k]->grav.multipole = NULL;
        c->progeny[k] = NULL;
      }
}

void space_rebuild_recycle_mapper(void *map_data, int num_elements,
                                  void *extra_data) {

  struct space *s = (struct space *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int k = 0; k < num_elements; k++) {
    struct cell *c = &cells[k];
    struct cell *cell_rec_begin = NULL, *cell_rec_end = NULL;
    struct gravity_tensors *multipole_rec_begin = NULL,
                           *multipole_rec_end = NULL;
    space_rebuild_recycle_rec(s, c, &cell_rec_begin, &cell_rec_end,
                              &multipole_rec_begin, &multipole_rec_end);
    if (cell_rec_begin != NULL)
      space_recycle_list(s, cell_rec_begin, cell_rec_end, multipole_rec_begin,
                         multipole_rec_end);
    c->hydro.sorts = NULL;
    c->stars.sorts = NULL;
    c->nr_tasks = 0;
    c->grav.nr_mm_tasks = 0;
    c->hydro.density = NULL;
    c->hydro.gradient = NULL;
    c->hydro.force = NULL;
    c->hydro.limiter = NULL;
    c->grav.grav = NULL;
    c->grav.mm = NULL;
    c->hydro.dx_max_part = 0.0f;
    c->hydro.dx_max_sort = 0.0f;
    c->sinks.dx_max_part = 0.f;
    c->stars.dx_max_part = 0.f;
    c->stars.dx_max_sort = 0.f;
    c->black_holes.dx_max_part = 0.f;
    c->hydro.sorted = 0;
    c->hydro.sort_allocated = 0;
    c->stars.sorted = 0;
    c->hydro.count = 0;
    c->hydro.count_total = 0;
    c->hydro.updated = 0;
    c->grav.count = 0;
    c->grav.count_total = 0;
    c->grav.updated = 0;
    c->sinks.count = 0;
    c->stars.count = 0;
    c->stars.count_total = 0;
    c->stars.updated = 0;
    c->black_holes.count = 0;
    c->black_holes.count_total = 0;
    c->black_holes.updated = 0;
    c->grav.init = NULL;
    c->grav.init_out = NULL;
    c->hydro.extra_ghost = NULL;
    c->hydro.ghost_in = NULL;
    c->hydro.ghost_out = NULL;
    c->hydro.ghost = NULL;
    c->hydro.prep1_ghost = NULL;
    c->hydro.star_formation = NULL;
    c->sinks.sink_formation = NULL;
    c->sinks.star_formation_sink = NULL;
    c->hydro.stars_resort = NULL;
    c->stars.density_ghost = NULL;
    c->stars.prep1_ghost = NULL;
    c->stars.prep2_ghost = NULL;
    c->stars.density = NULL;
    c->stars.feedback = NULL;
    c->stars.prepare1 = NULL;
    c->stars.prepare2 = NULL;
    c->sinks.swallow = NULL;
    c->sinks.do_sink_swallow = NULL;
    c->sinks.do_gas_swallow = NULL;
    c->black_holes.density_ghost = NULL;
    c->black_holes.swallow_ghost_1 = NULL;
    c->black_holes.swallow_ghost_2 = NULL;
    c->black_holes.swallow_ghost_3 = NULL;
    c->black_holes.density = NULL;
    c->black_holes.swallow = NULL;
    c->black_holes.do_gas_swallow = NULL;
    c->black_holes.do_bh_swallow = NULL;
    c->black_holes.feedback = NULL;
#ifdef WITH_CSDS
    c->csds = NULL;
#endif
    c->kick1 = NULL;
    c->kick2 = NULL;
    c->timestep = NULL;
    c->timestep_limiter = NULL;
    c->timestep_sync = NULL;
    c->timestep_collect = NULL;
    c->hydro.end_force = NULL;
    c->hydro.drift = NULL;
    c->sinks.drift = NULL;
    c->stars.drift = NULL;
    c->stars.stars_in = NULL;
    c->stars.stars_out = NULL;
    c->black_holes.drift = NULL;
    c->black_holes.black_holes_in = NULL;
    c->black_holes.black_holes_out = NULL;
    c->sinks.sink_in = NULL;
    c->sinks.sink_ghost1 = NULL;
    c->sinks.sink_ghost2 = NULL;
    c->sinks.sink_out = NULL;
    c->grav.drift = NULL;
    c->grav.drift_out = NULL;
    c->hydro.cooling_in = NULL;
    c->hydro.cooling_out = NULL;
    c->hydro.cooling = NULL;
    c->grav.long_range = NULL;
    c->grav.down_in = NULL;
    c->grav.down = NULL;
    c->grav.end_force = NULL;
    c->grav.neutrino_weight = NULL;
    c->top = c;
    c->super = c;
    c->hydro.super = c;
    c->grav.super = c;
    c->hydro.parts = NULL;
    c->hydro.xparts = NULL;
    c->grav.parts = NULL;
    c->grav.parts_rebuild = NULL;
    c->sinks.parts = NULL;
    c->stars.parts = NULL;
    c->stars.parts_rebuild = NULL;
    c->black_holes.parts = NULL;
    c->flags = 0;
    c->hydro.ti_end_min = -1;
    c->grav.ti_end_min = -1;
    c->sinks.ti_end_min = -1;
    c->stars.ti_end_min = -1;
    c->black_holes.ti_end_min = -1;
    c->rt.rt_in = NULL;
    c->rt.rt_ghost1 = NULL;
    c->rt.rt_gradient = NULL;
    c->rt.rt_ghost2 = NULL;
    c->rt.rt_transport = NULL;
    c->rt.rt_transport_out = NULL;
    c->rt.rt_tchem = NULL;
    c->rt.rt_advance_cell_time = NULL;
    c->rt.rt_sorts = NULL;
    c->rt.rt_out = NULL;
    c->rt.rt_collect_times = NULL;
    c->rt.ti_rt_end_min = -1;
    c->rt.ti_rt_min_step_size = -1;
    c->rt.updated = 0;
#ifdef SWIFT_RT_DEBUG_CHECKS
    c->rt.advanced_time = 0;
#endif

    star_formation_logger_init(&c->stars.sfh);
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
    c->cellID = 0;
#endif
    if (s->with_self_gravity)
      bzero(c->grav.multipole, sizeof(struct gravity_tensors));

    cell_free_hydro_sorts(c);
    cell_free_stars_sorts(c);
#if WITH_MPI
    c->mpi.tag = -1;
    c->mpi.recv = NULL;
    c->mpi.send = NULL;
#endif
  }
}

/**
 * @brief Return a used cell to the buffer of unused sub-cells.
 *
 * @param s The #space.
 * @param c The #cell.
 */
void space_recycle(struct space *s, struct cell *c) {

  /* Clear the cell. */
  if (lock_destroy(&c->hydro.lock) != 0 || lock_destroy(&c->grav.plock) != 0 ||
      lock_destroy(&c->grav.mlock) != 0 || lock_destroy(&c->stars.lock) != 0 ||
      lock_destroy(&c->sinks.lock) != 0 ||
      lock_destroy(&c->sinks.sink_formation_lock) != 0 ||
      lock_destroy(&c->black_holes.lock) != 0 ||
      lock_destroy(&c->grav.star_formation_lock) != 0 ||
      lock_destroy(&c->stars.star_formation_lock) != 0)
    error("Failed to destroy spinlocks.");

  /* Hook the multipole back in the buffer */
  if (s->with_self_gravity) {
    c->grav.multipole->next = s->multipoles_sub[c->tpid];
    s->multipoles_sub[c->tpid] = c->grav.multipole;
  }

  /* Hook this cell into the buffer. */
  c->next = s->cells_sub[c->tpid];
  s->cells_sub[c->tpid] = c;
  atomic_dec(&s->tot_cells);
}

/**
 * @brief Return a list of used cells to the buffer of unused sub-cells.
 *
 * @param s The #space.
 * @param cell_list_begin Pointer to the first #cell in the linked list of
 *        cells joined by their @c next pointers.
 * @param cell_list_end Pointer to the last #cell in the linked list of
 *        cells joined by their @c next pointers. It is assumed that this
 *        cell's @c next pointer is @c NULL.
 * @param multipole_list_begin Pointer to the first #multipole in the linked
 * list of
 *        multipoles joined by their @c next pointers.
 * @param multipole_list_end Pointer to the last #multipole in the linked list
 * of
 *        multipoles joined by their @c next pointers. It is assumed that this
 *        multipole's @c next pointer is @c NULL.
 */
void space_recycle_list(struct space *s, struct cell *cell_list_begin,
                        struct cell *cell_list_end,
                        struct gravity_tensors *multipole_list_begin,
                        struct gravity_tensors *multipole_list_end) {

  int count = 0;

  /* Clean up the list of cells. */
  for (struct cell *c = cell_list_begin; c != NULL; c = c->next) {
    /* Clear the cell. */
    if (lock_destroy(&c->hydro.lock) != 0 ||
        lock_destroy(&c->grav.plock) != 0 ||
        lock_destroy(&c->grav.mlock) != 0 ||
        lock_destroy(&c->stars.lock) != 0 ||
        lock_destroy(&c->sinks.lock) != 0 ||
        lock_destroy(&c->sinks.sink_formation_lock) != 0 ||
        lock_destroy(&c->black_holes.lock) != 0 ||
        lock_destroy(&c->stars.star_formation_lock) != 0 ||
        lock_destroy(&c->grav.star_formation_lock) != 0)
      error("Failed to destroy spinlocks.");

    /* Count this cell. */
    count += 1;
  }

  /* Lock the space. */
  lock_lock(&s->lock);

  /* Hook the cells into the buffer keeping tpid if we can. */
  short int tpid = cell_list_begin->tpid;
  if (tpid < 0) tpid = 0;
  cell_list_end->next = s->cells_sub[tpid];
  s->cells_sub[tpid] = cell_list_begin;
  atomic_sub(&s->tot_cells, count);

  /* Hook the multipoles into the buffer. */
  if (s->with_self_gravity) {
    multipole_list_end->next = s->multipoles_sub[tpid];
    s->multipoles_sub[tpid] = multipole_list_begin;
  }

  /* Unlock the space. */
  lock_unlock_blind(&s->lock);
}

/**
 * @brief Free up any allocated cells.
 *
 * @param s The #space.
 */
void space_free_cells(struct space *s) {

  ticks tic = getticks();

  threadpool_map(&s->e->threadpool, space_rebuild_recycle_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), threadpool_auto_chunk_size,
                 s);
  s->maxdepth = 0;

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
