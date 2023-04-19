/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "engine.h"
#include "lightcone/lightcone_array.h"

/**
 * @brief Mapper function to drift *all* the #part to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_part_mapper(void *map_data, int num_elements,
                                     void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Drift all the particles */
      cell_drift_part(c, e, /* force the drift=*/1, NULL);
    }
  }
}

/**
 * @brief Mapper function to drift *all* the #gpart to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_gpart_mapper(void *map_data, int num_elements,
                                      void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Drift all the particles */
      cell_drift_gpart(c, e, /* force the drift=*/1, /*replication_list=*/NULL);
    }
  }
}

/**
 * @brief Mapper function to drift *all* the #spart to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_spart_mapper(void *map_data, int num_elements,
                                      void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Drift all the particles */
      cell_drift_spart(c, e, /* force the drift=*/1, NULL);
    }
  }
}

/**
 * @brief Mapper function to drift *all* the #bpart to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_bpart_mapper(void *map_data, int num_elements,
                                      void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Drift all the particles */
      cell_drift_bpart(c, e, /* force the drift=*/1, NULL);
    }
  }
}

/**
 * @brief Mapper function to drift *all* the #sink to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_sink_mapper(void *map_data, int num_elements,
                                     void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_top[ind]];

    if (c->nodeID == e->nodeID) {

      /* Drift all the particles */
      cell_drift_sink(c, e, /* force the drift=*/1);
    }
  }
}

/**
 * @brief Mapper function to drift *all* the multipoles to the current time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_multipole_mapper(void *map_data, int num_elements,
                                          void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;
  const int restarting = e->restarting;
  struct space *s = e->s;
  struct cell *cells_top;
  int *local_cells_with_tasks_top;

  if (restarting) {

    /* When restarting, we loop over all top-level cells */
    cells_top = (struct cell *)map_data;
    local_cells_with_tasks_top = NULL;

  } else {

    /* In any other case, we use the list of local cells with tasks */
    cells_top = s->cells_top;
    local_cells_with_tasks_top = (int *)map_data;
  }

  for (int ind = 0; ind < num_elements; ind++) {

    struct cell *c;

    /* When restarting, the list of local cells does not
       yet exist. We use the raw list of top-level cells instead */
    if (restarting)
      c = &cells_top[ind];
    else
      c = &cells_top[local_cells_with_tasks_top[ind]];

    cell_drift_all_multipoles(c, e);
  }
}

/**
 * @brief Drift *all* particles and multipoles at all levels
 * forward to the current time.
 *
 * @param e The #engine.
 * @param drift_mpoles Do we want to drift all the multipoles as well?
 */
void engine_drift_all(struct engine *e, const int drift_mpoles) {

  const ticks tic = getticks();

  if (e->nodeID == 0 && e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Drifting all to a=%15.12e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Drifting all to t=%15.12e",
              e->ti_current * e->time_base + e->time_begin);
  }

#ifdef WITH_LIGHTCONE
  /* Determine which periodic replications could contribute to the lightcone
     during this time step */
  lightcone_array_prepare_for_step(e->lightcone_array_properties, e->cosmology,
                                   e->ti_earliest_undrifted, e->ti_current);
#endif

  if (!e->restarting) {

    /* Normal case: We have a list of local cells with tasks to play with */

    if (e->s->nr_parts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_part_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_gparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_gpart_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_sparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_spart_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_sinks > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_sink_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_bparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_bpart_mapper,
                     e->s->local_cells_top, e->s->nr_local_cells, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }
    if (drift_mpoles && (e->policy & engine_policy_self_gravity)) {
      threadpool_map(&e->threadpool, engine_do_drift_all_multipole_mapper,
                     e->s->local_cells_with_tasks_top,
                     e->s->nr_local_cells_with_tasks, sizeof(int),
                     threadpool_auto_chunk_size, e);
    }

  } else {

    /* When restarting, the list of local cells with tasks does not yet
       exist. We use the raw list of top-level cells instead */

    if (e->s->nr_parts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_part_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_sparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_spart_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_sinks > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_sink_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_bparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_bpart_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
    if (e->s->nr_gparts > 0) {
      threadpool_map(&e->threadpool, engine_do_drift_all_gpart_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
    if (e->policy & engine_policy_self_gravity) {
      threadpool_map(&e->threadpool, engine_do_drift_all_multipole_mapper,
                     e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                     threadpool_auto_chunk_size, e);
    }
  }

  /* Synchronize particle positions */
  space_synchronize_particle_positions(e->s);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time. */
  space_check_drift_point(
      e->s, e->ti_current,
      drift_mpoles && (e->policy & engine_policy_self_gravity));
  part_verify_links(e->s->parts, e->s->gparts, e->s->sinks, e->s->sparts,
                    e->s->bparts, e->s->nr_parts, e->s->nr_gparts,
                    e->s->nr_sinks, e->s->nr_sparts, e->s->nr_bparts,
                    e->verbose);
#endif

  /* All particles have now been drifted to ti_current */
  e->ti_earliest_undrifted = e->ti_current;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#ifdef WITH_LIGHTCONE
  /* Drifting all of the particles can cause many particles to cross
     the lightcone, so flush buffers now to reduce peak memory use . */
  lightcone_array_flush(e->lightcone_array_properties, &e->threadpool,
                        e->cosmology, e->internal_units, e->snapshot_units,
                        /*flush_map_updates=*/1, /*flush_particles=*/1,
                        /*end_file=*/0, /*dump_all_shells=*/0);
#endif
}

/**
 * @brief Mapper function to drift *all* top-level multipoles forward in
 * time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_top_multipoles_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &cells[ind];
    if (c != NULL) {

      /* Drift the multipole at this level only */
      if (c->grav.ti_old_multipole != e->ti_current) cell_drift_multipole(c, e);
    }
  }
}

/**
 * @brief Drift *all* top-level multipoles forward to the current time.
 *
 * @param e The #engine.
 */
void engine_drift_top_multipoles(struct engine *e) {

  const ticks tic = getticks();

  threadpool_map(&e->threadpool, engine_do_drift_top_multipoles_mapper,
                 e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                 threadpool_auto_chunk_size, e);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time. */
  space_check_top_multipoles_drift_point(e->s, e->ti_current);
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
