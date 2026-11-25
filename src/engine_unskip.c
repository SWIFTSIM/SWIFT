/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "memswap.h"

/* Load the profiler header, if needed. */
#ifdef WITH_PROFILER
#include <gperftools/profiler.h>
#endif

/**
 * @brief Broad categories of tasks.
 *
 * Each category is unskipped independently
 * of the others.
 */
enum task_broad_types {
  task_broad_types_hydro = 1,
  task_broad_types_gravity,
  task_broad_types_stars,
  task_broad_types_sinks,
  task_broad_types_black_holes,
  task_broad_types_rt,
  task_broad_types_count,
};

/**
 * @brief Meta-data for the unskipping
 */
struct unskip_data {

  /*! The #engine */
  struct engine *e;

  /*! Pointer to the start of the list of cells to unskip */
  int *list_base;

  /*! Number of times the list has been duplicated */
  int multiplier;

  /*! The number of active cells (without dulication) */
  int num_active_cells;

  /*! The #task_broad_types corresponding to each copy of the list */
  enum task_broad_types task_types[task_broad_types_count];
};

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Check that the chemistry flux prep pass is active wherever the force
 * loop is.
 *
 * A force task waits for its prep task through the chemistry flux ghosts of
 * its local supers. The scheduler ignores a skipped dependency, so a skipped
 * prep task or ghost would let the force loop read an outflow that is not
 * complete, with no error.
 *
 * @param e The #engine.
 */
static void engine_check_chemistry_flux_prep_activation(
    const struct engine *e) {

  if (!chemistry_need_flux_prep_loop(e->chemistry)) return;

  const struct scheduler *s = &e->sched;
  for (int k = 0; k < s->nr_tasks; k++) {
    const struct task *t = &s->tasks[k];
    if (t->skip || t->subtype != task_subtype_force) continue;
    if (t->type != task_type_self && t->type != task_type_pair) continue;

    const struct task *prep = NULL;
    for (const struct link *l = t->ci->hydro.chemistry_flux_prep; l != NULL;
         l = l->next) {
      if (l->t->ci == t->ci && l->t->cj == t->cj) {
        prep = l->t;
        break;
      }
    }
    const struct cell *ci = t->ci;
    const struct cell *cj = t->cj != NULL ? t->cj : t->ci;
    if (prep == NULL || prep->skip)
      error(
          "An active %s force task has no active chemistry flux prep task "
          "(ci: loc=[%g %g %g] depth=%d node=%d; cj: loc=[%g %g %g] depth=%d "
          "node=%d).",
          taskID_names[t->type], ci->loc[0], ci->loc[1], ci->loc[2], ci->depth,
          ci->nodeID, cj->loc[0], cj->loc[1], cj->loc[2], cj->depth,
          cj->nodeID);

    const struct cell *cells[2] = {t->ci, t->cj};
    for (int n = 0; n < 2; n++) {
      const struct cell *c = cells[n];
      if (c == NULL || c->nodeID != e->nodeID) continue;
      const struct task *ghost = c->hydro.super->hydro.chemistry_flux_ghost;
      if (ghost != NULL && ghost->skip)
        error(
            "An active %s force task has a skipped chemistry flux ghost on its "
            "local cell (loc=[%g %g %g] depth=%d, super depth=%d).",
            taskID_names[t->type], c->loc[0], c->loc[1], c->loc[2], c->depth,
            c->hydro.super->depth);
    }
  }
}
#endif

#if defined(SWIFT_DEBUG_CHECKS) && defined(WITH_MPI) && \
    defined(EXTRA_HYDRO_LOOP)
/**
 * @brief Find the communication task of the chemistry flux prep theta in a
 * list of links.
 *
 * @param l The first link of the list.
 * @param type The task type (send, recv, pack or unpack).
 * @param nodeID The rank the task sends to, or -1 for a recv or an unpack.
 */
static const struct task *engine_find_flux_theta_task(
    const struct link *l, const enum task_types type, const int nodeID) {
  for (; l != NULL; l = l->next) {
    const struct task *t = l->t;
    if (t->type != type || t->subtype != task_subtype_chemistry_flux_prep)
      continue;
    if (nodeID < 0 || t->cj->nodeID == nodeID) return t;
  }
  return NULL;
}
#endif

/**
 * @brief Check that the theta messages of the chemistry flux prep loop are
 * active wherever the force loop needs them, on both sides.
 *
 * A force pair with a foreign cell reads the theta of the foreign donors.
 * The scheduler ignores a skipped dependency, so a skipped receive, unpack,
 * pack or chemistry flux ghost would let it read an old theta, with no error.
 * A send without its receive on the other rank would hang. This compares the
 * counts of the two ranks, so it must be called on all ranks.
 *
 * @param e The #engine.
 */
void engine_check_chemistry_flux_theta_activation(const struct engine *e) {

#if defined(SWIFT_DEBUG_CHECKS) && defined(WITH_MPI) && \
    defined(EXTRA_HYDRO_LOOP)
  if (!(e->policy & engine_policy_hydro)) return;
  if (!chemistry_need_flux_prep_loop(e->chemistry)) return;

  const struct scheduler *s = &e->sched;
  const int nr_nodes = e->nr_nodes;
  int *n_send = (int *)calloc(nr_nodes, sizeof(int));
  int *n_recv = (int *)calloc(nr_nodes, sizeof(int));
  int *n_sent_here = (int *)calloc(nr_nodes, sizeof(int));
  if (n_send == NULL || n_recv == NULL || n_sent_here == NULL)
    error("Failed to allocate the theta message counts.");

  for (int k = 0; k < s->nr_tasks; k++) {
    const struct task *t = &s->tasks[k];
    if (t->skip) continue;

    if (t->subtype == task_subtype_chemistry_flux_prep &&
        t->type == task_type_send) {
      const struct cell *c = t->ci;
      n_send[t->cj->nodeID]++;
      const struct task *pack = engine_find_flux_theta_task(
          c->mpi.pack, task_type_pack, t->cj->nodeID);
      const struct task *ghost = c->hydro.super->hydro.chemistry_flux_ghost;
      if (pack == NULL || pack->skip || ghost == NULL || ghost->skip)
        error(
            "An active theta send to rank %d has a skipped pack or chemistry "
            "flux ghost (cell loc=[%g %g %g] depth=%d).",
            t->cj->nodeID, c->loc[0], c->loc[1], c->loc[2], c->depth);

    } else if (t->subtype == task_subtype_chemistry_flux_prep &&
               t->type == task_type_recv) {
      n_recv[t->ci->nodeID]++;
      if (t->nr_unlock_tasks != 1 || t->unlock_tasks[0]->skip)
        error(
            "An active theta recv from rank %d has no active unpack (cell "
            "loc=[%g %g %g] depth=%d).",
            t->ci->nodeID, t->ci->loc[0], t->ci->loc[1], t->ci->loc[2],
            t->ci->depth);

    } else if (t->subtype == task_subtype_force && t->type == task_type_pair &&
               (t->ci->nodeID != e->nodeID || t->cj->nodeID != e->nodeID)) {
      const int ci_local = t->ci->nodeID == e->nodeID;
      const struct cell *local = ci_local ? t->ci : t->cj;
      const struct cell *foreign = ci_local ? t->cj : t->ci;
      const struct task *send = engine_find_flux_theta_task(
          local->mpi.send, task_type_send, foreign->nodeID);
      const struct task *recv =
          engine_find_flux_theta_task(foreign->mpi.recv, task_type_recv, -1);
      if (send == NULL || send->skip || recv == NULL || recv->skip)
        error(
            "An active cross-rank force pair has a skipped theta send or recv "
            "(local cell loc=[%g %g %g] depth=%d; foreign cell loc=[%g %g %g] "
            "depth=%d node=%d).",
            local->loc[0], local->loc[1], local->loc[2], local->depth,
            foreign->loc[0], foreign->loc[1], foreign->loc[2], foreign->depth,
            foreign->nodeID);
    }
  }

  /* Each rank tells every other one how many theta messages it sends it. */
  MPI_Alltoall(n_send, 1, MPI_INT, n_sent_here, 1, MPI_INT, MPI_COMM_WORLD);
  for (int n = 0; n < nr_nodes; n++)
    if (n_sent_here[n] != n_recv[n])
      error(
          "Rank %d activated %d theta sends to rank %d, which activated %d "
          "theta recvs from it.",
          n, n_sent_here[n], e->nodeID, n_recv[n]);

  free(n_send);
  free(n_recv);
  free(n_sent_here);
#endif
}

/**
 * @brief Unskip any hydro tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_hydro(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->hydro.count == 0) return;

#ifndef MPI_SYMMETRIC_FORCE_INTERACTION
  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;
#endif

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_hydro(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_hydro_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any stars tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param with_star_formation Are we running with star formation switched on?
 */
static void engine_do_unskip_stars(struct cell *c, struct engine *e,
                                   const int with_star_formation,
                                   const int with_star_formation_sink) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  const int non_empty =
      c->stars.count > 0 || (with_star_formation && c->hydro.count > 0) ||
      (with_star_formation_sink && (c->hydro.count > 0 || c->sinks.count > 0));

  /* Ignore empty cells. */
  if (!non_empty) return;

  const int ci_active = cell_need_activating_stars(c, e, with_star_formation,
                                                   with_star_formation_sink);

  /* Skip inactive cells. */
  if (!ci_active) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_stars(cp, e, with_star_formation,
                               with_star_formation_sink);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_stars_tasks(
      c, &e->sched, with_star_formation, with_star_formation_sink);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any black hole tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_black_holes(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_black_holes(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_black_holes_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any sink tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_sinks(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->sinks.count == 0 && c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_sinks(c, e) && !cell_is_active_hydro(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_sinks(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_sinks_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any gravity tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_gravity(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->grav.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* At the top level we need to recursively check particles haven't moved
   * too far for the mesh gravity (if using the mesh). */
  if (c->depth == 0 && e->s->periodic) {
    cell_check_grav_mesh_pairs(c, e);
  }

  /* Recurse */
  if (c->split && ((c->maxdepth - c->depth) >= space_subdepth_diff_grav)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_gravity(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief Unskip any radiative transfer tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param sub_cycle 1 if this is a call for a sub cycle, 0 otherwise
 */
static void engine_do_unskip_rt(struct cell *c, struct engine *e,
                                const int sub_cycle) {

  /* Note: we only get this far if engine_policy_rt is flagged. */
#ifdef SWIFT_DEBUG_CHECKS
  if (!(e->policy & engine_policy_rt))
    error("Unksipping RT stuff without the policy being on");
#endif

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Do we have work to do? */
  if (c->hydro.count == 0) return;
  if (!cell_is_rt_active(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        engine_do_unskip_rt(c->progeny[k], e, sub_cycle);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_rt_tasks(c, &e->sched, sub_cycle);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Mapper function to unskip active tasks.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an unskip_data structure.
 */
void engine_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  /* Unpack the meta data */
  struct unskip_data *data = (struct unskip_data *)extra_data;
  const int num_active_cells = data->num_active_cells;
  const enum task_broad_types *const task_types = data->task_types;
  const int *const list_base = data->list_base;
  struct engine *e = data->e;
  struct cell *const cells_top = e->s->cells_top;

  /* What policies are we running? */
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_star_formation = e->policy & engine_policy_star_formation;
  const int with_star_formation_sink = with_sinks && with_stars;

  /* The current chunk of active cells */
  const int *const local_cells = (int *)map_data;

  /* Loop over this thread's chunk of cells to unskip */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Handle on the cell */
    struct cell *const c = &cells_top[local_cells[ind]];

    /* In what copy of the global list are we?
     * This gives us the broad type of task we are working on. */
    const ptrdiff_t delta = &local_cells[ind] - list_base;
    const int type = delta / num_active_cells;

#ifdef SWIFT_DEBUG_CHECKS
    if (type >= data->multiplier) error("Invalid broad task type!");
    if (c == NULL) error("Got an invalid cell index!");
#endif

    /* What broad type of tasks are we unskipping? */
    switch (task_types[type]) {
      case task_broad_types_hydro:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_hydro))
          error("Trying to unskip hydro tasks in a non-hydro run!");
#endif
        engine_do_unskip_hydro(c, e);
        break;
      case task_broad_types_gravity:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_self_gravity) &&
            !(e->policy & engine_policy_external_gravity))
          error("Trying to unskip gravity tasks in a non-gravity run!");
#endif
        engine_do_unskip_gravity(c, e);
        break;
      case task_broad_types_stars:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_stars))
          error("Trying to unskip star tasks in a non-stars run!");
#endif
        engine_do_unskip_stars(c, e, with_star_formation,
                               with_star_formation_sink);
        break;
      case task_broad_types_sinks:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_sinks))
          error("Trying to unskip sink tasks in a non-sinks run!");
#endif
        engine_do_unskip_sinks(c, e);
        break;
      case task_broad_types_black_holes:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_black_holes))
          error("Trying to unskip black holes tasks in a non-BH run!");
#endif
        engine_do_unskip_black_holes(c, e);
        break;
      case task_broad_types_rt:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_rt))
          error("Trying to unskip radiative transfer tasks in a non-rt run!");
#endif
        engine_do_unskip_rt(c, e, /*sub_cycle=*/0);
        break;
      default:
#ifdef SWIFT_DEBUG_CHECKS
        error("Invalid broad task type!");
#endif
        continue;
    }
  }
}

/**
 * @brief Unskip all the tasks that act on active cells at this time.
 *
 * @param e The #engine.
 */
void engine_unskip(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_self_grav = e->policy & engine_policy_self_gravity;
  const int with_ext_grav = e->policy & engine_policy_external_gravity;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_black_holes = e->policy & engine_policy_black_holes;
  const int with_rt = e->policy & engine_policy_rt;

#ifdef WITH_PROFILER
  static int count = 0;
  char filename[100];
  sprintf(filename, "/tmp/swift_engine_do_usnkip_mapper_%06i.prof", count++);
  ProfilerStart(filename);
#endif  // WITH_PROFILER

  /* Move the active local cells to the top of the list. */
  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if (cell_is_empty(c)) continue;

    if ((with_hydro && cell_is_active_hydro(c, e)) ||
        (with_self_grav && cell_is_active_gravity(c, e)) ||
        (with_ext_grav && c->nodeID == nodeID &&
         cell_is_active_gravity(c, e)) ||
        (with_feedback && cell_is_active_stars(c, e)) ||
        (with_stars && c->nodeID == nodeID && cell_is_active_stars(c, e)) ||
        (with_sinks && cell_is_active_sinks(c, e)) ||
        (with_black_holes && cell_is_active_black_holes(c, e)) ||
        (with_rt && cell_is_rt_active(c, e))) {

      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

  /* What kind of tasks do we have? */
  struct unskip_data data;
  bzero(&data, sizeof(struct unskip_data));
  int multiplier = 0;
  if (with_hydro) {
    data.task_types[multiplier] = task_broad_types_hydro;
    multiplier++;
  }
  if (with_self_grav || with_ext_grav) {
    data.task_types[multiplier] = task_broad_types_gravity;
    multiplier++;
  }
  if (with_feedback || with_stars) {
    data.task_types[multiplier] = task_broad_types_stars;
    multiplier++;
  }
  if (with_sinks) {
    data.task_types[multiplier] = task_broad_types_sinks;
    multiplier++;
  }
  if (with_black_holes) {
    data.task_types[multiplier] = task_broad_types_black_holes;
    multiplier++;
  }
  if (with_rt) {
    data.task_types[multiplier] = task_broad_types_rt;
    multiplier++;
  }

  /* Should we duplicate the list of active cells to better parallelise the
     unskip over the threads ? */
  int *local_active_cells;
  if (multiplier > 1) {

    /* Make space for copies of the list */
    local_active_cells =
        (int *)malloc(multiplier * num_active_cells * sizeof(int));
    if (local_active_cells == NULL)
      error(
          "Couldn't allocate memory for duplicated list of local active "
          "cells.");

    /* Make blind copies of the list */
    for (int m = 0; m < multiplier; m++) {
      memcpy(local_active_cells + m * num_active_cells, local_cells,
             num_active_cells * sizeof(int));
    }
  } else {
    local_active_cells = local_cells;
  }

  /* We now have a list of local active cells duplicated as many times as
   * we have broad task types. We can now release all the threads on the list */

  data.e = e;
  data.list_base = local_active_cells;
  data.num_active_cells = num_active_cells;
  data.multiplier = multiplier;

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, engine_do_unskip_mapper, local_active_cells,
                 num_active_cells * multiplier, sizeof(int), /*chunk=*/1,
                 &data);

#ifdef SWIFT_DEBUG_CHECKS
  if (with_hydro) engine_check_chemistry_flux_prep_activation(e);
#endif

#ifdef WITH_PROFILER
  ProfilerStop();
#endif  // WITH_PROFILER

  /* Free stuff? */
  if (multiplier > 1) {
    free(local_active_cells);
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void engine_do_unskip_sub_cycle_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *const cells_top = e->s->cells_top;

  /* The current chunk of active cells */
  const int *const local_cells = (int *)map_data;

  /* Loop over this thread's chunk of cells to unskip */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Handle on the cell */
    struct cell *const c = &cells_top[local_cells[ind]];

    engine_do_unskip_rt(c, e, /*sub_cycle=*/1);
  }
}

/**
 * @brief Unskip all the RT tasks that are active during this sub-cycle.
 *
 * @param e The #engine.
 */
void engine_unskip_rt_sub_cycle(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;
  const int with_rt = e->policy & engine_policy_rt;

  if (!with_rt) error("Unskipping sub-cycles when running without RT!");

  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if (c->hydro.count == 0) continue;

    if (cell_is_rt_active(c, e)) {

      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, engine_do_unskip_sub_cycle_mapper, local_cells,
                 num_active_cells, sizeof(int), /*chunk=*/1, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
