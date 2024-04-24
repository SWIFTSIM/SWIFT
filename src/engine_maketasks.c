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

/* Some standard headers. */
#include <stdlib.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Load the profiler header, if needed. */
#ifdef WITH_PROFILER
#include <gperftools/profiler.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "adaptive_softening.h"
#include "atomic.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "feedback.h"
#include "neutrino_properties.h"
#include "proxy.h"
#include "rt_properties.h"
#include "timers.h"

extern int engine_max_parts_per_ghost;
extern int engine_max_sparts_per_ghost;
extern int engine_star_resort_task_depth;
extern int engine_max_parts_per_cooling;

/**
 * @brief Add send tasks for the gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_grav The send_grav #task, if it has already been created.
 */
void engine_addtasks_send_gravity(struct engine *e, struct cell *ci,
                                  struct cell *cj, struct task *t_grav) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(ci, cell_flag_has_tasks)) return;

  /* Check if any of the gravity tasks are for the target node. */
  for (l = ci->grav.grav; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_grav == NULL) {

      /* Make sure this cell is tagged. */
      cell_ensure_tagged(ci);

      t_grav = scheduler_addtask(s, task_type_send, task_subtype_gpart,
                                 ci->mpi.tag, 0, ci, cj);

      /* The sends should unlock the down pass. */
      scheduler_addunlock(s, t_grav, ci->grav.super->grav.down);

      /* Drift before you send */
      scheduler_addunlock(s, ci->grav.super->grav.drift, t_grav);

      if (gravity_after_hydro_density)
        scheduler_addunlock(s, ci->grav.super->grav.init_out, t_grav);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.send, t_grav);
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_gravity(e, ci->progeny[k], cj, t_grav);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the hydro pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_xv The send_xv #task, if it has already been created.
 * @param t_rho The send_rho #task, if it has already been created.
 * @param t_gradient The send_gradient #task, if already created.
 * @param t_prep1 The send_prep1 #task, if it has already been created.
 * @param t_limiter The send_limiter #task, if it has already been created.
 * @param t_rt_gradient The send_rt_gradient #task, if it has already been
 * created.
 * @param t_rt_transport The send_rt_transport #task, if it has already been
 * @param with_feedback Are we running with stellar feedback?
 * @param with_limiter Are we running with the time-step limiter?
 * @param with_sync Are we running with time-step synchronization?
 * @param with_rt Are we running with radiative transfer?
 */
void engine_addtasks_send_hydro(struct engine *e, struct cell *ci,
                                struct cell *cj, struct task *t_xv,
                                struct task *t_rho, struct task *t_gradient,
                                struct task *t_prep1, struct task *t_limiter,
                                struct task *t_pack_limiter,
                                struct task *t_rt_gradient,
                                struct task *t_rt_transport,
                                const int with_feedback, const int with_limiter,
                                const int with_sync, const int with_rt) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(ci, cell_flag_has_tasks)) return;

  /* Check if any of the density tasks are for the target node. */
  for (l = ci->hydro.density; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_xv == NULL) {

      /* Make sure this cell is tagged. */
      cell_ensure_tagged(ci);

      t_xv = scheduler_addtask(s, task_type_send, task_subtype_xv, ci->mpi.tag,
                               0, ci, cj);
      t_rho = scheduler_addtask(s, task_type_send, task_subtype_rho,
                                ci->mpi.tag, 0, ci, cj);
      scheduler_addunlock(s, t_xv, t_rho);

#ifdef EXTRA_HYDRO_LOOP
      t_gradient = scheduler_addtask(s, task_type_send, task_subtype_gradient,
                                     ci->mpi.tag, 0, ci, cj);
      scheduler_addunlock(s, t_rho, t_gradient);
#endif

      if (with_limiter) {
        t_limiter = scheduler_addtask(s, task_type_send, task_subtype_limiter,
                                      ci->mpi.tag, 0, ci, cj);
        t_pack_limiter = scheduler_addtask(s, task_type_pack,
                                           task_subtype_limiter, 0, 0, ci, cj);

        scheduler_addunlock(s, t_pack_limiter, t_limiter);
      }

#ifdef EXTRA_STAR_LOOPS
      if (with_feedback) {
        t_prep1 = scheduler_addtask(s, task_type_send, task_subtype_part_prep1,
                                    ci->mpi.tag, 0, ci, cj);
      }
#endif

      if (with_rt) {
        /* Add the RT sends */
        t_rt_gradient =
            scheduler_addtask(s, task_type_send, task_subtype_rt_gradient,
                              ci->mpi.tag, 0, ci, cj);

        t_rt_transport =
            scheduler_addtask(s, task_type_send, task_subtype_rt_transport,
                              ci->mpi.tag, 0, ci, cj);
      }

#ifdef EXTRA_HYDRO_LOOP

      scheduler_addunlock(s, t_gradient, ci->hydro.super->hydro.end_force);

      scheduler_addunlock(s, ci->hydro.super->hydro.extra_ghost, t_gradient);

      /* The send_rho task should unlock the super_hydro-cell's extra_ghost
       * task. */
      scheduler_addunlock(s, t_rho, ci->hydro.super->hydro.extra_ghost);

      /* The send_rho task depends on the cell's ghost task. */
      scheduler_addunlock(s, ci->hydro.super->hydro.ghost_out, t_rho);

      /* The send_xv task should unlock the super_hydro-cell's ghost task. */
      scheduler_addunlock(s, t_xv, ci->hydro.super->hydro.ghost_in);

#else
      /* The send_rho task should unlock the super_hydro-cell's kick task. */
      scheduler_addunlock(s, t_rho, ci->hydro.super->hydro.end_force);

      /* The send_rho task depends on the cell's ghost task. */
      scheduler_addunlock(s, ci->hydro.super->hydro.ghost_out, t_rho);

      /* The send_xv task should unlock the super_hydro-cell's ghost task. */
      scheduler_addunlock(s, t_xv, ci->hydro.super->hydro.ghost_in);

#endif

      scheduler_addunlock(s, ci->hydro.super->hydro.drift, t_rho);

      /* Drift before you send */
      scheduler_addunlock(s, ci->hydro.super->hydro.drift, t_xv);

      if (with_limiter)
        scheduler_addunlock(s, ci->super->timestep, t_pack_limiter);

#ifdef EXTRA_STAR_LOOPS
      /* In stellar feedback, send gas parts only after they have finished their
       * hydro ghosts */
      if (with_feedback) {
        scheduler_addunlock(s, ci->hydro.super->hydro.prep1_ghost, t_prep1);
        scheduler_addunlock(s, t_prep1, ci->hydro.super->stars.prep2_ghost);
      }
#endif

      if (with_rt) {
        /* Don't send the transport stuff before the gradient stuff */
        scheduler_addunlock(s, t_rt_gradient, t_rt_transport);

        /* The send_gradient task depends on the cell's ghost1 task. */
        scheduler_addunlock(s, ci->hydro.super->rt.rt_ghost1, t_rt_gradient);

        /* The send_transport task depends on the cell's ghost2 task. */
        scheduler_addunlock(s, ci->hydro.super->rt.rt_ghost2, t_rt_transport);

        /* Safety measure: collect dependencies and make sure data is sent
         * before modifying it */
        scheduler_addunlock(s, t_rt_gradient, ci->hydro.super->rt.rt_ghost2);

        /* Safety measure: collect dependencies and make sure data is sent
         * before modifying it */
        scheduler_addunlock(s, t_rt_transport,
                            ci->hydro.super->rt.rt_transport_out);

        /* Drift before you send. Especially intended to cover inactive cells
         * being sent. */
        scheduler_addunlock(s, ci->hydro.super->hydro.drift, t_rt_gradient);
        scheduler_addunlock(s, ci->hydro.super->hydro.drift, t_rt_transport);

        /* Make sure the gradient sends don't run before the xv is finished.
         * This can occur when a cell itself is inactive for both hydro and
         * RT, but needs to be sent over for some other cell's pair task.
         * The rt_gradient - xv dependency is special because when received,
         * these two tasks will/may activate the sorts.*/
        scheduler_addunlock(s, t_xv, t_rt_gradient);
      }
    } /* if t_xv == NULL */

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.send, t_xv);
    engine_addlink(e, &ci->mpi.send, t_rho);
#ifdef EXTRA_HYDRO_LOOP
    engine_addlink(e, &ci->mpi.send, t_gradient);
#endif
    if (with_limiter) {
      engine_addlink(e, &ci->mpi.send, t_limiter);
      engine_addlink(e, &ci->mpi.pack, t_pack_limiter);
    }
#ifdef EXTRA_STAR_LOOPS
    if (with_feedback) engine_addlink(e, &ci->mpi.send, t_prep1);
#endif

    if (with_rt) {
      /* Add them to the local cell. */
      engine_addlink(e, &ci->mpi.send, t_rt_gradient);
      engine_addlink(e, &ci->mpi.send, t_rt_transport);
    }
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_hydro(
            e, ci->progeny[k], cj, t_xv, t_rho, t_gradient, t_prep1, t_limiter,
            t_pack_limiter, t_rt_gradient, t_rt_transport, with_feedback,
            with_limiter, with_sync, with_rt);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the stars pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_density The send_density #task, if it has already been created.
 * @param t_prep2 The send_prep2 #task, if it has already been created.
 * @param t_sf_counts The send_sf_counts, if it has been created.
 * @param with_star_formation Are we running with star formation on?
 */
void engine_addtasks_send_stars(struct engine *e, struct cell *ci,
                                struct cell *cj, struct task *t_density,
                                struct task *t_prep2, struct task *t_sf_counts,
                                const int with_star_formation) {
#ifdef SWIFT_DEBUG_CHECKS
  if (e->policy & engine_policy_sinks && e->policy & engine_policy_stars) {
    error("TODO");
  }
#endif

#ifdef WITH_MPI

  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(ci, cell_flag_has_tasks)) return;

  if (t_sf_counts == NULL && with_star_formation && ci->hydro.count > 0) {
#ifdef SWIFT_DEBUG_CHECKS
    if (ci->depth != 0)
      error(
          "Attaching a sf_count task at a non-top level c->depth=%d "
          "c->count=%d",
          ci->depth, ci->hydro.count);
#endif
    t_sf_counts = scheduler_addtask(s, task_type_send, task_subtype_sf_counts,
                                    ci->mpi.tag, 0, ci, cj);
    scheduler_addunlock(s, ci->hydro.star_formation, t_sf_counts);
  }

  /* Check if any of the density tasks are for the target node. */
  for (l = ci->stars.density; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    if (t_density == NULL) {

      /* Make sure this cell is tagged. */
      cell_ensure_tagged(ci);

      /* Create the tasks and their dependencies? */
      t_density =
          scheduler_addtask(s, task_type_send, task_subtype_spart_density,
                            ci->mpi.tag, 0, ci, cj);

#ifdef EXTRA_STAR_LOOPS
      t_prep2 = scheduler_addtask(s, task_type_send, task_subtype_spart_prep2,
                                  ci->mpi.tag, 0, ci, cj);
#endif

#ifdef EXTRA_STAR_LOOPS
      /* The first send_stars task should unlock prep1 ghost */
      scheduler_addunlock(s, t_density, ci->hydro.super->stars.prep1_ghost);

      /* Prep2 ghost before second send */
      scheduler_addunlock(s, ci->hydro.super->stars.prep2_ghost, t_prep2);

      /* The second send_stars task should unlock the super_cell's "end of star
       * block" task. */
      scheduler_addunlock(s, t_prep2, ci->hydro.super->stars.stars_out);
#else
      /* The send_stars task should unlock the super_cell's "end of star block"
       * task. */
      scheduler_addunlock(s, t_density, ci->hydro.super->stars.stars_out);
#endif

      /* Density ghost before first send */
      scheduler_addunlock(s, ci->hydro.super->stars.density_ghost, t_density);

      /* Drift before first send */
      scheduler_addunlock(s, ci->hydro.super->stars.drift, t_density);

      if (with_star_formation && ci->hydro.count > 0) {
        scheduler_addunlock(s, t_sf_counts, t_density);
#ifdef EXTRA_STAR_LOOPS
        scheduler_addunlock(s, t_sf_counts, t_prep2);
#endif
      }
    }

    engine_addlink(e, &ci->mpi.send, t_density);
#ifdef EXTRA_STAR_LOOPS
    engine_addlink(e, &ci->mpi.send, t_prep2);
#endif
    if (with_star_formation && ci->hydro.count > 0) {
      engine_addlink(e, &ci->mpi.send, t_sf_counts);
    }
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_stars(e, ci->progeny[k], cj, t_density, t_prep2,
                                   t_sf_counts, with_star_formation);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the black holes pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_rho The density comm. task, if it has already been created.
 * @param t_bh_merger The BH swallow comm. task, if it has already been created.
 * @param t_gas_swallow The gas swallow comm. task, if it has already been
 * created.
 * @param t_feedback The send_feed #task, if it has already been created.
 */
void engine_addtasks_send_black_holes(struct engine *e, struct cell *ci,
                                      struct cell *cj, struct task *t_rho,
                                      struct task *t_bh_merger,
                                      struct task *t_gas_swallow,
                                      struct task *t_feedback) {

#ifdef WITH_MPI

  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(ci, cell_flag_has_tasks)) return;

  /* Check if any of the density tasks are for the target node. */
  for (l = ci->black_holes.density; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    if (t_rho == NULL) {

      /* Make sure this cell is tagged. */
      cell_ensure_tagged(ci);

      /* Create the tasks and their dependencies? */
      t_rho = scheduler_addtask(s, task_type_send, task_subtype_bpart_rho,
                                ci->mpi.tag, 0, ci, cj);

      t_bh_merger = scheduler_addtask(
          s, task_type_send, task_subtype_bpart_merger, ci->mpi.tag, 0, ci, cj);

      t_gas_swallow = scheduler_addtask(
          s, task_type_send, task_subtype_part_swallow, ci->mpi.tag, 0, ci, cj);

      t_feedback =
          scheduler_addtask(s, task_type_send, task_subtype_bpart_feedback,
                            ci->mpi.tag, 0, ci, cj);

      /* The send_black_holes task should unlock the super_cell's BH exit point
       * task. */
      scheduler_addunlock(s, t_feedback,
                          ci->hydro.super->black_holes.black_holes_out);

      scheduler_addunlock(s, ci->hydro.super->black_holes.swallow_ghost_3,
                          t_feedback);

      /* Ghost before you send */
      scheduler_addunlock(s, ci->hydro.super->black_holes.drift, t_rho);
      scheduler_addunlock(s, ci->hydro.super->black_holes.density_ghost, t_rho);
      scheduler_addunlock(s, t_rho,
                          ci->hydro.super->black_holes.swallow_ghost_1);

      scheduler_addunlock(s, ci->hydro.super->black_holes.swallow_ghost_1,
                          t_bh_merger);
      scheduler_addunlock(s, t_bh_merger,
                          ci->hydro.super->black_holes.swallow_ghost_3);

      scheduler_addunlock(s, ci->hydro.super->black_holes.swallow_ghost_1,
                          t_gas_swallow);
      scheduler_addunlock(s, t_gas_swallow,
                          ci->hydro.super->black_holes.swallow_ghost_2);
    }

    engine_addlink(e, &ci->mpi.send, t_rho);
    engine_addlink(e, &ci->mpi.send, t_bh_merger);
    engine_addlink(e, &ci->mpi.send, t_gas_swallow);
    engine_addlink(e, &ci->mpi.send, t_feedback);
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_black_holes(e, ci->progeny[k], cj, t_rho,
                                         t_bh_merger, t_gas_swallow,
                                         t_feedback);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for hydro pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_xv The recv_xv #task, if it has already been created.
 * @param t_rho The recv_rho #task, if it has already been created.
 * @param t_gradient The recv_gradient #task, if it has already been created.
 * @param t_prep1 The recv_prep1 #task, if it has already been created.
 * @param t_limiter The recv_limiter #task, if it has already been created.
 * @param t_unpack_limiter The unpack_limiter #task, if it has already been
 * created.
 * @param t_rt_gradient The recv_rt_gradient #task, if it has already been
 * created.
 * @param t_rt_transport The recv_rt_transport #task, if it has already been
 * created.
 * @param t_rt_sorts The rt_sort #task, if it has already been created.
 * @param tend The top-level time-step communication #task.
 * @param with_feedback Are we running with stellar feedback?
 * @param with_black_holes Are we running with black holes?
 * @param with_limiter Are we running with the time-step limiter?
 * @param with_sync Are we running with time-step synchronization?
 * @param with_rt Are we running with radiative transfer?
 */
void engine_addtasks_recv_hydro(
    struct engine *e, struct cell *c, struct task *t_xv, struct task *t_rho,
    struct task *t_gradient, struct task *t_prep1, struct task *t_limiter,
    struct task *t_unpack_limiter, struct task *t_rt_gradient,
    struct task *t_rt_transport, struct task *t_rt_sorts,
    struct task *const tend, const int with_feedback,
    const int with_black_holes, const int with_limiter, const int with_sync,
    const int with_rt) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Have we reached a level where there are any hydro tasks ? */
  if (t_xv == NULL && c->hydro.density != NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif /* SWIFT_DEBUG_CHECKS */

    /* Create the tasks. */
    t_xv = scheduler_addtask(s, task_type_recv, task_subtype_xv, c->mpi.tag, 0,
                             c, NULL);
    t_rho = scheduler_addtask(s, task_type_recv, task_subtype_rho, c->mpi.tag,
                              0, c, NULL);

    scheduler_addunlock(s, t_xv, t_rho);

#ifdef EXTRA_HYDRO_LOOP
    t_gradient = scheduler_addtask(s, task_type_recv, task_subtype_gradient,
                                   c->mpi.tag, 0, c, NULL);
    scheduler_addunlock(s, t_xv, t_gradient);
    scheduler_addunlock(s, t_rho, t_gradient);
#endif

    if (with_limiter) {
      t_limiter = scheduler_addtask(s, task_type_recv, task_subtype_limiter,
                                    c->mpi.tag, 0, c, NULL);
      t_unpack_limiter = scheduler_addtask(s, task_type_unpack,
                                           task_subtype_limiter, 0, 0, c, NULL);

      scheduler_addunlock(s, t_limiter, t_unpack_limiter);
    }

#ifdef EXTRA_STAR_LOOPS
    if (with_feedback) {
      t_prep1 = scheduler_addtask(s, task_type_recv, task_subtype_part_prep1,
                                  c->mpi.tag, 0, c, NULL);
    }
#endif

    if (with_rt) {
      /* Create the tasks. */
      t_rt_gradient = scheduler_addtask(
          s, task_type_recv, task_subtype_rt_gradient, c->mpi.tag, 0, c, NULL);
      t_rt_transport = scheduler_addtask(
          s, task_type_recv, task_subtype_rt_transport, c->mpi.tag, 0, c, NULL);
      /* Also create the rt_advance_cell_time tasks for the foreign cells
       * for the sub-cycling. */

#ifdef SWIFT_RT_DEBUG_CHECKS
      if (c->super == NULL)
        error("trying to add rt_advance_cell_time above super level...");
      if (c->top->rt.rt_collect_times == NULL) {
        error("rt_collect_times should exist already");
      }
      if (c->super->rt.rt_advance_cell_time == NULL) {
        error("rt_advance_cell_times should exist already");
      }
#endif

      /* Make sure we sort after receiving RT data. The hydro sorts may or may
       * not be active. Blocking them with dependencies deadlocks with MPI. So
       * add a new sort task instead, which will just do nothing if the cell is
       * already sorted. */
      t_rt_sorts = scheduler_addtask(s, task_type_rt_sort, task_subtype_none, 0,
                                     0, c, NULL);
      c->rt.rt_sorts = t_rt_sorts;
      if (c->hydro.sorts != NULL) {
        /* Copy task flags. While these should always be empty for sorts, better
         * be safe than spend hours looking for this. */
        t_rt_sorts->flags = c->hydro.sorts->flags;
        /* Make sure the normal hydro sorts run before the RT sorts run. */
        scheduler_addunlock(s, c->hydro.sorts, t_rt_sorts);
        /* Don't run gradients on unsorted cells. */
        scheduler_addunlock(s, c->hydro.sorts, t_rt_gradient);
      }

      /* Make sure the second receive doesn't get enqueued before the first one
       * is done */
      scheduler_addunlock(s, t_rt_gradient, t_rt_sorts);
      scheduler_addunlock(s, t_rt_gradient, t_rt_transport);
      /* Avoid situation where we receive while the sort hasn't finished yet. */
      scheduler_addunlock(s, t_rt_sorts, t_rt_transport);
      /* If one or both recv tasks are active, make sure the
       * rt_advance_cell_time tasks doesn't run before them */
      scheduler_addunlock(s, t_rt_gradient, c->super->rt.rt_advance_cell_time);
      scheduler_addunlock(s, t_rt_transport, c->super->rt.rt_advance_cell_time);
      /* Make sure the gradient recv don't run before the xv is finished.
       * This can occur when a cell itself is inactive for both hydro and
       * RT, but needs to be sent over for some other cell's pair task.
       * For active cells, you must make sure that t_rho and t_gradient have
       * been received first. As there is no guarantee which message will
       * arrive first, you might overwrite data otherwise. */

      scheduler_addunlock(s, t_xv, t_rt_gradient);
      scheduler_addunlock(s, t_rho, t_rt_gradient);
#ifdef EXTRA_HYDRO_LOOP
      scheduler_addunlock(s, t_gradient, t_rt_gradient);
#endif
    }
  }

  if (t_xv != NULL) {
    engine_addlink(e, &c->mpi.recv, t_xv);
    engine_addlink(e, &c->mpi.recv, t_rho);
#ifdef EXTRA_HYDRO_LOOP
    engine_addlink(e, &c->mpi.recv, t_gradient);
#endif
    if (with_limiter) {
      engine_addlink(e, &c->mpi.recv, t_limiter);
      engine_addlink(e, &c->mpi.unpack, t_unpack_limiter);
    }
#ifdef EXTRA_STAR_LOOPS
    if (with_feedback) engine_addlink(e, &c->mpi.recv, t_prep1);
#endif

    /* Add dependencies. */
    if (c->hydro.sorts != NULL) {
      scheduler_addunlock(s, t_xv, c->hydro.sorts);
      scheduler_addunlock(s, c->hydro.sorts, t_rho);
#if defined(MPI_SYMMETRIC_FORCE_INTERACTION) && defined(EXTRA_HYDRO_LOOP)
      scheduler_addunlock(s, c->hydro.sorts, t_gradient);
#endif
    }

    for (struct link *l = c->hydro.density; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_xv, l->t);
      scheduler_addunlock(s, l->t, t_rho);
    }
#ifdef EXTRA_HYDRO_LOOP
    for (struct link *l = c->hydro.gradient; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_rho, l->t);
      scheduler_addunlock(s, l->t, t_gradient);
    }
    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_gradient, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
#else
    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_rho, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
#endif

    if (with_limiter) {
      for (struct link *l = c->hydro.limiter; l != NULL; l = l->next) {
        scheduler_addunlock(s, t_unpack_limiter, l->t);
      }
    }

    /* Make sure the gas density has been computed before the
     * stars compute theirs. */
    if (with_feedback) {
      for (struct link *l = c->stars.density; l != NULL; l = l->next) {
        scheduler_addunlock(s, t_rho, l->t);
      }
    }
#ifdef EXTRA_STAR_LOOPS
    if (with_feedback) {
      /* Receive gas parts after everything is finished in prep1 loop */
      for (struct link *l = c->stars.prepare1; l != NULL; l = l->next) {
        scheduler_addunlock(s, l->t, t_prep1);
      }

      /* Start updating stars in prep2 only after the updated gas parts have
       * been received */
      for (struct link *l = c->stars.prepare2; l != NULL; l = l->next) {
        scheduler_addunlock(s, t_prep1, l->t);
      }
    }
#endif
    /* Make sure the part have been received before the BHs compute their
     * accretion rates (depends on particles' rho). */
    if (with_black_holes) {
      for (struct link *l = c->black_holes.density; l != NULL; l = l->next) {
        /* t_rho is not activated for cells with no active hydro, so we need
           to add an additional dependency on t_xv for these cells */
        scheduler_addunlock(s, t_xv, l->t);
        scheduler_addunlock(s, t_rho, l->t);
      }
    }

    if (with_rt) {
      engine_addlink(e, &c->mpi.recv, t_rt_gradient);
      engine_addlink(e, &c->mpi.recv, t_rt_transport);

      /* RT recvs mustn't run before hydro force has completed. */
      for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
        scheduler_addunlock(s, l->t, t_rt_gradient);
      }

      for (struct link *l = c->rt.rt_gradient; l != NULL; l = l->next) {
        /* RT gradient tasks mustn't run before we receive necessary data */
        scheduler_addunlock(s, t_rt_gradient, l->t);
        /* Don't run gradient tasks without sorting */
        scheduler_addunlock(s, t_rt_sorts, l->t);
        /* Don't update local particles before gradient tasks are finished */
        scheduler_addunlock(s, l->t, t_rt_transport);
      }

      for (struct link *l = c->rt.rt_transport; l != NULL; l = l->next) {
        /* RT transport tasks (iact, not comm tasks!!) mustn't run before we
         * receive necessary data */
        scheduler_addunlock(s, t_rt_transport, l->t);
        /* add dependency for the timestep communication tasks. In cases where
         * RT is inactive, rt_advance_cell_time won't run, so we need to make
         * sure we don't receive data before we're done with all the work. */
        scheduler_addunlock(s, l->t, tend);
        /* advance cell time mustn't run before transport is done */
        scheduler_addunlock(s, l->t, c->super->rt.rt_advance_cell_time);
      }
    }
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_hydro(
            e, c->progeny[k], t_xv, t_rho, t_gradient, t_prep1, t_limiter,
            t_unpack_limiter, t_rt_gradient, t_rt_transport, t_rt_sorts, tend,
            with_feedback, with_black_holes, with_limiter, with_sync, with_rt);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add time rt_advance_cell_time tasks to super levels of
 * foreign cells. This function recurses down to the super level
 * and creates the required tasks, and adds a dependency between
 * rt_advance_cell_time, rt_collect_times, and tend tasks.
 *
 * In normal steps, tend mustn't run before rt_advance_cell_time or the
 * cell's ti_rt_end_min will be updated wrongly. In sub-cycles, we don't
 * have the tend tasks, so there's no worry about that. (Them missing is
 * the reason we need the rt_advanced_cell_time to complete the
 * sub-cycles in the first place)
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param tend The top-level time-step communication #task.
 */

void engine_addtasks_recv_rt_advance_cell_time(struct engine *e, struct cell *c,
                                               struct task *const tend) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Have we reached the super level? */
  if (c->super == c) {

#ifdef SWIFT_RT_DEBUG_CHECKS
    if (c->super == NULL)
      error("trying to add rt_advance_cell_time above super level...");
    if (c->top->rt.rt_collect_times == NULL)
      error("rt_collect_times should have been created already????");
#endif

    /* Create the rt advance times task at the super level, if it hasn't
     * already. also set all the dependencies */
    if (c->rt.rt_advance_cell_time == NULL) {

      c->rt.rt_advance_cell_time = scheduler_addtask(
          s, task_type_rt_advance_cell_time, task_subtype_none, 0, 0, c, NULL);

      /* don't run collect times before you run advance cell time */
      scheduler_addunlock(s, c->rt.rt_advance_cell_time,
                          c->top->rt.rt_collect_times);

      /* Add the dependency */
      scheduler_addunlock(s, c->super->rt.rt_advance_cell_time, tend);
    }

    /* we're done. */
    return;
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_rt_advance_cell_time(e, c->progeny[k], tend);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for stars pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_density The recv_density #task, if it has already been created.
 * @param t_prep2 The recv_prep2 #task, if it has already been created.
 * @param t_sf_counts The recv_sf_counts, if it has been created.
 * @param tend The top-level time-step communication #task.
 * @param with_star_formation Are we running with star formation on?
 */
void engine_addtasks_recv_stars(struct engine *e, struct cell *c,
                                struct task *t_density, struct task *t_prep2,
                                struct task *t_sf_counts,
                                struct task *const tend,
                                const int with_star_formation) {
#ifdef SWIFT_DEBUG_CHECKS
  if (e->policy & engine_policy_sinks && e->policy & engine_policy_stars) {
    error("TODO");
  }
#endif

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  if (t_sf_counts == NULL && with_star_formation && c->hydro.count > 0) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->depth != 0)
      error(
          "Attaching a sf_count task at a non-top level c->depth=%d "
          "c->count=%d",
          c->depth, c->hydro.count);
#endif
    t_sf_counts = scheduler_addtask(s, task_type_recv, task_subtype_sf_counts,
                                    c->mpi.tag, 0, c, NULL);
  }

  /* Have we reached a level where there are any stars tasks ? */
  if (t_density == NULL && c->stars.density != NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    /* Create the tasks. */
    t_density = scheduler_addtask(s, task_type_recv, task_subtype_spart_density,
                                  c->mpi.tag, 0, c, NULL);

#ifdef EXTRA_STAR_LOOPS
    t_prep2 = scheduler_addtask(s, task_type_recv, task_subtype_spart_prep2,
                                c->mpi.tag, 0, c, NULL);
#endif
    if (with_star_formation && c->hydro.count > 0) {

      /* Receive the stars only once the counts have been received */
      scheduler_addunlock(s, t_sf_counts, c->stars.sorts);
      scheduler_addunlock(s, t_sf_counts, t_density);
#ifdef EXTRA_STAR_LOOPS
      scheduler_addunlock(s, t_sf_counts, t_prep2);
#endif
    }
  }

  if (t_density != NULL) {
    engine_addlink(e, &c->mpi.recv, t_density);
#ifdef EXTRA_STAR_LOOPS
    engine_addlink(e, &c->mpi.recv, t_prep2);
#endif
    if (with_star_formation && c->hydro.count > 0) {
      engine_addlink(e, &c->mpi.recv, t_sf_counts);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (c->nodeID == e->nodeID) error("Local cell!");
#endif
    if (c->stars.sorts != NULL) {
      scheduler_addunlock(s, t_density, c->stars.sorts);
#ifdef EXTRA_STAR_LOOPS
      scheduler_addunlock(s, c->stars.sorts, t_prep2);
#endif
    }

    /* Receive stars after the density loop */
    for (struct link *l = c->stars.density; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_density);
    }

#ifdef EXTRA_STAR_LOOPS
    /* Start updating local gas only after sparts have been received */
    for (struct link *l = c->stars.prepare1; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_density, l->t);
      scheduler_addunlock(s, l->t, t_prep2);
    }

    /* Receive stars for the second time after the prep2 loop */
    for (struct link *l = c->stars.prepare2; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_prep2);
    }

    /* Start updating local gas only after sparts have been received */
    for (struct link *l = c->stars.feedback; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_prep2, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
#else
    /* Start updating local gas only after sparts have been received */
    for (struct link *l = c->stars.feedback; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_density, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
#endif
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_stars(e, c->progeny[k], t_density, t_prep2,
                                   t_sf_counts, tend, with_star_formation);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for black_holes pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_rho The density comm. task, if it has already been created.
 * @param t_bh_merger The BH swallow comm. task, if it has already been created.
 * @param t_gas_swallow The gas swallow comm. task, if it has already been
 * created.
 * @param t_feedback The recv_feed #task, if it has already been created.
 * @param tend The top-level time-step communication #task.
 */
void engine_addtasks_recv_black_holes(struct engine *e, struct cell *c,
                                      struct task *t_rho,
                                      struct task *t_bh_merger,
                                      struct task *t_gas_swallow,
                                      struct task *t_feedback,
                                      struct task *const tend) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Have we reached a level where there are any black_holes tasks ? */
  if (t_rho == NULL && c->black_holes.density != NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    /* Create the tasks. */
    t_rho = scheduler_addtask(s, task_type_recv, task_subtype_bpart_rho,
                              c->mpi.tag, 0, c, NULL);

    t_bh_merger = scheduler_addtask(
        s, task_type_recv, task_subtype_bpart_merger, c->mpi.tag, 0, c, NULL);

    t_gas_swallow = scheduler_addtask(
        s, task_type_recv, task_subtype_part_swallow, c->mpi.tag, 0, c, NULL);

    t_feedback = scheduler_addtask(
        s, task_type_recv, task_subtype_bpart_feedback, c->mpi.tag, 0, c, NULL);
  }

  if (t_rho != NULL) {
    engine_addlink(e, &c->mpi.recv, t_rho);
    engine_addlink(e, &c->mpi.recv, t_bh_merger);
    engine_addlink(e, &c->mpi.recv, t_gas_swallow);
    engine_addlink(e, &c->mpi.recv, t_feedback);

#ifdef SWIFT_DEBUG_CHECKS
    if (c->nodeID == e->nodeID) error("Local cell!");
#endif

    for (struct link *l = c->black_holes.density; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_rho);
    }

    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_gas_swallow);
    }

    for (struct link *l = c->black_holes.swallow; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_rho, l->t);
      scheduler_addunlock(s, l->t, t_gas_swallow);
      scheduler_addunlock(s, l->t, t_bh_merger);
    }
    for (struct link *l = c->black_holes.do_gas_swallow; l != NULL;
         l = l->next) {
      scheduler_addunlock(s, t_gas_swallow, l->t);
    }
    for (struct link *l = c->black_holes.do_bh_swallow; l != NULL;
         l = l->next) {
      scheduler_addunlock(s, t_bh_merger, l->t);
      scheduler_addunlock(s, l->t, t_feedback);
    }
    for (struct link *l = c->black_holes.feedback; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_feedback, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_black_holes(e, c->progeny[k], t_rho, t_bh_merger,
                                         t_gas_swallow, t_feedback, tend);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_grav The recv_gpart #task, if it has already been created.
 * @param tend The top-level time-step communication #task.
 */
void engine_addtasks_recv_gravity(struct engine *e, struct cell *c,
                                  struct task *t_grav,
                                  struct task *const tend) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Have we reached a level where there are any gravity tasks ? */
  if (t_grav == NULL && c->grav.grav != NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    /* Create the tasks. */
    t_grav = scheduler_addtask(s, task_type_recv, task_subtype_gpart,
                               c->mpi.tag, 0, c, NULL);
  }

  /* If we have tasks, link them. */
  if (t_grav != NULL) {
    engine_addlink(e, &c->mpi.recv, t_grav);

    for (struct link *l = c->grav.grav; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_grav, l->t);
      scheduler_addunlock(s, l->t, tend);
    }
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_gravity(e, c->progeny[k], t_grav, tend);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- timestep version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_common(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;
  const int with_sinks = (e->policy & engine_policy_sinks);
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_star_formation = (e->policy & engine_policy_star_formation);
  const int with_star_formation_sink = with_sinks && with_stars;
  const int with_timestep_limiter =
      (e->policy & engine_policy_timestep_limiter);
  const int with_timestep_sync = (e->policy & engine_policy_timestep_sync);
  const int with_rt = (e->policy & engine_policy_rt);
#ifdef WITH_CSDS
  const int with_csds = e->policy & engine_policy_csds;
#endif

  /* Are we at the top-level? */
  if (c->top == c && c->nodeID == e->nodeID) {

    if (c->hydro.count > 0 || c->grav.count > 0 || c->stars.count > 0 ||
        c->black_holes.count > 0 || c->sinks.count > 0) {
      c->timestep_collect = scheduler_addtask(s, task_type_collect,
                                              task_subtype_none, 0, 0, c, NULL);
    }

    if (with_star_formation && c->hydro.count > 0) {
      c->hydro.star_formation = scheduler_addtask(
          s, task_type_star_formation, task_subtype_none, 0, 0, c, NULL);
    }

    if (with_star_formation_sink &&
        (c->hydro.count > 0 || c->sinks.count > 0)) {
      c->sinks.star_formation_sink = scheduler_addtask(
          s, task_type_star_formation_sink, task_subtype_none, 0, 0, c, NULL);
    }

    if (with_sinks) {
      /* sinks.sink_formation plays the role of a ghost => always created when
       * playing with sinks*/
      c->sinks.sink_formation = scheduler_addtask(
          s, task_type_sink_formation, task_subtype_none, 0, 0, c, NULL);
    }

    if (with_rt) {
      c->rt.rt_collect_times = scheduler_addtask(
          s, task_type_rt_collect_times, task_subtype_none, 0, 0, c, NULL);
    }
  }

  /* Are we in a super-cell ? */
  if (c->super == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the two half kicks */
      c->kick1 = scheduler_addtask(s, task_type_kick1, task_subtype_none, 0, 0,
                                   c, NULL);

      c->kick2 = scheduler_addtask(s, task_type_kick2, task_subtype_none, 0, 0,
                                   c, NULL);

      /* Weighting task for neutrinos after the last kick */
      if (e->neutrino_properties->use_delta_f) {
        c->grav.neutrino_weight = scheduler_addtask(
            s, task_type_neutrino_weight, task_subtype_none, 0, 0, c, NULL);
        scheduler_addunlock(s, c->kick1, c->grav.neutrino_weight);
      }

#if defined(WITH_CSDS)
      struct task *kick2_or_csds;
      if (with_csds) {
        /* Add the hydro csds task. */
        c->csds = scheduler_addtask(s, task_type_csds, task_subtype_none, 0, 0,
                                    c, NULL);

        /* Add the kick2 dependency */
        scheduler_addunlock(s, c->kick2, c->csds);

        /* Create a variable in order to avoid to many ifdef */
        kick2_or_csds = c->csds;
      } else {
        kick2_or_csds = c->kick2;
      }
#else
      struct task *kick2_or_csds = c->kick2;
#endif

      /* Add the time-step calculation task and its dependency */
      c->timestep = scheduler_addtask(s, task_type_timestep, task_subtype_none,
                                      0, 0, c, NULL);

      scheduler_addunlock(s, kick2_or_csds, c->timestep);
      scheduler_addunlock(s, c->timestep, c->kick1);
      scheduler_addunlock(s, c->timestep, c->top->timestep_collect);

      /* Subgrid tasks: star formation */
      if (with_star_formation && c->hydro.count > 0) {
        scheduler_addunlock(s, kick2_or_csds, c->top->hydro.star_formation);
        scheduler_addunlock(s, c->top->hydro.star_formation, c->timestep);
      }

      /* Subgrid tasks: star formation from sinks */
      if (with_star_formation_sink &&
          (c->hydro.count > 0 || c->sinks.count > 0)) {
        scheduler_addunlock(s, kick2_or_csds,
                            c->top->sinks.star_formation_sink);
        scheduler_addunlock(s, c->top->sinks.star_formation_sink, c->timestep);
      }

      /* Time-step limiter */
      if (with_timestep_limiter) {

        c->timestep_limiter = scheduler_addtask(
            s, task_type_timestep_limiter, task_subtype_none, 0, 0, c, NULL);

        scheduler_addunlock(s, c->timestep, c->timestep_limiter);
        scheduler_addunlock(s, c->timestep_limiter, c->kick1);
        scheduler_addunlock(s, c->timestep_limiter, c->top->timestep_collect);
      }

      /* Time-step synchronization */
      if (with_timestep_sync) {

        c->timestep_sync = scheduler_addtask(s, task_type_timestep_sync,
                                             task_subtype_none, 0, 0, c, NULL);

        scheduler_addunlock(s, c->timestep, c->timestep_sync);
        scheduler_addunlock(s, c->timestep_sync, c->kick1);
        scheduler_addunlock(s, c->timestep_sync, c->top->timestep_collect);
      }

      if (with_timestep_limiter && with_timestep_sync) {
        scheduler_addunlock(s, c->timestep_limiter, c->timestep_sync);
      }
    }
  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_common(e, c->progeny[k]);
  }
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- gravity version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_gravity(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;
  const int is_self_gravity = (e->policy & engine_policy_self_gravity);
  const int stars_only_gravity =
      (e->policy & engine_policy_stars) && !(e->policy & engine_policy_hydro);

  /* Are we in a super-cell ? */
  if (c->grav.super == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      if (stars_only_gravity) {

        /* In the special case where we have stars that just act under gravity
         * we must create their drift task here and not just copy over the hydro
         * behaviour. */
        c->stars.drift = scheduler_addtask(s, task_type_drift_spart,
                                           task_subtype_none, 0, 0, c, NULL);

        scheduler_addunlock(s, c->stars.drift, c->super->kick2);
      }

      c->grav.drift = scheduler_addtask(s, task_type_drift_gpart,
                                        task_subtype_none, 0, 0, c, NULL);

      c->grav.end_force = scheduler_addtask(s, task_type_end_grav_force,
                                            task_subtype_none, 0, 0, c, NULL);

      scheduler_addunlock(s, c->grav.end_force, c->super->kick2);

      if (is_self_gravity) {

        /* Initialisation of the multipoles */
        c->grav.init = scheduler_addtask(s, task_type_init_grav,
                                         task_subtype_none, 0, 0, c, NULL);

        /* Gravity non-neighbouring pm calculations */
        c->grav.long_range = scheduler_addtask(
            s, task_type_grav_long_range, task_subtype_none, 0, 0, c, NULL);

        /* Gravity recursive down-pass */
        c->grav.down = scheduler_addtask(s, task_type_grav_down,
                                         task_subtype_none, 0, 0, c, NULL);

        /* Implicit tasks for the up and down passes */
        c->grav.drift_out = scheduler_addtask(s, task_type_drift_gpart_out,
                                              task_subtype_none, 0, 1, c, NULL);
        c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                             task_subtype_none, 0, 1, c, NULL);
        c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                            task_subtype_none, 0, 1, c, NULL);

        /* Long-range gravity forces (not the mesh ones!) */
        scheduler_addunlock(s, c->grav.init, c->grav.long_range);
        scheduler_addunlock(s, c->grav.long_range, c->grav.down);
        scheduler_addunlock(s, c->grav.down, c->grav.super->grav.end_force);

        /* With adaptive softening, force the hydro density to complete first */
        if (gravity_after_hydro_density && c->hydro.super == c) {
          scheduler_addunlock(s, c->hydro.ghost_out, c->grav.init_out);
        }

        /* Link in the implicit tasks */
        scheduler_addunlock(s, c->grav.init, c->grav.init_out);
        scheduler_addunlock(s, c->grav.drift, c->grav.drift_out);
        scheduler_addunlock(s, c->grav.down_in, c->grav.down);
      }
    }
  }

  /* We are below the super-cell but not below the maximal splitting depth */
  else if ((c->grav.super != NULL) &&
           ((c->maxdepth - c->depth) >= space_subdepth_diff_grav)) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      if (is_self_gravity) {

        c->grav.drift_out = scheduler_addtask(s, task_type_drift_gpart_out,
                                              task_subtype_none, 0, 1, c, NULL);

        c->grav.init_out = scheduler_addtask(s, task_type_init_grav_out,
                                             task_subtype_none, 0, 1, c, NULL);

        c->grav.down_in = scheduler_addtask(s, task_type_grav_down_in,
                                            task_subtype_none, 0, 1, c, NULL);

        scheduler_addunlock(s, c->parent->grav.init_out, c->grav.init_out);
        scheduler_addunlock(s, c->parent->grav.drift_out, c->grav.drift_out);
        scheduler_addunlock(s, c->grav.down_in, c->parent->grav.down_in);
      }
    }
  }

  /* Recurse but not below the maximal splitting depth */
  if (c->split && ((c->maxdepth - c->depth) >= space_subdepth_diff_grav))
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_make_hierarchical_tasks_gravity(e, c->progeny[k]);
}

/**
 * @brief Recursively add non-implicit ghost tasks to a cell hierarchy.
 */
void engine_add_ghosts(struct engine *e, struct cell *c, struct task *ghost_in,
                       struct task *ghost_out) {

  /* Abort as there are no hydro particles here? */
  if (c->hydro.count_total == 0) return;

  /* If we have reached the leaf OR have to few particles to play with*/
  if (!c->split || c->hydro.count_total < engine_max_parts_per_ghost) {

    /* Add the ghost task and its dependencies */
    struct scheduler *s = &e->sched;
    c->hydro.ghost =
        scheduler_addtask(s, task_type_ghost, task_subtype_none, 0, 0, c, NULL);
    scheduler_addunlock(s, ghost_in, c->hydro.ghost);
    scheduler_addunlock(s, c->hydro.ghost, ghost_out);

  } else {
    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_add_ghosts(e, c->progeny[k], ghost_in, ghost_out);
  }
}

/**
 * @brief Recursively add non-implicit cooling tasks to a cell hierarchy.
 */
void engine_add_cooling(struct engine *e, struct cell *c,
                        struct task *cooling_in, struct task *cooling_out) {

  /* Abort as there are no hydro particles here? */
  if (c->hydro.count_total == 0) return;

  /* If we have reached the leaf OR have to few particles to play with*/
  if (!c->split || c->hydro.count_total < engine_max_parts_per_cooling) {

    /* Add the cooling task and its dependencies */
    struct scheduler *s = &e->sched;
    c->hydro.cooling = scheduler_addtask(s, task_type_cooling,
                                         task_subtype_none, 0, 0, c, NULL);
    scheduler_addunlock(s, cooling_in, c->hydro.cooling);
    scheduler_addunlock(s, c->hydro.cooling, cooling_out);

  } else {
    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_add_cooling(e, c->progeny[k], cooling_in, cooling_out);
  }
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- hydro version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param star_resort_cell Pointer to the cell where the star_resort task has
 * been created. NULL above that level or if not running with star formation.
 */
void engine_make_hierarchical_tasks_hydro(struct engine *e, struct cell *c,
                                          struct cell *star_resort_cell) {

  struct scheduler *s = &e->sched;
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_sinks = (e->policy & engine_policy_sinks);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_cooling = (e->policy & engine_policy_cooling);
  const int with_star_formation = (e->policy & engine_policy_star_formation);
  const int with_star_formation_sink = (with_sinks && with_stars);
  const int with_black_holes = (e->policy & engine_policy_black_holes);
  const int with_rt = (e->policy & engine_policy_rt);
#ifdef WITH_CSDS
  const int with_csds = (e->policy & engine_policy_csds);
#endif

  /* Are we are the level where we create the stars' resort tasks?
   * If the tree is shallow, we need to do this at the super-level if the
   * super-level is above the level we want */
  if ((c->nodeID == e->nodeID) && (star_resort_cell == NULL) &&
      (c->depth == engine_star_resort_task_depth || c->hydro.super == c)) {

    /* Star formation */
    if (with_feedback && c->hydro.count > 0 && with_star_formation) {

      /* Record this is the level where we re-sort */
      star_resort_cell = c;

      c->hydro.stars_resort = scheduler_addtask(
          s, task_type_stars_resort, task_subtype_none, 0, 0, c, NULL);

      scheduler_addunlock(s, c->top->hydro.star_formation,
                          c->hydro.stars_resort);
    }

    /* Star formation from sinks */
    if (with_feedback && with_star_formation_sink &&
        (c->hydro.count > 0 || c->sinks.count > 0)) {

      /* Record this is the level where we re-sort */
      star_resort_cell = c;

      c->hydro.stars_resort = scheduler_addtask(
          s, task_type_stars_resort, task_subtype_none, 0, 0, c, NULL);

      scheduler_addunlock(s, c->top->sinks.star_formation_sink,
                          c->hydro.stars_resort);
    }
  }

  /* Are we in a super-cell ? */
  if (c->hydro.super == c) {

    /* Add the sort task. */
    c->hydro.sorts =
        scheduler_addtask(s, task_type_sort, task_subtype_none, 0, 0, c, NULL);

    if (with_feedback) {
      c->stars.sorts = scheduler_addtask(s, task_type_stars_sort,
                                         task_subtype_none, 0, 0, c, NULL);
    }

    if (with_black_holes) {
      c->black_holes.swallow_ghost_1 =
          scheduler_addtask(s, task_type_bh_swallow_ghost1, task_subtype_none,
                            0, /* implicit =*/1, c, NULL);
    }

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the drift task. */
      c->hydro.drift = scheduler_addtask(s, task_type_drift_part,
                                         task_subtype_none, 0, 0, c, NULL);

      /* Add the task finishing the force calculation */
      c->hydro.end_force = scheduler_addtask(s, task_type_end_hydro_force,
                                             task_subtype_none, 0, 0, c, NULL);

      /* Generate the ghost tasks. */
      c->hydro.ghost_in =
          scheduler_addtask(s, task_type_ghost_in, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      c->hydro.ghost_out =
          scheduler_addtask(s, task_type_ghost_out, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      engine_add_ghosts(e, c, c->hydro.ghost_in, c->hydro.ghost_out);

      /* Generate the extra ghost task. */
#ifdef EXTRA_HYDRO_LOOP
      c->hydro.extra_ghost = scheduler_addtask(
          s, task_type_extra_ghost, task_subtype_none, 0, 0, c, NULL);
#endif

      /* Stars */
      if (with_stars) {
        c->stars.drift = scheduler_addtask(s, task_type_drift_spart,
                                           task_subtype_none, 0, 0, c, NULL);
        scheduler_addunlock(s, c->stars.drift, c->super->kick2);

        if (with_star_formation && c->top->hydro.count > 0)
          scheduler_addunlock(s, c->stars.drift, c->top->hydro.star_formation);
      }

      /* Sinks */
      if (with_sinks) {
        c->sinks.drift = scheduler_addtask(s, task_type_drift_sink,
                                           task_subtype_none, 0, 0, c, NULL);
        scheduler_addunlock(s, c->sinks.drift, c->super->kick2);

        c->sinks.sink_in =
            scheduler_addtask(s, task_type_sink_in, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->sinks.sink_ghost1 =
            scheduler_addtask(s, task_type_sink_ghost1, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->sinks.sink_ghost2 =
            scheduler_addtask(s, task_type_sink_ghost2, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->sinks.sink_out =
            scheduler_addtask(s, task_type_sink_out, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        /* Link to the main tasks */
        scheduler_addunlock(s, c->super->kick2, c->sinks.sink_in);
        scheduler_addunlock(s, c->sinks.sink_out, c->super->timestep);

        if (with_stars &&
            (c->top->hydro.count > 0 || c->top->sinks.count > 0)) {
          scheduler_addunlock(s, c->hydro.super->sinks.sink_out,
                              c->top->sinks.star_formation_sink);
        }
      }

      /* Black holes */
      if (with_black_holes) {
        c->black_holes.drift = scheduler_addtask(
            s, task_type_drift_bpart, task_subtype_none, 0, 0, c, NULL);
        scheduler_addunlock(s, c->black_holes.drift, c->super->kick2);
      }

      /* Subgrid tasks: cooling */
      if (with_cooling) {

        c->hydro.cooling_in =
            scheduler_addtask(s, task_type_cooling_in, task_subtype_none, 0,
                              /*implicit=*/1, c, NULL);
        c->hydro.cooling_out =
            scheduler_addtask(s, task_type_cooling_out, task_subtype_none, 0,
                              /*implicit=*/1, c, NULL);

        engine_add_cooling(e, c, c->hydro.cooling_in, c->hydro.cooling_out);

        scheduler_addunlock(s, c->hydro.end_force, c->hydro.cooling_in);
        scheduler_addunlock(s, c->hydro.cooling_out, c->super->kick2);

      } else {
        scheduler_addunlock(s, c->hydro.end_force, c->super->kick2);
      }

      /* Subgrid tasks: feedback */
      if (with_feedback) {

        c->stars.stars_in =
            scheduler_addtask(s, task_type_stars_in, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->stars.stars_out =
            scheduler_addtask(s, task_type_stars_out, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->stars.density_ghost = scheduler_addtask(
            s, task_type_stars_ghost, task_subtype_none, 0, 0, c, NULL);

#ifdef EXTRA_STAR_LOOPS
        c->stars.prep1_ghost =
            scheduler_addtask(s, task_type_stars_prep_ghost1, task_subtype_none,
                              0, /* implicit = */ 1, c, NULL);

        c->hydro.prep1_ghost =
            scheduler_addtask(s, task_type_hydro_prep_ghost1, task_subtype_none,
                              0, /* implicit = */ 1, c, NULL);

        c->stars.prep2_ghost =
            scheduler_addtask(s, task_type_stars_prep_ghost2, task_subtype_none,
                              0, /* implicit = */ 1, c, NULL);
#endif

#ifdef WITH_CSDS
        if (with_csds) {
          scheduler_addunlock(s, c->super->csds, c->stars.stars_in);
        } else {
          scheduler_addunlock(s, c->super->kick2, c->stars.stars_in);
        }
#else
        scheduler_addunlock(s, c->super->kick2, c->stars.stars_in);
#endif
        scheduler_addunlock(s, c->stars.stars_out, c->super->timestep);

        /* Star formation*/
        if (with_feedback && c->hydro.count > 0 && with_star_formation) {
          scheduler_addunlock(s, star_resort_cell->hydro.stars_resort,
                              c->stars.stars_in);
        }
        /* Star formation from sinks */
        if (with_feedback && with_star_formation_sink &&
            (c->hydro.count > 0 || c->sinks.count > 0)) {
          scheduler_addunlock(s, star_resort_cell->hydro.stars_resort,
                              c->stars.stars_in);
        }
      }

      /* Radiative Transfer */
      if (with_rt) {
        /* RT ghost in task */
        c->rt.rt_in =
            scheduler_addtask(s, task_type_rt_in, task_subtype_none, 0,
                              /* implicit= */ 1, c, NULL);
        scheduler_addunlock(s, c->super->kick2, c->rt.rt_in);
        /* Star formation */
        if (c->top->hydro.count > 0 && with_star_formation)
          scheduler_addunlock(s, c->top->hydro.star_formation, c->rt.rt_in);
        /* Star formation from sinks */
        if (with_star_formation_sink &&
            (c->top->hydro.count > 0 || c->top->sinks.count > 0))
          scheduler_addunlock(s, c->top->sinks.star_formation_sink,
                              c->rt.rt_in);
        if (with_feedback)
          scheduler_addunlock(s, c->stars.stars_out, c->rt.rt_in);
        /* TODO: check/add dependencies from Loic's new sink SF tasks */

        /* RT ghost out task */
        c->rt.rt_out =
            scheduler_addtask(s, task_type_rt_out, task_subtype_none, 0,
                              /* implicit= */ 1, c, NULL);
        scheduler_addunlock(s, c->rt.rt_out, c->super->timestep);

        /* In cases where nothing but RT is active, don't allow the timestep
         * collect to run before we've finished */
        scheduler_addunlock(s, c->rt.rt_out, c->top->timestep_collect);

        /* non-implicit ghost 1 */
        c->rt.rt_ghost1 = scheduler_addtask(s, task_type_rt_ghost1,
                                            task_subtype_none, 0, 0, c, NULL);
        scheduler_addunlock(s, c->rt.rt_in, c->rt.rt_ghost1);

        /* non-implicit ghost 2 */
        c->rt.rt_ghost2 = scheduler_addtask(s, task_type_rt_ghost2,
                                            task_subtype_none, 0, 0, c, NULL);

        /* implicit transport out */
        c->rt.rt_transport_out =
            scheduler_addtask(s, task_type_rt_transport_out, task_subtype_none,
                              0, /*implicit= */ 1, c, NULL);

        /* thermochemistry */
        c->rt.rt_tchem = scheduler_addtask(s, task_type_rt_tchem,
                                           task_subtype_none, 0, 0, c, NULL);

        /* Advance cell time for subcycling */
        /* We need to make sure that rt_advance_cell_time is at the same level
         * as the timestep task, not below. Otherwise, the updated cell times
         * won't propagate up the hierarchy enough, and the cell_is_rt_active
         * will return bogus results. Note that c->super is not necessarily
         * c->hydro.super in general. */
        /* Create the task only once ! */
        if (c->super->rt.rt_advance_cell_time == NULL) {
          c->super->rt.rt_advance_cell_time =
              scheduler_addtask(s, task_type_rt_advance_cell_time,
                                task_subtype_none, 0, 0, c->super, NULL);
          /* Don't run the rt_collect_times before the rt_advance_cell_time */
          scheduler_addunlock(s, c->super->rt.rt_advance_cell_time,
                              c->top->rt.rt_collect_times);
        }

        scheduler_addunlock(s, c->rt.rt_transport_out, c->rt.rt_tchem);
        scheduler_addunlock(s, c->rt.rt_tchem,
                            c->super->rt.rt_advance_cell_time);
        scheduler_addunlock(s, c->super->rt.rt_advance_cell_time, c->rt.rt_out);
      }

      /* Subgrid tasks: black hole feedback */
      if (with_black_holes) {

        c->black_holes.black_holes_in =
            scheduler_addtask(s, task_type_bh_in, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->black_holes.black_holes_out =
            scheduler_addtask(s, task_type_bh_out, task_subtype_none, 0,
                              /* implicit = */ 1, c, NULL);

        c->black_holes.density_ghost = scheduler_addtask(
            s, task_type_bh_density_ghost, task_subtype_none, 0, 0, c, NULL);

        c->black_holes.swallow_ghost_2 =
            scheduler_addtask(s, task_type_bh_swallow_ghost2, task_subtype_none,
                              0, /* implicit =*/1, c, NULL);

        c->black_holes.swallow_ghost_3 = scheduler_addtask(
            s, task_type_bh_swallow_ghost3, task_subtype_none, 0, 0, c, NULL);

#ifdef WITH_CSDS
        if (with_csds) {
          scheduler_addunlock(s, c->super->csds, c->black_holes.black_holes_in);
        } else {
          scheduler_addunlock(s, c->super->kick2,
                              c->black_holes.black_holes_in);
        }
#else
        scheduler_addunlock(s, c->super->kick2, c->black_holes.black_holes_in);
#endif
        scheduler_addunlock(s, c->black_holes.black_holes_out,
                            c->super->timestep);
      }

      if (with_black_holes && with_feedback) {

        /* Make sure we don't start swallowing gas particles before the stars
           have converged on their smoothing lengths. */
        scheduler_addunlock(s, c->stars.density_ghost,
                            c->black_holes.swallow_ghost_1);
      }
    }
  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_hydro(e, c->progeny[k],
                                               star_resort_cell);
  }
}

void engine_make_hierarchical_tasks_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_self_gravity = (e->policy & engine_policy_self_gravity);
  const int with_ext_gravity = (e->policy & engine_policy_external_gravity);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &((struct cell *)map_data)[ind];
    /* Make the common tasks (time integration) */
    engine_make_hierarchical_tasks_common(e, c);
    /* Add the hydro stuff */
    if (with_hydro)
      engine_make_hierarchical_tasks_hydro(e, c, /*star_resort_cell=*/NULL);
    /* And the gravity stuff */
    if (with_self_gravity || with_ext_gravity)
      engine_make_hierarchical_tasks_gravity(e, c);
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 */
void engine_make_self_gravity_tasks_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

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
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

#ifdef WITH_MPI
          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");
#endif

          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (periodic && min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm(ci, cj, e, s, /*use_rebuild_data=*/1,
                                    /*is_tree_walk=*/0)) {

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
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
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == cj) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cjd);
            } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[ci->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", ci->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == ci) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cid);
            }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the external gravity.
 *
 * @param e The #engine.
 */
void engine_make_external_gravity_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;

  for (int cid = 0; cid < nr_cells; ++cid) {

    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* Is that neighbour local ? */
    if (ci->nodeID != nodeID) continue;

    /* If the cell is local, build a self-interaction */
    scheduler_addtask(sched, task_type_self, task_subtype_external_grav, 0, 0,
                      ci, NULL);
  }
}

/**
 * @brief Counts the tasks associated with one cell and constructs the links
 *
 * For each hydrodynamic and gravity task, construct the links with
 * the corresponding cell.  Similarly, construct the dependencies for
 * all the sorting tasks.
 */
void engine_count_and_link_tasks_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *const sched = &e->sched;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &((struct task *)map_data)[ind];

    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

    /* Link sort tasks to all the higher sort task. */
    if (t_type == task_type_sort) {
      for (struct cell *finger = t->ci->parent; finger != NULL;
           finger = finger->parent) {
        if (finger->hydro.sorts != NULL)
          scheduler_addunlock(sched, t, finger->hydro.sorts);
      }

      /* Link stars sort tasks to all the higher sort task. */
    } else if (t_type == task_type_stars_sort) {
      for (struct cell *finger = t->ci->parent; finger != NULL;
           finger = finger->parent) {
        if (finger->stars.sorts != NULL)
          scheduler_addunlock(sched, t, finger->stars.sorts);
      }

      /* Link self tasks to cells. */
    } else if (t_type == task_type_self) {
      atomic_inc(&ci->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->hydro.density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      } else if (t_subtype == task_subtype_external_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      }

      /* Link pair tasks to cells. */
    } else if (t_type == task_type_pair) {
      atomic_inc(&ci->nr_tasks);
      atomic_inc(&cj->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->hydro.density, t);
        engine_addlink(e, &cj->hydro.density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav.grav, t);
        engine_addlink(e, &cj->grav.grav, t);
      }
#ifdef SWIFT_DEBUG_CHECKS
      else if (t_subtype == task_subtype_external_grav) {
        error("Found a pair/external-gravity task...");
      }
#endif

      /* Link sub-self tasks to cells. */
    } else if (t_type == task_type_sub_self) {
      atomic_inc(&ci->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->hydro.density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      } else if (t_subtype == task_subtype_external_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      }

      /* Link sub-pair tasks to cells. */
    } else if (t_type == task_type_sub_pair) {
      atomic_inc(&ci->nr_tasks);
      atomic_inc(&cj->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->hydro.density, t);
        engine_addlink(e, &cj->hydro.density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav.grav, t);
        engine_addlink(e, &cj->grav.grav, t);
      }
#ifdef SWIFT_DEBUG_CHECKS
      else if (t_subtype == task_subtype_external_grav) {
        error("Found a sub-pair/external-gravity task...");
      }
#endif

      /* Multipole-multipole interaction of progenies */
    } else if (t_type == task_type_grav_mm) {

      atomic_inc(&ci->grav.nr_mm_tasks);
      atomic_inc(&cj->grav.nr_mm_tasks);
      engine_addlink(e, &ci->grav.mm, t);
      engine_addlink(e, &cj->grav.mm, t);
    }
  }
}

/**
 * @brief Creates all the task dependencies for the gravity
 *
 * @param e The #engine
 */
void engine_link_gravity_tasks(struct engine *e) {

  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int nr_tasks = sched->nr_tasks;

  for (int k = 0; k < nr_tasks; k++) {

    /* Get a pointer to the task. */
    struct task *t = &sched->tasks[k];

    if (t->type == task_type_none) continue;

    /* Get the cells we act on */
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

    /* Pointers to the parent cells for tasks going up and down the tree
     * In the case where we are at the super-level we don't
     * want the parent as no tasks are defined above that level. */
    struct cell *ci_parent, *cj_parent;
    if (ci->parent != NULL && ci->grav.super != ci)
      ci_parent = ci->parent;
    else
      ci_parent = ci;

    if (cj != NULL && cj->parent != NULL && cj->grav.super != cj)
      cj_parent = cj->parent;
    else
      cj_parent = cj;

/* Node ID (if running with MPI) */
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Self-interaction for self-gravity? */
    if (t_type == task_type_self && t_subtype == task_subtype_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local self task");
#endif

      /* drift ---+-> gravity --> grav_down */
      /* init  --/    */
      scheduler_addunlock(sched, ci_parent->grav.drift_out, t);
      scheduler_addunlock(sched, ci_parent->grav.init_out, t);
      scheduler_addunlock(sched, t, ci_parent->grav.down_in);
    }

    /* Self-interaction for external gravity ? */
    if (t_type == task_type_self && t_subtype == task_subtype_external_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local self task");
#endif

      /* drift -----> gravity --> end_gravity_force */
      scheduler_addunlock(sched, ci->grav.super->grav.drift, t);
      scheduler_addunlock(sched, t, ci->grav.super->grav.end_force);
    }

    /* Otherwise, pair interaction? */
    else if (t_type == task_type_pair && t_subtype == task_subtype_grav) {

      if (ci_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        scheduler_addunlock(sched, ci_parent->grav.drift_out, t);
        scheduler_addunlock(sched, ci_parent->grav.init_out, t);
        scheduler_addunlock(sched, t, ci_parent->grav.down_in);
      }
      if (cj_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        if (ci_parent != cj_parent) { /* Avoid double unlock */
          scheduler_addunlock(sched, cj_parent->grav.drift_out, t);
          scheduler_addunlock(sched, cj_parent->grav.init_out, t);
          scheduler_addunlock(sched, t, cj_parent->grav.down_in);
        }
      }
    }

    /* Otherwise, sub-self interaction? */
    else if (t_type == task_type_sub_self && t_subtype == task_subtype_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local sub-self task");
#endif
      /* drift ---+-> gravity --> grav_down */
      /* init  --/    */
      scheduler_addunlock(sched, ci_parent->grav.drift_out, t);
      scheduler_addunlock(sched, ci_parent->grav.init_out, t);
      scheduler_addunlock(sched, t, ci_parent->grav.down_in);
    }

    /* Sub-self-interaction for external gravity ? */
    else if (t_type == task_type_sub_self &&
             t_subtype == task_subtype_external_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local sub-self task");
#endif

      /* drift -----> gravity --> end_force */
      scheduler_addunlock(sched, ci->grav.super->grav.drift, t);
      scheduler_addunlock(sched, t, ci->grav.super->grav.end_force);
    }

    /* Otherwise, sub-pair interaction? */
    else if (t_type == task_type_sub_pair && t_subtype == task_subtype_grav) {

      if (ci_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        scheduler_addunlock(sched, ci_parent->grav.drift_out, t);
        scheduler_addunlock(sched, ci_parent->grav.init_out, t);
        scheduler_addunlock(sched, t, ci_parent->grav.down_in);
      }
      if (cj_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        if (ci_parent != cj_parent) { /* Avoid double unlock */
          scheduler_addunlock(sched, cj_parent->grav.drift_out, t);
          scheduler_addunlock(sched, cj_parent->grav.init_out, t);
          scheduler_addunlock(sched, t, cj_parent->grav.down_in);
        }
      }

    }

    /* Otherwise M-M interaction? */
    else if (t_type == task_type_grav_mm) {

      if (ci_nodeID == nodeID) {

        /* init -----> gravity --> grav_down */
        scheduler_addunlock(sched, ci_parent->grav.init_out, t);
        scheduler_addunlock(sched, t, ci_parent->grav.down_in);
      }
      if (cj_nodeID == nodeID) {

        /* init -----> gravity --> grav_down */
        if (ci_parent != cj_parent) { /* Avoid double unlock */
          scheduler_addunlock(sched, cj_parent->grav.init_out, t);
          scheduler_addunlock(sched, t, cj_parent->grav.down_in);
        }
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP

/**
 * @brief Creates the dependency network for the hydro tasks of a given cell.
 *
 * @param sched The #scheduler.
 * @param density The density task to link.
 * @param gradient The gradient task to link.
 * @param force The force task to link.
 * @param limiter The limiter task to link.
 * @param c The cell.
 * @param with_cooling Do we have a cooling task ?
 * @param with_limiter Do we have a time-step limiter ?
 */
static inline void engine_make_hydro_loops_dependencies(
    struct scheduler *sched, struct task *density, struct task *gradient,
    struct task *force, struct task *limiter, struct cell *c, int with_cooling,
    int with_limiter) {

  /* density loop --> ghost --> gradient loop --> extra_ghost */
  /* extra_ghost --> force loop  */
  scheduler_addunlock(sched, density, c->hydro.super->hydro.ghost_in);
  scheduler_addunlock(sched, c->hydro.super->hydro.ghost_out, gradient);
  scheduler_addunlock(sched, gradient, c->hydro.super->hydro.extra_ghost);
  scheduler_addunlock(sched, c->hydro.super->hydro.extra_ghost, force);
}

#else

/**
 * @brief Creates the dependency network for the hydro tasks of a given cell.
 *
 * @param sched The #scheduler.
 * @param density The density task to link.
 * @param force The force task to link.
 * @param limiter The limiter task to link.
 * @param c The cell.
 * @param with_cooling Are we running with cooling switched on?
 * @param with_limiter Are we running with limiter switched on?
 */
static inline void engine_make_hydro_loops_dependencies(
    struct scheduler *sched, struct task *density, struct task *force,
    struct task *limiter, struct cell *c, int with_cooling, int with_limiter) {

  /* density loop --> ghost --> force loop */
  scheduler_addunlock(sched, density, c->hydro.super->hydro.ghost_in);
  scheduler_addunlock(sched, c->hydro.super->hydro.ghost_out, force);
}

#endif

/**
 * @brief Duplicates the first hydro loop and construct all the
 * dependencies for the hydro part
 *
 * This is done by looping over all the previously constructed tasks
 * and adding another task involving the same cells but this time
 * corresponding to the second hydro loop over neighbours.
 * With all the relevant tasks for a given cell available, we construct
 * all the dependencies for that cell.
 */
void engine_make_extra_hydroloop_tasks_mapper(void *map_data, int num_elements,
                                              void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int with_cooling = (e->policy & engine_policy_cooling);
  const int with_timestep_limiter =
      (e->policy & engine_policy_timestep_limiter);
  const int with_timestep_sync = (e->policy & engine_policy_timestep_sync);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_black_holes = (e->policy & engine_policy_black_holes);
  const int with_rt = (e->policy & engine_policy_rt);
  const int with_sink = (e->policy & engine_policy_sinks);
#ifdef EXTRA_HYDRO_LOOP
  struct task *t_gradient = NULL;
#endif
#ifdef EXTRA_STAR_LOOPS
  struct task *t_star_prep1 = NULL;
  struct task *t_star_prep2 = NULL;
#endif
  struct task *t_force = NULL;
  struct task *t_limiter = NULL;
  struct task *t_star_density = NULL;
  struct task *t_star_feedback = NULL;
  struct task *t_bh_density = NULL;
  struct task *t_bh_swallow = NULL;
  struct task *t_do_gas_swallow = NULL;
  struct task *t_do_bh_swallow = NULL;
  struct task *t_bh_feedback = NULL;
  struct task *t_sink_swallow = NULL;
  struct task *t_rt_gradient = NULL;
  struct task *t_rt_transport = NULL;
  struct task *t_sink_do_sink_swallow = NULL;
  struct task *t_sink_do_gas_swallow = NULL;

  for (int ind = 0; ind < num_elements; ind++) {

    struct task *t = &((struct task *)map_data)[ind];
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;
    const long long flags = t->flags;
    struct cell *const ci = t->ci;
    struct cell *const cj = t->cj;

    /* Escape early */
    if (t->type == task_type_none) continue;
    if (t->type == task_type_stars_resort) continue;
    if (t->type == task_type_star_formation) continue;
    if (t->type == task_type_star_formation_sink) continue;
    if (t->type == task_type_sink_formation) continue;

    /* Sort tasks depend on the drift of the cell (gas version). */
    if (t_type == task_type_sort && ci->nodeID == nodeID) {
      scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t);
    }

    /* Sort tasks depend on the drift of the cell (stars version). */
    else if (t_type == task_type_stars_sort && ci->nodeID == nodeID) {
      scheduler_addunlock(sched, ci->hydro.super->stars.drift, t);
    }

    /* Self-interaction? */
    else if (t_type == task_type_self && t_subtype == task_subtype_density) {

      const int bcount_i = ci->black_holes.count;

      /* Make the self-density tasks depend on the drift only. */
      scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t);

      /* Task for the second hydro loop, */
      t_force = scheduler_addtask(sched, task_type_self, task_subtype_force,
                                  flags, 0, ci, NULL);

      /* the task for the time-step limiter */
      if (with_timestep_limiter) {
        t_limiter = scheduler_addtask(sched, task_type_self,
                                      task_subtype_limiter, flags, 0, ci, NULL);
      }

      /* The stellar feedback tasks */
      if (with_feedback) {
        t_star_density =
            scheduler_addtask(sched, task_type_self, task_subtype_stars_density,
                              flags, 0, ci, NULL);
        t_star_feedback =
            scheduler_addtask(sched, task_type_self,
                              task_subtype_stars_feedback, flags, 0, ci, NULL);

#ifdef EXTRA_STAR_LOOPS
        t_star_prep1 =
            scheduler_addtask(sched, task_type_self, task_subtype_stars_prep1,
                              flags, 0, ci, NULL);
        t_star_prep2 =
            scheduler_addtask(sched, task_type_self, task_subtype_stars_prep2,
                              flags, 0, ci, NULL);
#endif
      }

      /* The sink tasks */
      if (with_sink) {
        t_sink_swallow =
            scheduler_addtask(sched, task_type_self, task_subtype_sink_swallow,
                              flags, 0, ci, NULL);
        t_sink_do_sink_swallow = scheduler_addtask(
            sched, task_type_self, task_subtype_sink_do_sink_swallow, flags, 0,
            ci, NULL);
        t_sink_do_gas_swallow = scheduler_addtask(
            sched, task_type_self, task_subtype_sink_do_gas_swallow, flags, 0,
            ci, NULL);
      }

      /* The black hole feedback tasks */
      if (with_black_holes && bcount_i > 0) {
        t_bh_density = scheduler_addtask(
            sched, task_type_self, task_subtype_bh_density, flags, 0, ci, NULL);
        t_bh_swallow = scheduler_addtask(
            sched, task_type_self, task_subtype_bh_swallow, flags, 0, ci, NULL);
        t_do_gas_swallow =
            scheduler_addtask(sched, task_type_self,
                              task_subtype_do_gas_swallow, flags, 0, ci, NULL);
        t_do_bh_swallow =
            scheduler_addtask(sched, task_type_self, task_subtype_do_bh_swallow,
                              flags, 0, ci, NULL);
        t_bh_feedback =
            scheduler_addtask(sched, task_type_self, task_subtype_bh_feedback,
                              flags, 0, ci, NULL);
      }

      if (with_rt) {
        t_rt_gradient =
            scheduler_addtask(sched, task_type_self, task_subtype_rt_gradient,
                              flags, 0, ci, NULL);
        t_rt_transport =
            scheduler_addtask(sched, task_type_self, task_subtype_rt_transport,
                              flags, 0, ci, NULL);
      }

      /* Link the tasks to the cells */
      engine_addlink(e, &ci->hydro.force, t_force);
      if (with_timestep_limiter) {
        engine_addlink(e, &ci->hydro.limiter, t_limiter);
      }
      if (with_feedback) {
        engine_addlink(e, &ci->stars.density, t_star_density);
        engine_addlink(e, &ci->stars.feedback, t_star_feedback);
#ifdef EXTRA_STAR_LOOPS
        engine_addlink(e, &ci->stars.prepare1, t_star_prep1);
        engine_addlink(e, &ci->stars.prepare2, t_star_prep2);
#endif
      }
      if (with_sink) {
        engine_addlink(e, &ci->sinks.swallow, t_sink_swallow);
        engine_addlink(e, &ci->sinks.do_sink_swallow, t_sink_do_sink_swallow);
        engine_addlink(e, &ci->sinks.do_gas_swallow, t_sink_do_gas_swallow);
      }
      if (with_black_holes && bcount_i > 0) {
        engine_addlink(e, &ci->black_holes.density, t_bh_density);
        engine_addlink(e, &ci->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &ci->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &ci->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &ci->black_holes.feedback, t_bh_feedback);
      }
      if (with_rt) {
        engine_addlink(e, &ci->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &ci->rt.rt_transport, t_rt_transport);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Same work for the additional hydro loop */
      t_gradient = scheduler_addtask(sched, task_type_self,
                                     task_subtype_gradient, flags, 0, ci, NULL);

      /* Add the link between the new loops and the cell */
      engine_addlink(e, &ci->hydro.gradient, t_gradient);

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                           t_limiter, ci, with_cooling,
                                           with_timestep_limiter);
#else

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, ci,
                                           with_cooling, with_timestep_limiter);
#endif

      /* Create the task dependencies */
      scheduler_addunlock(sched, t_force, ci->hydro.super->hydro.end_force);

      if (with_feedback) {

        if (with_cooling)
          scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                              t_star_density);

        scheduler_addunlock(sched, ci->hydro.super->stars.drift,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->stars.stars_in,
                            t_star_density);
        scheduler_addunlock(sched, t_star_density,
                            ci->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
        scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                            t_star_prep1);
        scheduler_addunlock(sched, t_star_prep1,
                            ci->hydro.super->stars.prep1_ghost);
        scheduler_addunlock(sched, t_star_prep1,
                            ci->hydro.super->hydro.prep1_ghost);
        scheduler_addunlock(sched, ci->hydro.super->stars.prep1_ghost,
                            t_star_prep2);
        scheduler_addunlock(sched, ci->hydro.super->hydro.prep1_ghost,
                            t_star_prep2);
        scheduler_addunlock(sched, t_star_prep2,
                            ci->hydro.super->stars.prep2_ghost);
        scheduler_addunlock(sched, ci->hydro.super->stars.prep2_ghost,
                            t_star_feedback);
#else
        scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                            t_star_feedback);
#endif
        scheduler_addunlock(sched, t_star_feedback,
                            ci->hydro.super->stars.stars_out);
      }

      /* The sink's tasks. */
      if (with_sink) {

        /* Do the sink_formation */
        scheduler_addunlock(sched, ci->hydro.super->sinks.drift,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_in,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->top->sinks.sink_formation,
                            t_sink_swallow);

        /* Do the sink_swallow */
        scheduler_addunlock(sched, t_sink_swallow,
                            ci->hydro.super->sinks.sink_ghost1);

        /* Do the sink_do_gas_swallow */
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost1,
                            t_sink_do_gas_swallow);
        scheduler_addunlock(sched, t_sink_do_gas_swallow,
                            ci->hydro.super->sinks.sink_ghost2);

        /* Do the sink_do_sink_swallow */
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost2,
                            t_sink_do_sink_swallow);
        scheduler_addunlock(sched, t_sink_do_sink_swallow,
                            ci->hydro.super->sinks.sink_out);
      }

      if (with_black_holes && bcount_i > 0) {

        if (with_cooling)
          scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                              t_bh_density);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.drift,
                            t_bh_density);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_bh_density);
        scheduler_addunlock(sched, ci->hydro.super->black_holes.black_holes_in,
                            t_bh_density);
        scheduler_addunlock(sched, t_bh_density,
                            ci->hydro.super->black_holes.density_ghost);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.density_ghost,
                            t_bh_swallow);
        scheduler_addunlock(sched, t_bh_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_1);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_1,
                            t_do_gas_swallow);
        scheduler_addunlock(sched, t_do_gas_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_2);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_2,
                            t_do_bh_swallow);
        scheduler_addunlock(sched, t_do_bh_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_3);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_3,
                            t_bh_feedback);
        scheduler_addunlock(sched, t_bh_feedback,
                            ci->hydro.super->black_holes.black_holes_out);
      }

      if (with_timestep_limiter) {
        scheduler_addunlock(sched, ci->super->timestep, t_limiter);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_limiter);
        scheduler_addunlock(sched, t_limiter, ci->super->kick1);
        scheduler_addunlock(sched, t_limiter, ci->super->timestep_limiter);
      }

      if (with_timestep_sync && with_feedback) {
        scheduler_addunlock(sched, t_star_feedback, ci->super->timestep_sync);
      }
      if (with_timestep_sync && with_black_holes && bcount_i > 0) {
        scheduler_addunlock(sched, t_bh_feedback, ci->super->timestep_sync);
      }

      if (with_rt) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_rt_gradient);
        scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost1,
                            t_rt_gradient);
        scheduler_addunlock(sched, t_rt_gradient,
                            ci->hydro.super->rt.rt_ghost2);
        scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost2,
                            t_rt_transport);
        scheduler_addunlock(sched, t_rt_transport,
                            ci->hydro.super->rt.rt_transport_out);
      }
    }

    /* Otherwise, pair interaction? */
    else if (t_type == task_type_pair && t_subtype == task_subtype_density) {

      const int bcount_i = ci->black_holes.count;
      const int bcount_j = cj->black_holes.count;

      /* Make all density tasks depend on the drift */
      if (ci->nodeID == nodeID) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.drift, t);
      }

      /* Make all density tasks depend on the sorts */
      scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t);
      if (ci->hydro.super != cj->hydro.super) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.sorts, t);
      }

      /* New task for the force */
      t_force = scheduler_addtask(sched, task_type_pair, task_subtype_force,
                                  flags, 0, ci, cj);

#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
      /* The order of operations for an inactive local cell interacting
       * with an active foreign cell is not guaranteed because the density
       * (and gradient) iact loops don't exist in that case. So we need
       * an explicit dependency here to have sorted cells. */

      /* Make all force tasks depend on the sorts */
      scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t_force);
      if (ci->hydro.super != cj->hydro.super) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.sorts, t_force);
      }
#endif

      /* and the task for the time-step limiter */
      if (with_timestep_limiter) {
        t_limiter = scheduler_addtask(sched, task_type_pair,
                                      task_subtype_limiter, flags, 0, ci, cj);
      }

      /* The stellar feedback tasks */
      if (with_feedback) {
        t_star_density =
            scheduler_addtask(sched, task_type_pair, task_subtype_stars_density,
                              flags, 0, ci, cj);
        t_star_feedback =
            scheduler_addtask(sched, task_type_pair,
                              task_subtype_stars_feedback, flags, 0, ci, cj);
#ifdef EXTRA_STAR_LOOPS
        t_star_prep1 = scheduler_addtask(
            sched, task_type_pair, task_subtype_stars_prep1, flags, 0, ci, cj);
        t_star_prep2 = scheduler_addtask(
            sched, task_type_pair, task_subtype_stars_prep2, flags, 0, ci, cj);
#endif
      }

      /* The sink tasks */
      if (with_sink) {
        t_sink_swallow = scheduler_addtask(
            sched, task_type_pair, task_subtype_sink_swallow, flags, 0, ci, cj);
        t_sink_do_sink_swallow = scheduler_addtask(
            sched, task_type_pair, task_subtype_sink_do_sink_swallow, flags, 0,
            ci, cj);
        t_sink_do_gas_swallow = scheduler_addtask(
            sched, task_type_pair, task_subtype_sink_do_gas_swallow, flags, 0,
            ci, cj);
      }

      /* The black hole feedback tasks */
      if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
        t_bh_density = scheduler_addtask(
            sched, task_type_pair, task_subtype_bh_density, flags, 0, ci, cj);
        t_bh_swallow = scheduler_addtask(
            sched, task_type_pair, task_subtype_bh_swallow, flags, 0, ci, cj);
        t_do_gas_swallow =
            scheduler_addtask(sched, task_type_pair,
                              task_subtype_do_gas_swallow, flags, 0, ci, cj);
        t_do_bh_swallow =
            scheduler_addtask(sched, task_type_pair, task_subtype_do_bh_swallow,
                              flags, 0, ci, cj);
        t_bh_feedback = scheduler_addtask(
            sched, task_type_pair, task_subtype_bh_feedback, flags, 0, ci, cj);
      }

      if (with_rt) {
        t_rt_gradient = scheduler_addtask(
            sched, task_type_pair, task_subtype_rt_gradient, flags, 0, ci, cj);
        t_rt_transport = scheduler_addtask(
            sched, task_type_pair, task_subtype_rt_transport, flags, 0, ci, cj);
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
        /* The order of operations for an inactive local cell interacting
         * with an active foreign cell is not guaranteed because the gradient
         * iact loops don't exist in that case. So we need an explicit
         * dependency here to have sorted cells. */

        /* Make all force tasks depend on the sorts */
        if (ci->hydro.super->rt.rt_sorts != NULL)
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_sorts,
                              t_rt_transport);
        if (ci->hydro.super != cj->hydro.super) {
          if (cj->hydro.super->rt.rt_sorts != NULL)
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_sorts,
                                t_rt_transport);
        }
        /* We need to ensure that a local inactive cell is sorted before
         * the interaction in the transport loop. Local cells don't have an
         * rt_sorts task. */
        if (ci->hydro.super->hydro.sorts != NULL)
          scheduler_addunlock(sched, ci->hydro.super->hydro.sorts,
                              t_rt_transport);
        if ((ci->hydro.super != cj->hydro.super) &&
            (cj->hydro.super->hydro.sorts != NULL))
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_rt_transport);
#endif
      }

      engine_addlink(e, &ci->hydro.force, t_force);
      engine_addlink(e, &cj->hydro.force, t_force);
      if (with_timestep_limiter) {
        engine_addlink(e, &ci->hydro.limiter, t_limiter);
        engine_addlink(e, &cj->hydro.limiter, t_limiter);
      }
      if (with_feedback) {
        engine_addlink(e, &ci->stars.density, t_star_density);
        engine_addlink(e, &cj->stars.density, t_star_density);
        engine_addlink(e, &ci->stars.feedback, t_star_feedback);
        engine_addlink(e, &cj->stars.feedback, t_star_feedback);
#ifdef EXTRA_STAR_LOOPS
        engine_addlink(e, &ci->stars.prepare1, t_star_prep1);
        engine_addlink(e, &cj->stars.prepare1, t_star_prep1);
        engine_addlink(e, &ci->stars.prepare2, t_star_prep2);
        engine_addlink(e, &cj->stars.prepare2, t_star_prep2);
#endif
      }
      if (with_sink) {
        /* Formation */
        engine_addlink(e, &ci->sinks.swallow, t_sink_swallow);
        engine_addlink(e, &cj->sinks.swallow, t_sink_swallow);
        /* Merger */
        engine_addlink(e, &ci->sinks.do_sink_swallow, t_sink_do_sink_swallow);
        engine_addlink(e, &cj->sinks.do_sink_swallow, t_sink_do_sink_swallow);
        /* Accretion */
        engine_addlink(e, &ci->sinks.do_gas_swallow, t_sink_do_gas_swallow);
        engine_addlink(e, &cj->sinks.do_gas_swallow, t_sink_do_gas_swallow);
      }
      if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
        engine_addlink(e, &ci->black_holes.density, t_bh_density);
        engine_addlink(e, &cj->black_holes.density, t_bh_density);
        engine_addlink(e, &ci->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &cj->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &ci->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &cj->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &ci->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &cj->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &ci->black_holes.feedback, t_bh_feedback);
        engine_addlink(e, &cj->black_holes.feedback, t_bh_feedback);
      }
      if (with_rt) {
        engine_addlink(e, &ci->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &cj->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &ci->rt.rt_transport, t_rt_transport);
        engine_addlink(e, &cj->rt.rt_transport, t_rt_transport);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      t_gradient = scheduler_addtask(sched, task_type_pair,
                                     task_subtype_gradient, flags, 0, ci, cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &ci->hydro.gradient, t_gradient);
      engine_addlink(e, &cj->hydro.gradient, t_gradient);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, ci, with_cooling,
                                             with_timestep_limiter);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, cj, with_cooling,
                                             with_timestep_limiter);
      }
#else

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, ci,
                                             with_cooling,
                                             with_timestep_limiter);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, cj,
                                             with_cooling,
                                             with_timestep_limiter);
      }
#endif

      if (with_feedback) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts,
                            t_star_density);

        if (ci->hydro.super != cj->hydro.super) {
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_star_density);
        }
      }
      if (with_rt) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t_rt_gradient);

        if (ci->hydro.super != cj->hydro.super) {
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_rt_gradient);
        }
      }

      if (ci->nodeID == nodeID) {
        scheduler_addunlock(sched, t_force, ci->hydro.super->hydro.end_force);

        if (with_feedback) {

          if (with_cooling)
            scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                                t_star_density);

          scheduler_addunlock(sched, ci->hydro.super->stars.drift,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->stars.stars_in,
                              t_star_density);
          scheduler_addunlock(sched, t_star_density,
                              ci->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                              t_star_prep1);
          scheduler_addunlock(sched, t_star_prep1,
                              ci->hydro.super->stars.prep1_ghost);
          scheduler_addunlock(sched, t_star_prep1,
                              ci->hydro.super->hydro.prep1_ghost);
          scheduler_addunlock(sched, ci->hydro.super->stars.prep1_ghost,
                              t_star_prep2);
          scheduler_addunlock(sched, ci->hydro.super->hydro.prep1_ghost,
                              t_star_prep2);
          scheduler_addunlock(sched, t_star_prep2,
                              ci->hydro.super->stars.prep2_ghost);
          scheduler_addunlock(sched, ci->hydro.super->stars.prep2_ghost,
                              t_star_feedback);
#else
          scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                              t_star_feedback);
#endif
          scheduler_addunlock(sched, t_star_feedback,
                              ci->hydro.super->stars.stars_out);
        }

        if (with_sink) {

          /* Do the sink_formation */
          scheduler_addunlock(sched, ci->hydro.super->sinks.drift,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_in,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->top->sinks.sink_formation,
                              t_sink_swallow);

          /* Do the sink_swallow */
          scheduler_addunlock(sched, t_sink_swallow,
                              ci->hydro.super->sinks.sink_ghost1);

          /* Do the sink_do_gas_swallow */
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost1,
                              t_sink_do_gas_swallow);
          scheduler_addunlock(sched, t_sink_do_gas_swallow,
                              ci->hydro.super->sinks.sink_ghost2);

          /* Do the sink_do_sink_swallow */
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost2,
                              t_sink_do_sink_swallow);
          scheduler_addunlock(sched, t_sink_do_sink_swallow,
                              ci->hydro.super->sinks.sink_out);
        }

        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

          if (with_cooling)
            scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                                t_bh_density);

          scheduler_addunlock(sched, ci->hydro.super->black_holes.drift,
                              t_bh_density);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_bh_density);
          scheduler_addunlock(
              sched, ci->hydro.super->black_holes.black_holes_in, t_bh_density);
          scheduler_addunlock(sched, t_bh_density,
                              ci->hydro.super->black_holes.density_ghost);

          scheduler_addunlock(sched, ci->hydro.super->black_holes.density_ghost,
                              t_bh_swallow);
          scheduler_addunlock(sched, t_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_1);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_1,
                              t_do_gas_swallow);
          scheduler_addunlock(sched, t_do_gas_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_2);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_2,
                              t_do_bh_swallow);
          scheduler_addunlock(sched, t_do_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_3);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_3,
                              t_bh_feedback);
          scheduler_addunlock(sched, t_bh_feedback,
                              ci->hydro.super->black_holes.black_holes_out);
        }

        if (with_timestep_limiter) {
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_limiter);
          scheduler_addunlock(sched, ci->super->timestep, t_limiter);
          scheduler_addunlock(sched, t_limiter, ci->super->kick1);
          scheduler_addunlock(sched, t_limiter, ci->super->timestep_limiter);
        }

        if (with_timestep_sync && with_feedback) {
          scheduler_addunlock(sched, t_star_feedback, ci->super->timestep_sync);
        }
        if (with_timestep_sync && with_black_holes &&
            (bcount_i > 0 || bcount_j > 0)) {
          scheduler_addunlock(sched, t_bh_feedback, ci->super->timestep_sync);
        }

        if (with_rt) {
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_rt_gradient);
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost1,
                              t_rt_gradient);
          scheduler_addunlock(sched, t_rt_gradient,
                              ci->hydro.super->rt.rt_ghost2);
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost2,
                              t_rt_transport);
          scheduler_addunlock(sched, t_rt_transport,
                              ci->hydro.super->rt.rt_transport_out);
        }

      } else /*(ci->nodeID != nodeID) */ {
        if (with_feedback) {
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_prep1);
#endif
          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_feedback);
        }

        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
          scheduler_addunlock(sched, t_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_1);
        }
      }

      if (cj->nodeID == nodeID) {

        if (ci->hydro.super != cj->hydro.super) {

          scheduler_addunlock(sched, t_force, cj->hydro.super->hydro.end_force);

          if (with_feedback) {

            if (with_cooling)
              scheduler_addunlock(sched, cj->hydro.super->hydro.cooling_out,
                                  t_star_density);

            scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->stars.drift,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->stars.stars_in,
                                t_star_density);
            scheduler_addunlock(sched, t_star_density,
                                cj->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
            scheduler_addunlock(sched, cj->hydro.super->stars.density_ghost,
                                t_star_prep1);
            scheduler_addunlock(sched, t_star_prep1,
                                cj->hydro.super->stars.prep1_ghost);
            scheduler_addunlock(sched, t_star_prep1,
                                cj->hydro.super->hydro.prep1_ghost);
            scheduler_addunlock(sched, cj->hydro.super->stars.prep1_ghost,
                                t_star_prep2);
            scheduler_addunlock(sched, cj->hydro.super->hydro.prep1_ghost,
                                t_star_prep2);
            scheduler_addunlock(sched, t_star_prep2,
                                cj->hydro.super->stars.prep2_ghost);
            scheduler_addunlock(sched, cj->hydro.super->stars.prep2_ghost,
                                t_star_feedback);
#else
            scheduler_addunlock(sched, cj->hydro.super->stars.density_ghost,
                                t_star_feedback);
#endif
            scheduler_addunlock(sched, t_star_feedback,
                                cj->hydro.super->stars.stars_out);
          }

          if (with_sink) {

            /* Do the sink_formation */
            scheduler_addunlock(sched, cj->hydro.super->sinks.drift,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_in,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->top->sinks.sink_formation,
                                t_sink_swallow);

            /* Do the sink_swallow */
            scheduler_addunlock(sched, t_sink_swallow,
                                cj->hydro.super->sinks.sink_ghost1);

            /* Do the sink_do_gas_swallow */
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_ghost1,
                                t_sink_do_gas_swallow);
            scheduler_addunlock(sched, t_sink_do_gas_swallow,
                                cj->hydro.super->sinks.sink_ghost2);

            /* Do the sink_do_sink_swallow */
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_ghost2,
                                t_sink_do_sink_swallow);
            scheduler_addunlock(sched, t_sink_do_sink_swallow,
                                cj->hydro.super->sinks.sink_out);
          }

          if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

            if (with_cooling)
              scheduler_addunlock(sched, cj->hydro.super->hydro.cooling_out,
                                  t_bh_density);

            scheduler_addunlock(sched, cj->hydro.super->black_holes.drift,
                                t_bh_density);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_bh_density);
            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.black_holes_in,
                                t_bh_density);
            scheduler_addunlock(sched, t_bh_density,
                                cj->hydro.super->black_holes.density_ghost);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.density_ghost,
                                t_bh_swallow);
            scheduler_addunlock(sched, t_bh_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_1);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_1,
                                t_do_gas_swallow);
            scheduler_addunlock(sched, t_do_gas_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_2);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_2,
                                t_do_bh_swallow);
            scheduler_addunlock(sched, t_do_bh_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_3);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_3,
                                t_bh_feedback);
            scheduler_addunlock(sched, t_bh_feedback,
                                cj->hydro.super->black_holes.black_holes_out);
          }

          if (with_rt) {
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_rt_gradient);
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_ghost1,
                                t_rt_gradient);
            scheduler_addunlock(sched, t_rt_gradient,
                                cj->hydro.super->rt.rt_ghost2);
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_ghost2,
                                t_rt_transport);
            scheduler_addunlock(sched, t_rt_transport,
                                cj->hydro.super->rt.rt_transport_out);
          }

          if (with_timestep_limiter) {
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift, t_limiter);
          }
        }

        if (ci->super != cj->super) {

          if (with_timestep_limiter) {
            scheduler_addunlock(sched, cj->super->timestep, t_limiter);
            scheduler_addunlock(sched, t_limiter, cj->super->kick1);
            scheduler_addunlock(sched, t_limiter, cj->super->timestep_limiter);
          }

          if (with_timestep_sync && with_feedback) {
            scheduler_addunlock(sched, t_star_feedback,
                                cj->super->timestep_sync);
          }
          if (with_timestep_sync && with_black_holes &&
              (bcount_i > 0 || bcount_j > 0)) {
            scheduler_addunlock(sched, t_bh_feedback, cj->super->timestep_sync);
          }
        }

      } else /*(cj->nodeID != nodeID) */ {
        if (with_feedback) {
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                              t_star_prep1);
#endif
          scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                              t_star_feedback);
        }

        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

          scheduler_addunlock(sched, t_bh_swallow,
                              cj->hydro.super->black_holes.swallow_ghost_1);
        }
      }
    }

    /* Otherwise, sub-self interaction? */
    else if (t_type == task_type_sub_self &&
             t_subtype == task_subtype_density) {

      const int bcount_i = ci->black_holes.count;

      /* Make all density tasks depend on the drift and sorts. */
      scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t);

      /* Start by constructing the task for the second hydro loop */
      t_force = scheduler_addtask(sched, task_type_sub_self, task_subtype_force,
                                  flags, 0, ci, NULL);

      /* and the task for the time-step limiter */
      if (with_timestep_limiter) {
        t_limiter = scheduler_addtask(sched, task_type_sub_self,
                                      task_subtype_limiter, flags, 0, ci, NULL);
      }

      /* The stellar feedback tasks */
      if (with_feedback) {
        t_star_density =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_stars_density, flags, 0, ci, NULL);
        t_star_feedback =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_stars_feedback, flags, 0, ci, NULL);

#ifdef EXTRA_STAR_LOOPS
        t_star_prep1 =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_stars_prep1, flags, 0, ci, NULL);
        t_star_prep2 =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_stars_prep2, flags, 0, ci, NULL);
#endif
      }

      /* The sink tasks */
      if (with_sink) {
        t_sink_swallow =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_sink_swallow, flags, 0, ci, NULL);
        t_sink_do_sink_swallow = scheduler_addtask(
            sched, task_type_sub_self, task_subtype_sink_do_sink_swallow, flags,
            0, ci, NULL);
        t_sink_do_gas_swallow = scheduler_addtask(
            sched, task_type_sub_self, task_subtype_sink_do_gas_swallow, flags,
            0, ci, NULL);
      }

      /* The black hole feedback tasks */
      if (with_black_holes && bcount_i > 0) {
        t_bh_density =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_bh_density, flags, 0, ci, NULL);
        t_bh_swallow =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_bh_swallow, flags, 0, ci, NULL);

        t_do_gas_swallow =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_do_gas_swallow, flags, 0, ci, NULL);

        t_do_bh_swallow =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_do_bh_swallow, flags, 0, ci, NULL);

        t_bh_feedback =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_bh_feedback, flags, 0, ci, NULL);
      }

      if (with_rt) {
        t_rt_gradient =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_rt_gradient, flags, 0, ci, NULL);
        t_rt_transport =
            scheduler_addtask(sched, task_type_sub_self,
                              task_subtype_rt_transport, flags, 0, ci, NULL);
      }

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &ci->hydro.force, t_force);
      if (with_timestep_limiter) {
        engine_addlink(e, &ci->hydro.limiter, t_limiter);
      }
      if (with_feedback) {
        engine_addlink(e, &ci->stars.density, t_star_density);
        engine_addlink(e, &ci->stars.feedback, t_star_feedback);
#ifdef EXTRA_STAR_LOOPS
        engine_addlink(e, &ci->stars.prepare1, t_star_prep1);
        engine_addlink(e, &ci->stars.prepare2, t_star_prep2);
#endif
      }
      if (with_sink) {
        engine_addlink(e, &ci->sinks.swallow, t_sink_swallow);
        engine_addlink(e, &ci->sinks.do_sink_swallow, t_sink_do_sink_swallow);
        engine_addlink(e, &ci->sinks.do_gas_swallow, t_sink_do_gas_swallow);
      }
      if (with_black_holes && bcount_i > 0) {
        engine_addlink(e, &ci->black_holes.density, t_bh_density);
        engine_addlink(e, &ci->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &ci->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &ci->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &ci->black_holes.feedback, t_bh_feedback);
      }
      if (with_rt) {
        engine_addlink(e, &ci->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &ci->rt.rt_transport, t_rt_transport);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      t_gradient = scheduler_addtask(sched, task_type_sub_self,
                                     task_subtype_gradient, flags, 0, ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &ci->hydro.gradient, t_gradient);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                           t_limiter, ci, with_cooling,
                                           with_timestep_limiter);
#else

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, ci,
                                           with_cooling, with_timestep_limiter);
#endif

      /* Create the task dependencies */
      scheduler_addunlock(sched, t_force, ci->hydro.super->hydro.end_force);

      if (with_feedback) {

        if (with_cooling)
          scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                              t_star_density);

        scheduler_addunlock(sched, ci->hydro.super->stars.drift,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts,
                            t_star_density);
        scheduler_addunlock(sched, ci->hydro.super->stars.stars_in,
                            t_star_density);
        scheduler_addunlock(sched, t_star_density,
                            ci->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
        scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                            t_star_prep1);
        scheduler_addunlock(sched, t_star_prep1,
                            ci->hydro.super->stars.prep1_ghost);
        scheduler_addunlock(sched, t_star_prep1,
                            ci->hydro.super->hydro.prep1_ghost);
        scheduler_addunlock(sched, ci->hydro.super->stars.prep1_ghost,
                            t_star_prep2);
        scheduler_addunlock(sched, ci->hydro.super->hydro.prep1_ghost,
                            t_star_prep2);
        scheduler_addunlock(sched, t_star_prep2,
                            ci->hydro.super->stars.prep2_ghost);
        scheduler_addunlock(sched, ci->hydro.super->stars.prep2_ghost,
                            t_star_feedback);
#else
        scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                            t_star_feedback);
#endif
        scheduler_addunlock(sched, t_star_feedback,
                            ci->hydro.super->stars.stars_out);
      }

      if (with_sink) {

        /* Do the sink_formation */
        scheduler_addunlock(sched, ci->hydro.super->sinks.drift,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_in,
                            ci->top->sinks.sink_formation);
        scheduler_addunlock(sched, ci->top->sinks.sink_formation,
                            t_sink_swallow);

        /* Do the sink_swallow */
        scheduler_addunlock(sched, t_sink_swallow,
                            ci->hydro.super->sinks.sink_ghost1);

        /* Do the sink_do_gas_swallow */
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost1,
                            t_sink_do_gas_swallow);
        scheduler_addunlock(sched, t_sink_do_gas_swallow,
                            ci->hydro.super->sinks.sink_ghost2);

        /* Do the sink_do_sink_swallow */
        scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost2,
                            t_sink_do_sink_swallow);
        scheduler_addunlock(sched, t_sink_do_sink_swallow,
                            ci->hydro.super->sinks.sink_out);
      }

      if (with_black_holes && bcount_i > 0) {

        if (with_cooling)
          scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                              t_bh_density);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.drift,
                            t_bh_density);
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_bh_density);
        scheduler_addunlock(sched, ci->hydro.super->black_holes.black_holes_in,
                            t_bh_density);
        scheduler_addunlock(sched, t_bh_density,
                            ci->hydro.super->black_holes.density_ghost);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.density_ghost,
                            t_bh_swallow);
        scheduler_addunlock(sched, t_bh_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_1);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_1,
                            t_do_gas_swallow);
        scheduler_addunlock(sched, t_do_gas_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_2);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_2,
                            t_do_bh_swallow);
        scheduler_addunlock(sched, t_do_bh_swallow,
                            ci->hydro.super->black_holes.swallow_ghost_3);

        scheduler_addunlock(sched, ci->hydro.super->black_holes.swallow_ghost_3,
                            t_bh_feedback);
        scheduler_addunlock(sched, t_bh_feedback,
                            ci->hydro.super->black_holes.black_holes_out);
      }

      if (with_timestep_limiter) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_limiter);
        scheduler_addunlock(sched, ci->super->timestep, t_limiter);
        scheduler_addunlock(sched, t_limiter, ci->super->kick1);
        scheduler_addunlock(sched, t_limiter, ci->super->timestep_limiter);
      }

      if (with_timestep_sync && with_feedback) {
        scheduler_addunlock(sched, t_star_feedback, ci->super->timestep_sync);
      }
      if (with_timestep_sync && with_black_holes && bcount_i > 0) {
        scheduler_addunlock(sched, t_bh_feedback, ci->super->timestep_sync);
      }

      if (with_rt) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_rt_gradient);
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t_rt_gradient);
        scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost1,
                            t_rt_gradient);
        scheduler_addunlock(sched, t_rt_gradient,
                            ci->hydro.super->rt.rt_ghost2);
        scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost2,
                            t_rt_transport);
        scheduler_addunlock(sched, t_rt_transport,
                            ci->hydro.super->rt.rt_transport_out);
      }
    }

    /* Otherwise, sub-pair interaction? */
    else if (t_type == task_type_sub_pair &&
             t_subtype == task_subtype_density) {

      const int bcount_i = ci->black_holes.count;
      const int bcount_j = cj->black_holes.count;

      /* Make all density tasks depend on the drift */
      if (ci->nodeID == nodeID) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.drift, t);
      }

      /* Make all density tasks depend on the sorts */
      scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t);
      if (ci->hydro.super != cj->hydro.super) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.sorts, t);
      }

      /* New task for the force */
      t_force = scheduler_addtask(sched, task_type_sub_pair, task_subtype_force,
                                  flags, 0, ci, cj);

#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
      /* The order of operations for an inactive local cell interacting
       * with an active foreign cell is not guaranteed because the density
       * (and gradient) iact loops don't exist in that case. So we need
       * an explicit dependency here to have sorted cells. */

      /* Make all force tasks depend on the sorts */
      scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t_force);
      if (ci->hydro.super != cj->hydro.super) {
        scheduler_addunlock(sched, cj->hydro.super->hydro.sorts, t_force);
      }
#endif

      /* and the task for the time-step limiter */
      if (with_timestep_limiter) {
        t_limiter = scheduler_addtask(sched, task_type_sub_pair,
                                      task_subtype_limiter, flags, 0, ci, cj);
      }

      /* The stellar feedback tasks */
      if (with_feedback) {
        t_star_density =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_stars_density, flags, 0, ci, cj);
        t_star_feedback =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_stars_feedback, flags, 0, ci, cj);

#ifdef EXTRA_STAR_LOOPS
        t_star_prep1 =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_stars_prep1, flags, 0, ci, cj);
        t_star_prep2 =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_stars_prep2, flags, 0, ci, cj);
#endif
      }

      /* The sink tasks */
      if (with_sink) {
        t_sink_swallow =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_sink_swallow, flags, 0, ci, cj);
        t_sink_do_sink_swallow = scheduler_addtask(
            sched, task_type_sub_pair, task_subtype_sink_do_sink_swallow, flags,
            0, ci, cj);
        t_sink_do_gas_swallow = scheduler_addtask(
            sched, task_type_sub_pair, task_subtype_sink_do_gas_swallow, flags,
            0, ci, cj);
      }

      /* The black hole feedback tasks */
      if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
        t_bh_density =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_bh_density, flags, 0, ci, cj);
        t_bh_swallow =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_bh_swallow, flags, 0, ci, cj);
        t_do_gas_swallow =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_do_gas_swallow, flags, 0, ci, cj);
        t_do_bh_swallow =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_do_bh_swallow, flags, 0, ci, cj);
        t_bh_feedback =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_bh_feedback, flags, 0, ci, cj);
      }

      if (with_rt) {
        t_rt_gradient =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_rt_gradient, flags, 0, ci, cj);
        t_rt_transport =
            scheduler_addtask(sched, task_type_sub_pair,
                              task_subtype_rt_transport, flags, 0, ci, cj);
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
        /* The order of operations for an inactive local cell interacting
         * with an active foreign cell is not guaranteed because the gradient
         * iact loops don't exist in that case. So we need an explicit
         * dependency here to have sorted cells. */

        /* Make all force tasks depend on the sorts */
        if (ci->hydro.super->rt.rt_sorts != NULL)
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_sorts,
                              t_rt_transport);
        if (ci->hydro.super != cj->hydro.super) {
          if (cj->hydro.super->rt.rt_sorts != NULL)
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_sorts,
                                t_rt_transport);
        }
        /* We need to ensure that a local inactive cell is sorted before
         * the interaction in the transport loop. Local cells don't have
         * an rt_sort task. */
        if (ci->hydro.super->hydro.sorts != NULL)
          scheduler_addunlock(sched, ci->hydro.super->hydro.sorts,
                              t_rt_transport);
        if ((ci->hydro.super != cj->hydro.super) &&
            (cj->hydro.super->hydro.sorts != NULL))
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_rt_transport);
#endif
      }

      engine_addlink(e, &ci->hydro.force, t_force);
      engine_addlink(e, &cj->hydro.force, t_force);
      if (with_timestep_limiter) {
        engine_addlink(e, &ci->hydro.limiter, t_limiter);
        engine_addlink(e, &cj->hydro.limiter, t_limiter);
      }
      if (with_feedback) {
        engine_addlink(e, &ci->stars.density, t_star_density);
        engine_addlink(e, &cj->stars.density, t_star_density);
        engine_addlink(e, &ci->stars.feedback, t_star_feedback);
        engine_addlink(e, &cj->stars.feedback, t_star_feedback);
#ifdef EXTRA_STAR_LOOPS
        engine_addlink(e, &ci->stars.prepare1, t_star_prep1);
        engine_addlink(e, &cj->stars.prepare1, t_star_prep1);
        engine_addlink(e, &ci->stars.prepare2, t_star_prep2);
        engine_addlink(e, &cj->stars.prepare2, t_star_prep2);
#endif
      }
      if (with_sink) {
        /* Formation */
        engine_addlink(e, &ci->sinks.swallow, t_sink_swallow);
        engine_addlink(e, &cj->sinks.swallow, t_sink_swallow);

        /* Merger */
        engine_addlink(e, &ci->sinks.do_sink_swallow, t_sink_do_sink_swallow);
        engine_addlink(e, &cj->sinks.do_sink_swallow, t_sink_do_sink_swallow);

        /* Accretion */
        engine_addlink(e, &ci->sinks.do_gas_swallow, t_sink_do_gas_swallow);
        engine_addlink(e, &cj->sinks.do_gas_swallow, t_sink_do_gas_swallow);
      }
      if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
        engine_addlink(e, &ci->black_holes.density, t_bh_density);
        engine_addlink(e, &cj->black_holes.density, t_bh_density);
        engine_addlink(e, &ci->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &cj->black_holes.swallow, t_bh_swallow);
        engine_addlink(e, &ci->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &cj->black_holes.do_gas_swallow, t_do_gas_swallow);
        engine_addlink(e, &ci->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &cj->black_holes.do_bh_swallow, t_do_bh_swallow);
        engine_addlink(e, &ci->black_holes.feedback, t_bh_feedback);
        engine_addlink(e, &cj->black_holes.feedback, t_bh_feedback);
      }
      if (with_rt) {
        engine_addlink(e, &ci->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &cj->rt.rt_gradient, t_rt_gradient);
        engine_addlink(e, &ci->rt.rt_transport, t_rt_transport);
        engine_addlink(e, &cj->rt.rt_transport, t_rt_transport);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      t_gradient = scheduler_addtask(sched, task_type_sub_pair,
                                     task_subtype_gradient, flags, 0, ci, cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &ci->hydro.gradient, t_gradient);
      engine_addlink(e, &cj->hydro.gradient, t_gradient);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, ci, with_cooling,
                                             with_timestep_limiter);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, cj, with_cooling,
                                             with_timestep_limiter);
      }
#else

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, ci,
                                             with_cooling,
                                             with_timestep_limiter);
      }
      if ((cj->nodeID == nodeID) && (ci->hydro.super != cj->hydro.super)) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, cj,
                                             with_cooling,
                                             with_timestep_limiter);
      }
#endif

      if (with_feedback) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts,
                            t_star_density);
        if (ci->hydro.super != cj->hydro.super) {
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_star_density);
        }
      }

      if (with_rt) {
        scheduler_addunlock(sched, ci->hydro.super->hydro.sorts, t_rt_gradient);
        if (ci->hydro.super != cj->hydro.super) {
          scheduler_addunlock(sched, cj->hydro.super->hydro.sorts,
                              t_rt_gradient);
        }
      }

      if (ci->nodeID == nodeID) {
        scheduler_addunlock(sched, t_force, ci->hydro.super->hydro.end_force);

        if (with_feedback) {

          if (with_cooling)
            scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                                t_star_density);

          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->stars.drift,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_star_density);
          scheduler_addunlock(sched, ci->hydro.super->stars.stars_in,
                              t_star_density);
          scheduler_addunlock(sched, t_star_density,
                              ci->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                              t_star_prep1);
          scheduler_addunlock(sched, t_star_prep1,
                              ci->hydro.super->stars.prep1_ghost);
          scheduler_addunlock(sched, t_star_prep1,
                              ci->hydro.super->hydro.prep1_ghost);
          scheduler_addunlock(sched, ci->hydro.super->stars.prep1_ghost,
                              t_star_prep2);
          scheduler_addunlock(sched, ci->hydro.super->hydro.prep1_ghost,
                              t_star_prep2);
          scheduler_addunlock(sched, t_star_prep2,
                              ci->hydro.super->stars.prep2_ghost);
          scheduler_addunlock(sched, ci->hydro.super->stars.prep2_ghost,
                              t_star_feedback);
#else
          scheduler_addunlock(sched, ci->hydro.super->stars.density_ghost,
                              t_star_feedback);
#endif
          scheduler_addunlock(sched, t_star_feedback,
                              ci->hydro.super->stars.stars_out);
        }

        if (with_sink) {

          /* Do the sink_formation */
          scheduler_addunlock(sched, ci->hydro.super->sinks.drift,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_in,
                              ci->top->sinks.sink_formation);
          scheduler_addunlock(sched, ci->top->sinks.sink_formation,
                              t_sink_swallow);

          /* Do the sink_swallow */
          scheduler_addunlock(sched, t_sink_swallow,
                              ci->hydro.super->sinks.sink_ghost1);

          /* Do the sink_do_gas_swallow */
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost1,
                              t_sink_do_gas_swallow);
          scheduler_addunlock(sched, t_sink_do_gas_swallow,
                              ci->hydro.super->sinks.sink_ghost2);

          /* Do the sink_do_sink_swallow */
          scheduler_addunlock(sched, ci->hydro.super->sinks.sink_ghost2,
                              t_sink_do_sink_swallow);
          scheduler_addunlock(sched, t_sink_do_sink_swallow,
                              ci->hydro.super->sinks.sink_out);
        }

        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

          if (with_cooling)
            scheduler_addunlock(sched, ci->hydro.super->hydro.cooling_out,
                                t_bh_density);

          scheduler_addunlock(sched, ci->hydro.super->black_holes.drift,
                              t_bh_density);
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_bh_density);
          scheduler_addunlock(
              sched, ci->hydro.super->black_holes.black_holes_in, t_bh_density);
          scheduler_addunlock(sched, t_bh_density,
                              ci->hydro.super->black_holes.density_ghost);

          scheduler_addunlock(sched, ci->hydro.super->black_holes.density_ghost,
                              t_bh_swallow);
          scheduler_addunlock(sched, t_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_1);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_1,
                              t_do_gas_swallow);
          scheduler_addunlock(sched, t_do_gas_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_2);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_2,
                              t_do_bh_swallow);
          scheduler_addunlock(sched, t_do_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_3);

          scheduler_addunlock(sched,
                              ci->hydro.super->black_holes.swallow_ghost_3,
                              t_bh_feedback);
          scheduler_addunlock(sched, t_bh_feedback,
                              ci->hydro.super->black_holes.black_holes_out);
        }

        if (with_timestep_limiter) {
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift, t_limiter);
          scheduler_addunlock(sched, ci->super->timestep, t_limiter);
          scheduler_addunlock(sched, t_limiter, ci->super->kick1);
          scheduler_addunlock(sched, t_limiter, ci->super->timestep_limiter);
        }

        if (with_timestep_sync && with_feedback) {
          scheduler_addunlock(sched, t_star_feedback, ci->super->timestep_sync);
        }
        if (with_timestep_sync && with_black_holes &&
            (bcount_i > 0 || bcount_j > 0)) {
          scheduler_addunlock(sched, t_bh_feedback, ci->super->timestep_sync);
        }

        if (with_rt) {
          scheduler_addunlock(sched, ci->hydro.super->hydro.drift,
                              t_rt_gradient);
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost1,
                              t_rt_gradient);
          scheduler_addunlock(sched, t_rt_gradient,
                              ci->hydro.super->rt.rt_ghost2);
          scheduler_addunlock(sched, ci->hydro.super->rt.rt_ghost2,
                              t_rt_transport);
          scheduler_addunlock(sched, t_rt_transport,
                              ci->hydro.super->rt.rt_transport_out);
        }
      } else /* ci->nodeID != nodeID */ {

        if (with_feedback) {
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_prep1);
#endif
          scheduler_addunlock(sched, ci->hydro.super->stars.sorts,
                              t_star_feedback);
        }
        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

          scheduler_addunlock(sched, t_bh_swallow,
                              ci->hydro.super->black_holes.swallow_ghost_1);
        }
      }

      if (cj->nodeID == nodeID) {

        if (ci->hydro.super != cj->hydro.super) {

          scheduler_addunlock(sched, t_force, cj->hydro.super->hydro.end_force);

          if (with_feedback) {

            if (with_cooling)
              scheduler_addunlock(sched, cj->hydro.super->hydro.cooling_out,
                                  t_star_density);

            scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->stars.drift,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_star_density);
            scheduler_addunlock(sched, cj->hydro.super->stars.stars_in,
                                t_star_density);
            scheduler_addunlock(sched, t_star_density,
                                cj->hydro.super->stars.density_ghost);
#ifdef EXTRA_STAR_LOOPS
            scheduler_addunlock(sched, cj->hydro.super->stars.density_ghost,
                                t_star_prep1);
            scheduler_addunlock(sched, t_star_prep1,
                                cj->hydro.super->stars.prep1_ghost);
            scheduler_addunlock(sched, t_star_prep1,
                                cj->hydro.super->hydro.prep1_ghost);
            scheduler_addunlock(sched, cj->hydro.super->stars.prep1_ghost,
                                t_star_prep2);
            scheduler_addunlock(sched, cj->hydro.super->hydro.prep1_ghost,
                                t_star_prep2);
            scheduler_addunlock(sched, t_star_prep2,
                                cj->hydro.super->stars.prep2_ghost);
            scheduler_addunlock(sched, cj->hydro.super->stars.prep2_ghost,
                                t_star_feedback);
#else
            scheduler_addunlock(sched, cj->hydro.super->stars.density_ghost,
                                t_star_feedback);
#endif
            scheduler_addunlock(sched, t_star_feedback,
                                cj->hydro.super->stars.stars_out);
          }
          if (with_sink) {

            /* Do the sink_formation */
            scheduler_addunlock(sched, cj->hydro.super->sinks.drift,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_in,
                                cj->top->sinks.sink_formation);
            scheduler_addunlock(sched, cj->top->sinks.sink_formation,
                                t_sink_swallow);

            /* Do the sink_swallow */
            scheduler_addunlock(sched, t_sink_swallow,
                                cj->hydro.super->sinks.sink_ghost1);

            /* Do the sink_do_gas_swallow */
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_ghost1,
                                t_sink_do_gas_swallow);
            scheduler_addunlock(sched, t_sink_do_gas_swallow,
                                cj->hydro.super->sinks.sink_ghost2);

            /* Do the sink_do_sink_swallow */
            scheduler_addunlock(sched, cj->hydro.super->sinks.sink_ghost2,
                                t_sink_do_sink_swallow);
            scheduler_addunlock(sched, t_sink_do_sink_swallow,
                                cj->hydro.super->sinks.sink_out);
          }

          if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {

            if (with_cooling)
              scheduler_addunlock(sched, cj->hydro.super->hydro.cooling_out,
                                  t_bh_density);

            scheduler_addunlock(sched, cj->hydro.super->black_holes.drift,
                                t_bh_density);
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_bh_density);
            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.black_holes_in,
                                t_bh_density);
            scheduler_addunlock(sched, t_bh_density,
                                cj->hydro.super->black_holes.density_ghost);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.density_ghost,
                                t_bh_swallow);
            scheduler_addunlock(sched, t_bh_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_1);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_1,
                                t_do_gas_swallow);
            scheduler_addunlock(sched, t_do_gas_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_2);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_2,
                                t_do_bh_swallow);
            scheduler_addunlock(sched, t_do_bh_swallow,
                                cj->hydro.super->black_holes.swallow_ghost_3);

            scheduler_addunlock(sched,
                                cj->hydro.super->black_holes.swallow_ghost_3,
                                t_bh_feedback);
            scheduler_addunlock(sched, t_bh_feedback,
                                cj->hydro.super->black_holes.black_holes_out);
          }
          if (with_rt) {
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift,
                                t_rt_gradient);
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_ghost1,
                                t_rt_gradient);
            scheduler_addunlock(sched, t_rt_gradient,
                                cj->hydro.super->rt.rt_ghost2);
            scheduler_addunlock(sched, cj->hydro.super->rt.rt_ghost2,
                                t_rt_transport);
            scheduler_addunlock(sched, t_rt_transport,
                                cj->hydro.super->rt.rt_transport_out);
          }

          if (with_timestep_limiter) {
            scheduler_addunlock(sched, cj->hydro.super->hydro.drift, t_limiter);
          }
        }

        if (ci->super != cj->super) {

          if (with_timestep_limiter) {
            scheduler_addunlock(sched, cj->super->timestep, t_limiter);
            scheduler_addunlock(sched, t_limiter, cj->super->kick1);
            scheduler_addunlock(sched, t_limiter, cj->super->timestep_limiter);
          }

          if (with_timestep_sync && with_feedback) {
            scheduler_addunlock(sched, t_star_feedback,
                                cj->super->timestep_sync);
          }
          if (with_timestep_sync && with_black_holes &&
              (bcount_i > 0 || bcount_j > 0)) {
            scheduler_addunlock(sched, t_bh_feedback, cj->super->timestep_sync);
          }
        }
      } else /* cj->nodeID != nodeID */ {
        if (with_feedback) {
#ifdef EXTRA_STAR_LOOPS
          scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                              t_star_prep1);
#endif
          scheduler_addunlock(sched, cj->hydro.super->stars.sorts,
                              t_star_feedback);
        }

        if (with_black_holes && (bcount_i > 0 || bcount_j > 0)) {
          scheduler_addunlock(sched, t_bh_swallow,
                              cj->hydro.super->black_holes.swallow_ghost_1);
        }
      }
    }
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
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_hydroloop_tasks_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;
  const int periodic = e->s->periodic;
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_stars = (e->policy & engine_policy_stars);
  const int with_sinks = (e->policy & engine_policy_sinks);
  const int with_black_holes = (e->policy & engine_policy_black_holes);

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->cdim;
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
      if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Is that neighbour local and does it have gas or star particles ? */
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

struct cell_type_pair {
  struct cell *ci, *cj;
  int type;
};

/**
 * @brief Recurse down to the super level and add a dependency between
 * rt_advance_cell_time and tend tasks. Note: This function is intended
 * for the sending side, i.e. for local cells.
 *
 * If we're running with RT subcycling, we need to ensure that nothing
 * is sent before the advance cell time task has finished. This may
 * overwrite the correct cell times, particularly so when we're sending
 * over data for non-RT tasks, e.g. for gravity pair tasks. Therefore the
 * send/tend task needs to be unlocked by the rt_advance_cell_time task.
 *
 * The send/tend task is on the top level, while the rt_advance_cell_time
 * task is on the super level. This function simply recurses down to the
 * super level and adds the required dependency.
 *
 * @param c cell to check/recurse into
 * @param tend the send/tend task that needs to be unlocked.
 * @param e the engine
 */
void engine_addunlock_rt_advance_cell_time_tend(struct cell *c,
                                                struct task *tend,
                                                struct engine *e) {

  /* safety measure */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;
  if (cell_is_empty(c)) return;

  if (c->super == c) {
    /* Found the super level cell. Add dependency from rt_advance_cell_time, if
     * it exists. */
    if (c->super->rt.rt_advance_cell_time != NULL) {
      scheduler_addunlock(&e->sched, c->super->rt.rt_advance_cell_time, tend);
    }
#ifdef SWIFT_RT_DEBUG_CHECKS
    else {
      error("Got local super cell without rt_advance_cell_time task");
    }
#endif

  } else {
    /* descend the tree until you find the super level */
    if (c->split) {
      for (int k = 0; k < 8; k++) {
        if (c->progeny[k] != NULL) {
          engine_addunlock_rt_advance_cell_time_tend(c->progeny[k], tend, e);
        }
      }
    }
  }
}

void engine_addtasks_send_mapper(void *map_data, int num_elements,
                                 void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  const int with_star_formation = (e->policy & engine_policy_star_formation);
  const int with_limiter = (e->policy & engine_policy_timestep_limiter);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_sync = (e->policy & engine_policy_timestep_sync);
  const int with_rt = (e->policy & engine_policy_rt);
  struct cell_type_pair *cell_type_pairs = (struct cell_type_pair *)map_data;

#ifdef SWIFT_DEBUG_CHECKS
  if (e->policy & engine_policy_sinks) {
    error("TODO");
  }
#endif

  for (int k = 0; k < num_elements; k++) {
    struct cell *ci = cell_type_pairs[k].ci;
    struct cell *cj = cell_type_pairs[k].cj;
    const int type = cell_type_pairs[k].type;

#ifdef WITH_MPI

    if (!cell_is_empty(ci)) {
      /* Add the timestep exchange task */
      struct task *tend = scheduler_addtask(
          &e->sched, task_type_send, task_subtype_tend, ci->mpi.tag, 0, ci, cj);
      scheduler_addunlock(&e->sched, ci->timestep_collect, tend);
      engine_addlink(e, &ci->mpi.send, tend);

      if (with_rt && (type & proxy_cell_type_hydro))
        engine_addunlock_rt_advance_cell_time_tend(ci, tend, e);
    }
#endif

    /* Add the send tasks for the cells in the proxy that have a hydro
     * connection. */
    if ((e->policy & engine_policy_hydro) && (type & proxy_cell_type_hydro))
      engine_addtasks_send_hydro(e, ci, cj, /*t_xv=*/NULL,
                                 /*t_rho=*/NULL, /*t_gradient=*/NULL,
                                 /*t_prep1=*/NULL,
                                 /*t_limiter=*/NULL, /*t_pack_limiter=*/NULL,
                                 /*t_rt_gradient=*/NULL,
                                 /*t_rt_transport=*/NULL, with_feedback,
                                 with_limiter, with_sync, with_rt);

    /* Add the send tasks for the cells in the proxy that have a stars
     * connection. */
    if ((e->policy & engine_policy_feedback) && (type & proxy_cell_type_hydro))
      engine_addtasks_send_stars(e, ci, cj, /*t_density=*/NULL,
                                 /*t_prep2=*/NULL,
                                 /*t_sf_counts=*/NULL, with_star_formation);

    /* Add the send tasks for the cells in the proxy that have a black holes
     * connection. */
    if ((e->policy & engine_policy_black_holes) &&
        (type & proxy_cell_type_hydro))
      engine_addtasks_send_black_holes(e, ci, cj, /*t_rho=*/NULL,
                                       /*t_swallow=*/NULL,
                                       /*t_gas_swallow=*/NULL,
                                       /*t_feedback=*/NULL);

    /* Add the send tasks for the cells in the proxy that have a gravity
     * connection. */
    if ((e->policy & engine_policy_self_gravity) &&
        (type & proxy_cell_type_gravity))
      engine_addtasks_send_gravity(e, ci, cj, /*t_grav=*/NULL);
  }
}

void engine_addtasks_recv_mapper(void *map_data, int num_elements,
                                 void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  const int with_star_formation = (e->policy & engine_policy_star_formation);
  const int with_limiter = (e->policy & engine_policy_timestep_limiter);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_black_holes = (e->policy & engine_policy_black_holes);
  const int with_sync = (e->policy & engine_policy_timestep_sync);
  const int with_rt = (e->policy & engine_policy_rt);
  struct cell_type_pair *cell_type_pairs = (struct cell_type_pair *)map_data;

#ifdef SWIFT_DEBUG_CHECKS
  if (e->policy & engine_policy_sinks) {
    error("TODO");
  }
#endif

  for (int k = 0; k < num_elements; k++) {
    struct cell *ci = cell_type_pairs[k].ci;
    const int type = cell_type_pairs[k].type;
    struct task *tend = NULL;

#ifdef WITH_MPI
    /* Add the timestep exchange task */
    if (!cell_is_empty(ci)) {
      tend = scheduler_addtask(&e->sched, task_type_recv, task_subtype_tend,
                               ci->mpi.tag, 0, ci, NULL);
      engine_addlink(e, &ci->mpi.recv, tend);

      /* If we're running with RT, there may be foreign cells that
       * don't receive any actual hydro particles. These cells however
       * still need to have the "advance_cell_time" and "rt_collect_times"
       * tasks in order for their time variables to be correct, especially
       * during syb-cycles, where the cell times aren't communicated at the
       * end of the step. So we create them now. */
      if (with_rt) {
#ifdef SWIFT_RT_DEBUG_CHECKS
        if (ci->top == NULL) error("Working on a cell with top == NULL??");
#endif

        /* Create the RT collect times task at the top level, if it hasn't
         * already. */
        if (ci->top->rt.rt_collect_times == NULL) {
          ci->top->rt.rt_collect_times =
              scheduler_addtask(&e->sched, task_type_rt_collect_times,
                                task_subtype_none, 0, 0, ci->top, NULL);
        }
        /* We don't need rt_collect_times -> tend dependencies. They never
         * run at the same time. rt_collect_times runs in sub-cycles,
         * tend runs on normal steps. */

        /* Make sure the timestep task replacements, i.e. rt_advance_cell_time,
         * exists on the super levels regardless of proxy type.
         * This needs to be done before engine_addtasks_recv_hydro so we
         * can set appropriate unlocks there without re-creating tasks. */
        engine_addtasks_recv_rt_advance_cell_time(e, ci, tend);
      }
    }
#endif

    /* Add the recv tasks for the cells in the proxy that have a hydro
     * connection. */
    if ((e->policy & engine_policy_hydro) && (type & proxy_cell_type_hydro)) {
      engine_addtasks_recv_hydro(
          e, ci, /*t_xv=*/NULL, /*t_rho=*/NULL, /*t_gradient=*/NULL,
          /*t_prep1=*/NULL, /*t_limiter=*/NULL, /*t_unpack_limiter=*/NULL,
          /*t_rt_gradient=*/NULL, /*t_rt_transport=*/NULL,
          /*t_rt_sorts=*/NULL, tend, with_feedback, with_black_holes,
          with_limiter, with_sync, with_rt);
    }

    /* Add the recv tasks for the cells in the proxy that have a stars
     * connection. */
    if ((e->policy & engine_policy_feedback) && (type & proxy_cell_type_hydro))
      engine_addtasks_recv_stars(e, ci, /*t_density=*/NULL, /*t_prep2=*/NULL,
                                 /*t_sf_counts=*/NULL, tend,
                                 with_star_formation);

    /* Add the recv tasks for the cells in the proxy that have a black holes
     * connection. */
    if ((e->policy & engine_policy_black_holes) &&
        (type & proxy_cell_type_hydro))
      engine_addtasks_recv_black_holes(e, ci, /*t_rho=*/NULL,
                                       /*t_swallow=*/NULL,
                                       /*t_gas_swallow=*/NULL,
                                       /*t_feedback=*/NULL, tend);

    /* Add the recv tasks for the cells in the proxy that have a gravity
     * connection. */
    if ((e->policy & engine_policy_self_gravity) &&
        (type & proxy_cell_type_gravity))
      engine_addtasks_recv_gravity(e, ci, /*t_grav=*/NULL, tend);
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
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_fofloop_tasks_mapper(void *map_data, int num_elements,
                                      void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->cdim;
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
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_fof_self, task_subtype_none, 0, 0, ci,
                        NULL);
      scheduler_addtask(sched, task_type_fof_attach_self, task_subtype_none, 0,
                        0, ci, NULL);
    }

    /* Now loop over all the neighbours of this cell */
    for (int ii = -1; ii < 2; ii++) {
      int iii = i + ii;
      if (!s->periodic && (iii < 0 || iii >= cdim[0])) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -1; jj < 2; jj++) {
        int jjj = j + jj;
        if (!s->periodic && (jjj < 0 || jjj >= cdim[1])) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -1; kk < 2; kk++) {
          int kkk = k + kk;
          if (!s->periodic && (kkk < 0 || kkk >= cdim[2])) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the neighbouring cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Does that neighbour have particles ? */
          if (cid >= cjd || cj->grav.count == 0) continue;

          /* Construct the pair search task only for fully local pairs */
          if (ci->nodeID == nodeID && cj->nodeID == nodeID)
            scheduler_addtask(sched, task_type_fof_pair, task_subtype_none, 0,
                              0, ci, cj);

          /* Construct the pair search task for pairs overlapping with the node
           */
          if (ci->nodeID == nodeID || cj->nodeID == nodeID)
            scheduler_addtask(sched, task_type_fof_attach_pair,
                              task_subtype_none, 0, 0, ci, cj);
        }
      }
    }
  }
}

/**
 * @brief Fill the #space's task list with FOF tasks.
 *
 * @param e The #engine we are working with.
 */
void engine_make_fof_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  ticks tic = getticks();

  if (e->restarting) error("Running FOF on a restart step!");

  /* Construct a FOF loop over neighbours */
  if (e->policy & engine_policy_fof)
    threadpool_map(&e->threadpool, engine_make_fofloop_tasks_mapper, NULL,
                   s->nr_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making FOF tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Split the tasks. */
  scheduler_splittasks(sched, /*fof_tasks=*/1, e->verbose);

  if (e->verbose)
    message("Splitting FOF tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that we are not left with invalid tasks */
  for (int i = 0; i < e->sched.nr_tasks; ++i) {
    const struct task *t = &e->sched.tasks[i];
    if (t->ci == NULL && t->cj != NULL && !t->skip) error("Invalid task");
  }
#endif

  /* Report the number of tasks we actually used */
  if (e->verbose)
    message(
        "Nr. of tasks: %d allocated tasks: %d ratio: %f memory use: %zd MB.",
        e->sched.nr_tasks, e->sched.size,
        (float)e->sched.nr_tasks / (float)e->sched.size,
        e->sched.size * sizeof(struct task) / (1024 * 1024));

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Fill the #space's task list.
 *
 * @param e The #engine we are working with.
 */
void engine_maketasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;
  const ticks tic = getticks();

  /* Re-set the scheduler. */
  scheduler_reset(sched, engine_estimate_nr_tasks(e));

  ticks tic2 = getticks();

  /* Construct the first hydro loop over neighbours */
  if (e->policy & engine_policy_hydro)
    threadpool_map(&e->threadpool, engine_make_hydroloop_tasks_mapper, NULL,
                   s->nr_cells, 1, threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Making hydro tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Add the self gravity tasks. */
  if (e->policy & engine_policy_self_gravity) {
    threadpool_map(&e->threadpool, engine_make_self_gravity_tasks_mapper, NULL,
                   s->nr_cells, 1, threadpool_auto_chunk_size, e);
  }

  if (e->verbose)
    message("Making gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Add the external gravity tasks. */
  if (e->policy & engine_policy_external_gravity)
    engine_make_external_gravity_tasks(e);

  if (e->sched.nr_tasks == 0 && (s->nr_gparts > 0 || s->nr_parts > 0))
    error("We have particles but no hydro or gravity tasks were created.");

  tic2 = getticks();

  /* Split the tasks. */
  scheduler_splittasks(sched, /*fof_tasks=*/0, e->verbose);

  if (e->verbose)
    message("Splitting tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that we are not left with invalid tasks */
  for (int i = 0; i < e->sched.nr_tasks; ++i) {
    const struct task *t = &e->sched.tasks[i];
    if (t->ci == NULL && t->cj != NULL && !t->skip) error("Invalid task");
  }
#endif

  /* Free the old list of cell-task links. */
  if (e->links != NULL) swift_free("links", e->links);
  e->size_links = e->sched.nr_tasks * e->links_per_tasks;

  /* Make sure that we have space for more links than last time. */
  if (e->size_links < e->nr_links * engine_rebuild_link_alloc_margin)
    e->size_links = e->nr_links * engine_rebuild_link_alloc_margin;

  /* Allocate the new link list */
  if ((e->links = (struct link *)swift_malloc(
           "links", sizeof(struct link) * e->size_links)) == NULL)
    error("Failed to allocate cell-task links.");
  e->nr_links = 0;

  tic2 = getticks();

  /* Count the number of tasks associated with each cell and
     store the density tasks in each cell, and make each sort
     depend on the sorts of its progeny. */
  threadpool_map(&e->threadpool, engine_count_and_link_tasks_mapper,
                 sched->tasks, sched->nr_tasks, sizeof(struct task),
                 threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Counting and linking tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Re-set the tag counter. MPI tags are defined for top-level cells in
   * cell_set_super_mapper. */
#ifdef WITH_MPI
  cell_next_tag = 0;
#endif

  /* Now that the self/pair tasks are at the right level, set the super
   * pointers. */
  threadpool_map(&e->threadpool, cell_set_super_mapper, cells, nr_cells,
                 sizeof(struct cell), threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("Setting super-pointers took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Append hierarchical tasks to each cell. */
  threadpool_map(&e->threadpool, engine_make_hierarchical_tasks_mapper, cells,
                 nr_cells, sizeof(struct cell), threadpool_auto_chunk_size, e);

  tic2 = getticks();

  /* Run through the tasks and make force tasks for each density task.
     Each force task depends on the cell ghosts and unlocks the kick task
     of its super-cell. */
  if (e->policy & engine_policy_hydro) {

    /* Note that this does not scale well at all so we do not use the
     * threadpool version here until the reason for this is found.
     * We call the mapper function directly as if there was only 1 thread
     * in the pool. */
    engine_make_extra_hydroloop_tasks_mapper(sched->tasks, sched->nr_tasks, e);
    /* threadpool_map(&e->threadpool, engine_make_extra_hydroloop_tasks_mapper,
     *                sched->tasks, sched->nr_tasks, sizeof(struct task),
     *                threadpool_auto_chunk_size, e); */
  }

  if (e->verbose)
    message("Making extra hydroloop tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Add the dependencies for the gravity stuff */
  if (e->policy & (engine_policy_self_gravity | engine_policy_external_gravity))
    engine_link_gravity_tasks(e);

  if (e->verbose)
    message("Linking gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

#ifdef WITH_MPI
  /* Add the communication tasks if MPI is being used. */
  if (e->policy & engine_policy_mpi) {

    tic2 = getticks();

    /* Loop over the proxies and add the send tasks, which also generates the
     * cell tags for super-cells. */
    int max_num_send_cells = 0;
    for (int pid = 0; pid < e->nr_proxies; pid++)
      max_num_send_cells += e->proxies[pid].nr_cells_out;
    struct cell_type_pair *send_cell_type_pairs = NULL;
    if ((send_cell_type_pairs = (struct cell_type_pair *)malloc(
             sizeof(struct cell_type_pair) * max_num_send_cells)) == NULL)
      error("Failed to allocate temporary cell pointer list.");
    int num_send_cells = 0;

    for (int pid = 0; pid < e->nr_proxies; pid++) {

      /* Get a handle on the proxy. */
      struct proxy *p = &e->proxies[pid];

      for (int k = 0; k < p->nr_cells_out; k++) {
        send_cell_type_pairs[num_send_cells].ci = p->cells_out[k];
        send_cell_type_pairs[num_send_cells].cj = p->cells_in[0];
        send_cell_type_pairs[num_send_cells++].type = p->cells_out_type[k];
      }
    }

    threadpool_map(&e->threadpool, engine_addtasks_send_mapper,
                   send_cell_type_pairs, num_send_cells,
                   sizeof(struct cell_type_pair), threadpool_auto_chunk_size,
                   e);

    free(send_cell_type_pairs);

    if (e->verbose)
      message("Creating send tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic2), clocks_getunit());

    tic2 = getticks();

    /* Exchange the cell tags. */
    proxy_tags_exchange(e->proxies, e->nr_proxies, s);

    if (e->verbose)
      message("Exchanging cell tags took %.3f %s.",
              clocks_from_ticks(getticks() - tic2), clocks_getunit());

    tic2 = getticks();

    /* Loop over the proxies and add the recv tasks, which relies on having the
     * cell tags. */
    int max_num_recv_cells = 0;
    for (int pid = 0; pid < e->nr_proxies; pid++)
      max_num_recv_cells += e->proxies[pid].nr_cells_in;
    struct cell_type_pair *recv_cell_type_pairs = NULL;
    if ((recv_cell_type_pairs = (struct cell_type_pair *)malloc(
             sizeof(struct cell_type_pair) * max_num_recv_cells)) == NULL)
      error("Failed to allocate temporary cell pointer list.");
    int num_recv_cells = 0;
    for (int pid = 0; pid < e->nr_proxies; pid++) {

      /* Get a handle on the proxy. */
      struct proxy *p = &e->proxies[pid];
      for (int k = 0; k < p->nr_cells_in; k++) {
        recv_cell_type_pairs[num_recv_cells].ci = p->cells_in[k];
        recv_cell_type_pairs[num_recv_cells++].type = p->cells_in_type[k];
      }
    }
    threadpool_map(&e->threadpool, engine_addtasks_recv_mapper,
                   recv_cell_type_pairs, num_recv_cells,
                   sizeof(struct cell_type_pair), threadpool_auto_chunk_size,
                   e);
    free(recv_cell_type_pairs);

    if (e->verbose)
      message("Creating recv tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic2), clocks_getunit());
  }

  /* Allocate memory for foreign particles */
  engine_allocate_foreign_particles(e, /*fof=*/0);

#endif

  /* Report the number of tasks we actually used */
  if (e->verbose)
    message(
        "Nr. of tasks: %d allocated tasks: %d ratio: %f memory use: %zd MB.",
        e->sched.nr_tasks, e->sched.size,
        (float)e->sched.nr_tasks / (float)e->sched.size,
        e->sched.size * sizeof(struct task) / (1024 * 1024));

  /* Report the number of links we actually used */
  if (e->verbose)
    message(
        "Nr. of links: %zd allocated links: %zd ratio: %f memory use: %zd MB.",
        e->nr_links, e->size_links, (float)e->nr_links / (float)e->size_links,
        e->size_links * sizeof(struct link) / (1024 * 1024));

  /* Report the values that could have been used */
  if (e->verbose)
    message("Actual usage: tasks/cell: %f links/task: %f",
            (float)e->sched.nr_tasks / s->tot_cells,
            (float)e->nr_links / e->sched.nr_tasks);

  tic2 = getticks();

  /* Set the unlocks per task. */
  scheduler_set_unlocks(sched);

  if (e->verbose)
    message("Setting unlocks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Rank the tasks. */
  scheduler_ranktasks(sched);

  if (e->verbose)
    message("Ranking the tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Weight the tasks. */
  scheduler_reweight(sched, e->verbose);

  /* Set the tasks age. */
  e->tasks_age = 0;

  if (e->verbose)
    message("took %.3f %s (including reweight).",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}
