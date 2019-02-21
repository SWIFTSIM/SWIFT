/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

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
#include "atomic.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "proxy.h"
#include "timers.h"

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
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.grav.send, t_grav);
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
 */
void engine_addtasks_send_hydro(struct engine *e, struct cell *ci,
                                struct cell *cj, struct task *t_xv,
                                struct task *t_rho, struct task *t_gradient) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

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
#ifdef EXTRA_HYDRO_LOOP
      t_gradient = scheduler_addtask(s, task_type_send, task_subtype_gradient,
                                     ci->mpi.tag, 0, ci, cj);
#endif

#ifdef EXTRA_HYDRO_LOOP

      scheduler_addunlock(s, t_gradient, ci->super->kick2);

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
      scheduler_addunlock(s, t_rho, ci->super->end_force);

      /* The send_rho task depends on the cell's ghost task. */
      scheduler_addunlock(s, ci->hydro.super->hydro.ghost_out, t_rho);

      /* The send_xv task should unlock the super_hydro-cell's ghost task. */
      scheduler_addunlock(s, t_xv, ci->hydro.super->hydro.ghost_in);

#endif

      /* Drift before you send */
      scheduler_addunlock(s, ci->hydro.super->hydro.drift, t_xv);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.hydro.send_xv, t_xv);
    engine_addlink(e, &ci->mpi.hydro.send_rho, t_rho);
#ifdef EXTRA_HYDRO_LOOP
    engine_addlink(e, &ci->mpi.hydro.send_gradient, t_gradient);
#endif
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_hydro(e, ci->progeny[k], cj, t_xv, t_rho,
                                   t_gradient);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the time-step to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_ti The send_ti #task, if it has already been created.
 * @param t_limiter The send_limiter #task, if already created.
 * @param with_limiter Are we running with the time-step limiter?
 */
void engine_addtasks_send_timestep(struct engine *e, struct cell *ci,
                                   struct cell *cj, struct task *t_ti,
                                   struct task *t_limiter,
                                   const int with_limiter) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Check if any of the gravity tasks are for the target node. */
  for (l = ci->grav.grav; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* Check whether instead any of the hydro tasks are for the target node. */
  if (l == NULL)
    for (l = ci->hydro.density; l != NULL; l = l->next)
      if (l->t->ci->nodeID == nodeID ||
          (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
        break;

  /* If found anything, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_ti == NULL) {

      /* Make sure this cell is tagged. */
      cell_ensure_tagged(ci);

      t_ti = scheduler_addtask(s, task_type_send, task_subtype_tend,
                               ci->mpi.tag, 0, ci, cj);

      if (with_limiter)
        t_limiter = scheduler_addtask(s, task_type_send, task_subtype_limiter,
                                      ci->mpi.tag, 0, ci, cj);

      /* The super-cell's timestep task should unlock the send_ti task. */
      scheduler_addunlock(s, ci->super->timestep, t_ti);
      if (with_limiter) scheduler_addunlock(s, t_limiter, ci->super->timestep);
      if (with_limiter)
        scheduler_addunlock(s, t_limiter, ci->super->timestep_limiter);
      if (with_limiter) scheduler_addunlock(s, ci->super->kick2, t_limiter);
      if (with_limiter)
        scheduler_addunlock(s, ci->super->timestep_limiter, t_ti);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->mpi.send_ti, t_ti);
    if (with_limiter) engine_addlink(e, &ci->mpi.limiter.send, t_limiter);
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_timestep(e, ci->progeny[k], cj, t_ti, t_limiter,
                                      with_limiter);

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
 */
void engine_addtasks_recv_hydro(struct engine *e, struct cell *c,
                                struct task *t_xv, struct task *t_rho,
                                struct task *t_gradient) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Have we reached a level where there are any hydro tasks ? */
  if (t_xv == NULL && c->hydro.density != NULL) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    /* Create the tasks. */
    t_xv = scheduler_addtask(s, task_type_recv, task_subtype_xv, c->mpi.tag, 0,
                             c, NULL);
    t_rho = scheduler_addtask(s, task_type_recv, task_subtype_rho, c->mpi.tag,
                              0, c, NULL);
#ifdef EXTRA_HYDRO_LOOP
    t_gradient = scheduler_addtask(s, task_type_recv, task_subtype_gradient,
                                   c->mpi.tag, 0, c, NULL);
#endif
  }

  c->mpi.hydro.recv_xv = t_xv;
  c->mpi.hydro.recv_rho = t_rho;
  c->mpi.hydro.recv_gradient = t_gradient;

  /* Add dependencies. */
  if (c->hydro.sorts != NULL) scheduler_addunlock(s, t_xv, c->hydro.sorts);

  for (struct link *l = c->hydro.density; l != NULL; l = l->next) {
    scheduler_addunlock(s, t_xv, l->t);
    scheduler_addunlock(s, l->t, t_rho);
  }
#ifdef EXTRA_HYDRO_LOOP
  for (struct link *l = c->hydro.gradient; l != NULL; l = l->next) {
    scheduler_addunlock(s, t_rho, l->t);
    scheduler_addunlock(s, l->t, t_gradient);
  }
  for (struct link *l = c->hydro.force; l != NULL; l = l->next)
    scheduler_addunlock(s, t_gradient, l->t);
#else
  for (struct link *l = c->hydro.force; l != NULL; l = l->next)
    scheduler_addunlock(s, t_rho, l->t);
#endif

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_hydro(e, c->progeny[k], t_xv, t_rho, t_gradient);

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
 */
void engine_addtasks_recv_gravity(struct engine *e, struct cell *c,
                                  struct task *t_grav) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

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

  c->mpi.grav.recv = t_grav;

  for (struct link *l = c->grav.grav; l != NULL; l = l->next)
    scheduler_addunlock(s, t_grav, l->t);

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_gravity(e, c->progeny[k], t_grav);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_ti The recv_ti #task, if already been created.
 * @param t_limiter The recv_limiter #task, if already created.
 * @param with_limiter Are we running with the time-step limiter?
 */
void engine_addtasks_recv_timestep(struct engine *e, struct cell *c,
                                   struct task *t_ti, struct task *t_limiter,
                                   const int with_limiter) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Have we reached a level where there are any self/pair tasks ? */
  if (t_ti == NULL && (c->grav.grav != NULL || c->hydro.density != NULL)) {

#ifdef SWIFT_DEBUG_CHECKS
    /* Make sure this cell has a valid tag. */
    if (c->mpi.tag < 0) error("Trying to receive from untagged cell.");
#endif  // SWIFT_DEBUG_CHECKS

    t_ti = scheduler_addtask(s, task_type_recv, task_subtype_tend, c->mpi.tag,
                             0, c, NULL);

    if (with_limiter)
      t_limiter = scheduler_addtask(s, task_type_recv, task_subtype_limiter,
                                    c->mpi.tag, 0, c, NULL);
  }

  c->mpi.recv_ti = t_ti;

  for (struct link *l = c->grav.grav; l != NULL; l = l->next) {
    scheduler_addunlock(s, l->t, t_ti);
  }

  if (with_limiter) {

    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_limiter);
    }

    for (struct link *l = c->hydro.limiter; l != NULL; l = l->next) {
      scheduler_addunlock(s, t_limiter, l->t);
      scheduler_addunlock(s, l->t, t_ti);
    }

  } else {

    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_addunlock(s, l->t, t_ti);
    }
  }

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_timestep(e, c->progeny[k], t_ti, t_limiter,
                                      with_limiter);

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
  const int is_with_cooling = (e->policy & engine_policy_cooling);
  const int is_with_star_formation = (e->policy & engine_policy_star_formation);
  const int with_limiter = (e->policy & engine_policy_limiter);

  /* Are we in a super-cell ? */
  if (c->super == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the two half kicks */
      c->kick1 = scheduler_addtask(s, task_type_kick1, task_subtype_none, 0, 0,
                                   c, NULL);

#if defined(WITH_LOGGER)
      c->logger = scheduler_addtask(s, task_type_logger, task_subtype_none, 0,
                                    0, c, NULL);
#endif

      c->kick2 = scheduler_addtask(s, task_type_kick2, task_subtype_none, 0, 0,
                                   c, NULL);

      /* Add the time-step calculation task and its dependency */
      c->timestep = scheduler_addtask(s, task_type_timestep, task_subtype_none,
                                      0, 0, c, NULL);

      /* Add the task finishing the force calculation */
      c->end_force = scheduler_addtask(s, task_type_end_force,
                                       task_subtype_none, 0, 0, c, NULL);

      /* Subgrid tasks */
      if (is_with_cooling && c->hydro.count_total > 0) {

        c->hydro.cooling = scheduler_addtask(s, task_type_cooling,
                                             task_subtype_none, 0, 0, c, NULL);

        scheduler_addunlock(s, c->end_force, c->hydro.cooling);
        scheduler_addunlock(s, c->hydro.cooling, c->kick2);

      } else {
        scheduler_addunlock(s, c->end_force, c->kick2);
      }

      if (is_with_star_formation && c->hydro.count_total > 0) {

        c->hydro.star_formation = scheduler_addtask(
            s, task_type_star_formation, task_subtype_none, 0, 0, c, NULL);

        scheduler_addunlock(s, c->kick2, c->hydro.star_formation);
        scheduler_addunlock(s, c->hydro.star_formation, c->timestep);

      } else {
        scheduler_addunlock(s, c->kick2, c->timestep);
      }

      scheduler_addunlock(s, c->timestep, c->kick1);

      /* Time-step limiting */
      if (with_limiter) {
        c->timestep_limiter = scheduler_addtask(
            s, task_type_timestep_limiter, task_subtype_none, 0, 0, c, NULL);

        /* Make sure it is not run before kick2 */
        scheduler_addunlock(s, c->timestep, c->timestep_limiter);
        scheduler_addunlock(s, c->timestep_limiter, c->kick1);
      }

#if defined(WITH_LOGGER)
      scheduler_addunlock(s, c->kick1, c->logger);
#endif
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
  const int periodic = e->s->periodic;
  const int is_self_gravity = (e->policy & engine_policy_self_gravity);

  /* Are we in a super-cell ? */
  if (c->grav.super == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      c->grav.drift = scheduler_addtask(s, task_type_drift_gpart,
                                        task_subtype_none, 0, 0, c, NULL);

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

        /* Gravity mesh force propagation */
        if (periodic)
          c->grav.mesh = scheduler_addtask(s, task_type_grav_mesh,
                                           task_subtype_none, 0, 0, c, NULL);

        if (periodic) scheduler_addunlock(s, c->grav.drift, c->grav.mesh);
        if (periodic) scheduler_addunlock(s, c->grav.mesh, c->grav.down);
        scheduler_addunlock(s, c->grav.init, c->grav.long_range);
        scheduler_addunlock(s, c->grav.long_range, c->grav.down);
        scheduler_addunlock(s, c->grav.down, c->super->end_force);

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
 * @brief Recursively add non-implicit star ghost tasks to a cell hierarchy.
 */
void engine_add_stars_ghosts(struct engine *e, struct cell *c,
                             struct task *stars_ghost_in,
                             struct task *stars_ghost_out) {

  /* Abort as there are no star particles here? */
  if (c->stars.count_total == 0) return;

  /* If we have reached the leaf OR have to few particles to play with*/
  if (!c->split || c->stars.count_total < engine_max_sparts_per_ghost) {

    /* Add the ghost task and its dependencies */
    struct scheduler *s = &e->sched;
    c->stars.ghost = scheduler_addtask(s, task_type_stars_ghost,
                                       task_subtype_none, 0, 0, c, NULL);
    scheduler_addunlock(s, stars_ghost_in, c->stars.ghost);
    scheduler_addunlock(s, c->stars.ghost, stars_ghost_out);
  } else {
    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_add_stars_ghosts(e, c->progeny[k], stars_ghost_in,
                                stars_ghost_out);
  }
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
 */
void engine_make_hierarchical_tasks_hydro(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;

  /* Are we in a super-cell ? */
  if (c->hydro.super == c) {

    /* Add the sort task. */
    c->hydro.sorts =
        scheduler_addtask(s, task_type_sort, task_subtype_none, 0, 0, c, NULL);

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the drift task. */
      c->hydro.drift = scheduler_addtask(s, task_type_drift_part,
                                         task_subtype_none, 0, 0, c, NULL);

      /* Generate the ghost tasks. */
      c->hydro.ghost_in =
          scheduler_addtask(s, task_type_ghost_in, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      c->hydro.ghost_out =
          scheduler_addtask(s, task_type_ghost_out, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      engine_add_ghosts(e, c, c->hydro.ghost_in, c->hydro.ghost_out);

#ifdef EXTRA_HYDRO_LOOP
      /* Generate the extra ghost task. */
      c->hydro.extra_ghost = scheduler_addtask(
          s, task_type_extra_ghost, task_subtype_none, 0, 0, c, NULL);
#endif
    }

  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_hydro(e, c->progeny[k]);
  }
}

/**
 * @brief Generate the stars hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- star version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_stars(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;

  /* Are we in a super-cell ? */
  if (c->super == c) {

    /* Add the sort task. */
    c->stars.sorts = scheduler_addtask(s, task_type_stars_sort,
                                       task_subtype_none, 0, 0, c, NULL);

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Generate the ghost tasks. */
      c->stars.ghost_in =
          scheduler_addtask(s, task_type_stars_ghost_in, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      c->stars.ghost_out =
          scheduler_addtask(s, task_type_stars_ghost_out, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      engine_add_stars_ghosts(e, c, c->stars.ghost_in, c->stars.ghost_out);
    }
  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_stars(e, c->progeny[k]);
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
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;

  /* Compute how many cells away we need to walk */
  const double distance = 2.5 * cells[0].width[0] / theta_crit;
  int delta = (int)(distance / cells[0].width[0]) + 1;
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
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

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
    for (int ii = -delta_m; ii <= delta_p; ii++) {
      int iii = i + ii;
      if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
      iii = (iii + cdim[0]) % cdim[0];
      for (int jj = -delta_m; jj <= delta_p; jj++) {
        int jjj = j + jj;
        if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
        jjj = (jjj + cdim[1]) % cdim[1];
        for (int kk = -delta_m; kk <= delta_p; kk++) {
          int kkk = k + kk;
          if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
          kkk = (kkk + cdim[2]) % cdim[2];

          /* Get the cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          struct cell *cj = &cells[cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");

          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (periodic && min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm_rebuild(ci, cj, e, s)) {

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

void engine_make_hierarchical_tasks_mapper(void *map_data, int num_elements,
                                           void *extra_data) {
  struct engine *e = (struct engine *)extra_data;
  const int is_with_hydro = (e->policy & engine_policy_hydro);
  const int is_with_self_gravity = (e->policy & engine_policy_self_gravity);
  const int is_with_external_gravity =
      (e->policy & engine_policy_external_gravity);
  const int is_with_feedback = (e->policy & engine_policy_feedback);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &((struct cell *)map_data)[ind];
    /* Make the common tasks (time integration) */
    engine_make_hierarchical_tasks_common(e, c);
    /* Add the hydro stuff */
    if (is_with_hydro) engine_make_hierarchical_tasks_hydro(e, c);
    /* And the gravity stuff */
    if (is_with_self_gravity || is_with_external_gravity)
      engine_make_hierarchical_tasks_gravity(e, c);
    if (is_with_feedback) engine_make_hierarchical_tasks_stars(e, c);
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
           finger = finger->parent)
        if (finger->hydro.sorts != NULL)
          scheduler_addunlock(sched, t, finger->hydro.sorts);
    }

    /* Link stars sort tasks to all the higher sort task. */
    if (t_type == task_type_stars_sort) {
      for (struct cell *finger = t->ci->parent; finger != NULL;
           finger = finger->parent)
        if (finger->stars.sorts != NULL)
          scheduler_addunlock(sched, t, finger->stars.sorts);
    }

    /* Link self tasks to cells. */
    else if (t_type == task_type_self) {
      atomic_inc(&ci->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->hydro.density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      } else if (t_subtype == task_subtype_external_grav) {
        engine_addlink(e, &ci->grav.grav, t);
      } else if (t->subtype == task_subtype_stars_density) {
        engine_addlink(e, &ci->stars.density, t);
      } else if (t->subtype == task_subtype_stars_feedback) {
        engine_addlink(e, &ci->stars.feedback, t);
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
      } else if (t->subtype == task_subtype_stars_density) {
        engine_addlink(e, &ci->stars.density, t);
        engine_addlink(e, &cj->stars.density, t);
      } else if (t->subtype == task_subtype_stars_feedback) {
        engine_addlink(e, &ci->stars.feedback, t);
        engine_addlink(e, &cj->stars.feedback, t);
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
      } else if (t->subtype == task_subtype_stars_density) {
        engine_addlink(e, &ci->stars.density, t);
      } else if (t->subtype == task_subtype_stars_feedback) {
        engine_addlink(e, &ci->stars.feedback, t);
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
      } else if (t->subtype == task_subtype_stars_density) {
        engine_addlink(e, &ci->stars.density, t);
        engine_addlink(e, &cj->stars.density, t);
      } else if (t->subtype == task_subtype_stars_feedback) {
        engine_addlink(e, &ci->stars.feedback, t);
        engine_addlink(e, &cj->stars.feedback, t);
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

      /* drift -----> gravity --> end_force */
      scheduler_addunlock(sched, ci->grav.super->grav.drift, t);
      scheduler_addunlock(sched, t, ci->super->end_force);
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
      scheduler_addunlock(sched, t, ci->super->end_force);
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
 * @brief Creates the dependency network for the stars tasks of a given cell.
 *
 * @param sched The #scheduler.
 * @param density The star density task to link.
 * @param feedback The star feedback task to link.
 * @param c The cell.
 */
static inline void engine_make_stars_loops_dependencies(struct scheduler *sched,
                                                        struct task *density,
                                                        struct task *feedback,
                                                        struct cell *c) {
  /* density loop --> ghost --> feedback loop*/
  scheduler_addunlock(sched, density, c->super->stars.ghost_in);
  scheduler_addunlock(sched, c->super->stars.ghost_out, feedback);
}

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
  const int with_limiter = (e->policy & engine_policy_limiter);
#ifdef EXTRA_HYDRO_LOOP
  struct task *t_gradient = NULL;
#endif
  struct task *t_force = NULL;
  struct task *t_limiter = NULL;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &((struct task *)map_data)[ind];

    /* Sort tasks depend on the drift of the cell. */
    if (t->type == task_type_sort && t->ci->nodeID == engine_rank) {
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
    }

    /* Self-interaction? */
    else if (t->type == task_type_self && t->subtype == task_subtype_density) {

      /* Make the self-density tasks depend on the drift only. */
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);

#ifdef EXTRA_HYDRO_LOOP
      /* Start by constructing the task for the second  and third hydro loop. */
      t_gradient = scheduler_addtask(sched, task_type_self,
                                     task_subtype_gradient, 0, 0, t->ci, NULL);
      t_force = scheduler_addtask(sched, task_type_self, task_subtype_force, 0,
                                  0, t->ci, NULL);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter = scheduler_addtask(sched, task_type_self,
                                      task_subtype_limiter, 0, 0, t->ci, NULL);

      /* Add the link between the new loops and the cell */
      engine_addlink(e, &t->ci->hydro.gradient, t_gradient);
      engine_addlink(e, &t->ci->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                           t_limiter, t->ci, with_cooling,
                                           with_limiter);
      scheduler_addunlock(sched, t_force, t->ci->super->end_force);
      if (with_limiter)
        scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
      if (with_limiter)
        scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
#else

      /* Start by constructing the task for the second hydro loop */
      t_force = scheduler_addtask(sched, task_type_self, task_subtype_force, 0,
                                  0, t->ci, NULL);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter = scheduler_addtask(sched, task_type_self,
                                      task_subtype_limiter, 0, 0, t->ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter, t->ci,
                                           with_cooling, with_limiter);
      scheduler_addunlock(sched, t_force, t->ci->super->end_force);
      if (with_limiter)
        scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
      if (with_limiter)
        scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
#endif
    }

    /* Otherwise, pair interaction? */
    else if (t->type == task_type_pair && t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift and the sorts. */
      if (t->ci->nodeID == engine_rank)
        scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.sorts, t);
      if (t->ci->hydro.super != t->cj->hydro.super) {
        if (t->cj->nodeID == engine_rank)
          scheduler_addunlock(sched, t->cj->hydro.super->hydro.drift, t);
        scheduler_addunlock(sched, t->cj->hydro.super->hydro.sorts, t);
      }

#ifdef EXTRA_HYDRO_LOOP
      /* Start by constructing the task for the second and third hydro loop */
      t_gradient = scheduler_addtask(sched, task_type_pair,
                                     task_subtype_gradient, 0, 0, t->ci, t->cj);
      t_force = scheduler_addtask(sched, task_type_pair, task_subtype_force, 0,
                                  0, t->ci, t->cj);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter = scheduler_addtask(sched, task_type_pair,
                                      task_subtype_limiter, 0, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->hydro.gradient, t_gradient);
      engine_addlink(e, &t->cj->hydro.gradient, t_gradient);
      engine_addlink(e, &t->ci->hydro.force, t_force);
      engine_addlink(e, &t->cj->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);
      if (with_limiter) engine_addlink(e, &t->cj->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, t->ci, with_cooling,
                                             with_limiter);
        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->hydro.super != t->cj->hydro.super) {
          engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                               t_limiter, t->cj, with_cooling,
                                               with_limiter);
        }

        if (t->ci->super != t->cj->super) {
          scheduler_addunlock(sched, t_force, t->cj->super->end_force);
          if (with_limiter)
            scheduler_addunlock(sched, t->cj->super->kick2, t_limiter);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter, t->cj->super->timestep);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter,
                                t->cj->super->timestep_limiter);
        }
      }

#else

      /* Start by constructing the task for the second hydro loop */
      t_force = scheduler_addtask(sched, task_type_pair, task_subtype_force, 0,
                                  0, t->ci, t->cj);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter = scheduler_addtask(sched, task_type_pair,
                                      task_subtype_limiter, 0, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->hydro.force, t_force);
      engine_addlink(e, &t->cj->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);
      if (with_limiter) engine_addlink(e, &t->cj->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter,
                                             t->ci, with_cooling, with_limiter);
        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->hydro.super != t->cj->hydro.super) {
          engine_make_hydro_loops_dependencies(
              sched, t, t_force, t_limiter, t->cj, with_cooling, with_limiter);
        }

        if (t->ci->super != t->cj->super) {
          scheduler_addunlock(sched, t_force, t->cj->super->end_force);
          if (with_limiter)
            scheduler_addunlock(sched, t->cj->super->kick2, t_limiter);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter, t->cj->super->timestep);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter,
                                t->cj->super->timestep_limiter);
        }
      }

#endif

    }

    /* Otherwise, sub-self interaction? */
    else if (t->type == task_type_sub_self &&
             t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift and sorts. */
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.sorts, t);

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      t_gradient =
          scheduler_addtask(sched, task_type_sub_self, task_subtype_gradient,
                            t->flags, 0, t->ci, NULL);
      t_force = scheduler_addtask(sched, task_type_sub_self, task_subtype_force,
                                  t->flags, 0, t->ci, NULL);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter =
            scheduler_addtask(sched, task_type_sub_self, task_subtype_limiter,
                              t->flags, 0, t->ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->hydro.gradient, t_gradient);
      engine_addlink(e, &t->ci->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, t->ci, with_cooling,
                                             with_limiter);
        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }

#else
      /* Start by constructing the task for the second hydro loop */
      t_force = scheduler_addtask(sched, task_type_sub_self, task_subtype_force,
                                  t->flags, 0, t->ci, NULL);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter =
            scheduler_addtask(sched, task_type_sub_self, task_subtype_limiter,
                              t->flags, 0, t->ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter,
                                             t->ci, with_cooling, with_limiter);
        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }
#endif
    }

    /* Otherwise, sub-pair interaction? */
    else if (t->type == task_type_sub_pair &&
             t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift. */
      if (t->ci->nodeID == engine_rank)
        scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.sorts, t);
      if (t->ci->hydro.super != t->cj->hydro.super) {
        if (t->cj->nodeID == engine_rank)
          scheduler_addunlock(sched, t->cj->hydro.super->hydro.drift, t);
        scheduler_addunlock(sched, t->cj->hydro.super->hydro.sorts, t);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      t_gradient =
          scheduler_addtask(sched, task_type_sub_pair, task_subtype_gradient,
                            t->flags, 0, t->ci, t->cj);
      t_force = scheduler_addtask(sched, task_type_sub_pair, task_subtype_force,
                                  t->flags, 0, t->ci, t->cj);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter =
            scheduler_addtask(sched, task_type_sub_pair, task_subtype_limiter,
                              t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->hydro.gradient, t_gradient);
      engine_addlink(e, &t->cj->hydro.gradient, t_gradient);
      engine_addlink(e, &t->ci->hydro.force, t_force);
      engine_addlink(e, &t->cj->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);
      if (with_limiter) engine_addlink(e, &t->cj->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                             t_limiter, t->ci, with_cooling,
                                             with_limiter);
        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->hydro.super != t->cj->hydro.super) {
          engine_make_hydro_loops_dependencies(sched, t, t_gradient, t_force,
                                               t_limiter, t->cj, with_cooling,
                                               with_limiter);
        }

        if (t->ci->super != t->cj->super) {
          scheduler_addunlock(sched, t_force, t->cj->super->end_force);
          if (with_limiter)
            scheduler_addunlock(sched, t->cj->super->kick2, t_limiter);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter, t->cj->super->timestep);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter,
                                t->cj->super->timestep_limiter);
        }
      }

#else
      /* Start by constructing the task for the second hydro loop */
      t_force = scheduler_addtask(sched, task_type_sub_pair, task_subtype_force,
                                  t->flags, 0, t->ci, t->cj);

      /* and the task for the time-step limiter */
      if (with_limiter)
        t_limiter =
            scheduler_addtask(sched, task_type_sub_pair, task_subtype_limiter,
                              t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->hydro.force, t_force);
      engine_addlink(e, &t->cj->hydro.force, t_force);
      if (with_limiter) engine_addlink(e, &t->ci->hydro.limiter, t_limiter);
      if (with_limiter) engine_addlink(e, &t->cj->hydro.limiter, t_limiter);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t_force, t_limiter,
                                             t->ci, with_cooling, with_limiter);

        scheduler_addunlock(sched, t_force, t->ci->super->end_force);
        if (with_limiter)
          scheduler_addunlock(sched, t->ci->super->kick2, t_limiter);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep);
        if (with_limiter)
          scheduler_addunlock(sched, t_limiter, t->ci->super->timestep_limiter);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->hydro.super != t->cj->hydro.super) {
          engine_make_hydro_loops_dependencies(
              sched, t, t_force, t_limiter, t->cj, with_cooling, with_limiter);
        }

        if (t->ci->super != t->cj->super) {
          scheduler_addunlock(sched, t_force, t->cj->super->end_force);
          if (with_limiter)
            scheduler_addunlock(sched, t->cj->super->kick2, t_limiter);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter, t->cj->super->timestep);
          if (with_limiter)
            scheduler_addunlock(sched, t_limiter,
                                t->cj->super->timestep_limiter);
        }
      }
#endif
    }
  }
}

/**
 * @brief Duplicates the first stars loop and construct all the
 * dependencies for the stars part
 *
 * This is done by looping over all the previously constructed tasks
 * and adding another task involving the same cells but this time
 * corresponding to the second stars loop over neighbours.
 * With all the relevant tasks for a given cell available, we construct
 * all the dependencies for that cell.
 */
void engine_make_extra_starsloop_tasks_mapper(void *map_data, int num_elements,
                                              void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &((struct task *)map_data)[ind];

    /* Sort tasks depend on the drift and gravity drift of the cell. */
    if (t->type == task_type_stars_sort && t->ci->nodeID == engine_rank) {
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->super->grav.drift, t);
    }

    /* Self-interaction? */
    else if (t->type == task_type_self &&
             t->subtype == task_subtype_stars_density) {

      /* Make the self-density tasks depend on the drift and gravity drift. */
      scheduler_addunlock(sched, t->ci->hydro.super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->super->grav.drift, t);

      /* Start by constructing the task for the second stars loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_self, task_subtype_stars_feedback,
                            0, 0, t->ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->stars.feedback, t2);

      /* Now, build all the dependencies for the stars */
      engine_make_stars_loops_dependencies(sched, t, t2, t->ci);

      /* end_force depends on feedback tasks */
      scheduler_addunlock(sched, t2, t->ci->super->end_force);
    }

    /* Otherwise, pair interaction? */
    else if (t->type == task_type_pair &&
             t->subtype == task_subtype_stars_density) {

      /* Make all stars density tasks depend on the hydro drift and sorts,
       * gravity drift and star sorts. */
      if (t->ci->nodeID == engine_rank)
        scheduler_addunlock(sched, t->ci->super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->super->hydro.sorts, t);
      if (t->cj->nodeID == engine_rank)
        scheduler_addunlock(sched, t->cj->super->grav.drift, t);
      scheduler_addunlock(sched, t->ci->super->stars.sorts, t);

      if (t->ci->super != t->cj->super) {
        if (t->cj->nodeID == engine_rank)
          scheduler_addunlock(sched, t->cj->super->hydro.drift, t);
        scheduler_addunlock(sched, t->cj->super->hydro.sorts, t);
        if (t->ci->nodeID == engine_rank)
          scheduler_addunlock(sched, t->ci->super->grav.drift, t);
        scheduler_addunlock(sched, t->cj->super->stars.sorts, t);
      }

      /* Start by constructing the task for the second stars loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_pair, task_subtype_stars_feedback,
                            0, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->stars.feedback, t2);
      engine_addlink(e, &t->cj->stars.feedback, t2);

      /* Now, build all the dependencies for the stars for the cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_stars_loops_dependencies(sched, t, t2, t->ci);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super != t->cj->super)
          engine_make_stars_loops_dependencies(sched, t, t2, t->cj);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t2, t->cj->super->end_force);
      }
    }

    /* Otherwise, sub-self interaction? */
    else if (t->type == task_type_sub_self &&
             t->subtype == task_subtype_stars_density) {

      /* Make all stars density tasks depend on the hydro drift and sorts,
       * gravity drift and star sorts. */
      scheduler_addunlock(sched, t->ci->super->hydro.drift, t);
      scheduler_addunlock(sched, t->ci->super->hydro.sorts, t);
      scheduler_addunlock(sched, t->ci->super->grav.drift, t);
      scheduler_addunlock(sched, t->ci->super->stars.sorts, t);

      /* Start by constructing the task for the second stars loop */
      struct task *t2 = scheduler_addtask(sched, task_type_sub_self,
                                          task_subtype_stars_feedback, t->flags,
                                          0, t->ci, t->cj);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->stars.feedback, t2);

      /* Now, build all the dependencies for the stars for the cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_stars_loops_dependencies(sched, t, t2, t->ci);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
    }

    /* Otherwise, sub-pair interaction? */
    else if (t->type == task_type_sub_pair &&
             t->subtype == task_subtype_stars_density) {

      /* Make all stars density tasks depend on the hydro drift and sorts,
       * gravity drift and star sorts. */
      if (t->cj->nodeID == engine_rank)
        scheduler_addunlock(sched, t->cj->super->hydro.drift, t);
      scheduler_addunlock(sched, t->cj->super->hydro.sorts, t);
      if (t->cj->nodeID == engine_rank)
        scheduler_addunlock(sched, t->cj->super->grav.drift, t);
      scheduler_addunlock(sched, t->ci->super->stars.sorts, t);

      if (t->ci->super != t->cj->super) {
        if (t->ci->nodeID == engine_rank)
          scheduler_addunlock(sched, t->ci->super->hydro.drift, t);
        scheduler_addunlock(sched, t->ci->super->hydro.sorts, t);
        if (t->ci->nodeID == engine_rank)
          scheduler_addunlock(sched, t->ci->super->grav.drift, t);
        scheduler_addunlock(sched, t->cj->super->stars.sorts, t);
      }

      /* Start by constructing the task for the second stars loop */
      struct task *t2 = scheduler_addtask(sched, task_type_sub_pair,
                                          task_subtype_stars_feedback, t->flags,
                                          0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->stars.feedback, t2);
      engine_addlink(e, &t->cj->stars.feedback, t2);

      /* Now, build all the dependencies for the stars for the cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_stars_loops_dependencies(sched, t, t2, t->ci);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super != t->cj->super)
          engine_make_stars_loops_dependencies(sched, t, t2, t->cj);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t2, t->cj->super->end_force);
      }
    }
  }
}

/**
 * @brief Constructs the top-level pair tasks for the star loop over
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
void engine_make_starsloop_tasks_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;
  const int periodic = e->s->periodic;

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

    /* Skip cells without particles */
    if (ci->stars.count == 0 && ci->hydro.count == 0) continue;

    /* If the cells is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_stars_density, 0, 0,
                        ci, NULL);
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

          /* Is that neighbour local and does it have particles ? */
          if (cid >= cjd || (cj->stars.count == 0 && cj->hydro.count == 0) ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Construct the pair task */
          const int sid = sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
          scheduler_addtask(sched, task_type_pair, task_subtype_stars_density,
                            sid, 0, ci, cj);
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

    /* Skip cells without hydro particles */
    if (ci->hydro.count == 0) continue;

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

          /* Is that neighbour local and does it have particles ? */
          if (cid >= cjd || cj->hydro.count == 0 ||
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

void engine_addtasks_send_mapper(void *map_data, int num_elements,
                                 void *extra_data) {
  struct engine *e = (struct engine *)extra_data;
  const int with_limiter = (e->policy & engine_policy_limiter);
  struct cell_type_pair *cell_type_pairs = (struct cell_type_pair *)map_data;

  for (int k = 0; k < num_elements; k++) {
    struct cell *ci = cell_type_pairs[k].ci;
    struct cell *cj = cell_type_pairs[k].cj;
    const int type = cell_type_pairs[k].type;

    /* Add the send task for the particle timesteps. */
    engine_addtasks_send_timestep(e, ci, cj, NULL, NULL, with_limiter);

    /* Add the send tasks for the cells in the proxy that have a hydro
     * connection. */
    if ((e->policy & engine_policy_hydro) && (type & proxy_cell_type_hydro))
      engine_addtasks_send_hydro(e, ci, cj, /*t_xv=*/NULL,
                                 /*t_rho=*/NULL, /*t_gradient=*/NULL);

    /* Add the send tasks for the cells in the proxy that have a gravity
     * connection. */
    if ((e->policy & engine_policy_self_gravity) &&
        (type & proxy_cell_type_gravity))
      engine_addtasks_send_gravity(e, ci, cj, NULL);
  }
}

void engine_addtasks_recv_mapper(void *map_data, int num_elements,
                                 void *extra_data) {
  struct engine *e = (struct engine *)extra_data;
  const int with_limiter = (e->policy & engine_policy_limiter);
  struct cell_type_pair *cell_type_pairs = (struct cell_type_pair *)map_data;

  for (int k = 0; k < num_elements; k++) {
    struct cell *ci = cell_type_pairs[k].ci;
    const int type = cell_type_pairs[k].type;

    /* Add the recv task for the particle timesteps. */
    engine_addtasks_recv_timestep(e, ci, NULL, NULL, with_limiter);

    /* Add the recv tasks for the cells in the proxy that have a hydro
     * connection. */
    if ((e->policy & engine_policy_hydro) && (type & proxy_cell_type_hydro))
      engine_addtasks_recv_hydro(e, ci, NULL, NULL, NULL);

    /* Add the recv tasks for the cells in the proxy that have a gravity
     * connection. */
    if ((e->policy & engine_policy_self_gravity) &&
        (type & proxy_cell_type_gravity))
      engine_addtasks_recv_gravity(e, ci, NULL);
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
    if (ci->nodeID == nodeID)
      scheduler_addtask(sched, task_type_fof_self, task_subtype_none, 0, 0, ci,
                        NULL);
    else
      continue;

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
                   s->nr_cells, 1, 0, e);

  /* Construct a FOF loop over neighbours */
  if (e->policy & engine_policy_fof)
    threadpool_map(&e->threadpool, engine_make_fofloop_tasks_mapper, NULL,
                   s->nr_cells, 1, 0, e);

  if (e->verbose)
    message("Making hydro tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Construct the stars hydro loop over neighbours */
  if (e->policy & engine_policy_feedback) {
    threadpool_map(&e->threadpool, engine_make_starsloop_tasks_mapper, NULL,
                   s->nr_cells, 1, 0, e);
  }

  if (e->verbose)
    message("Making stellar feedback tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Add the self gravity tasks. */
  if (e->policy & engine_policy_self_gravity) {
    threadpool_map(&e->threadpool, engine_make_self_gravity_tasks_mapper, NULL,
                   s->nr_cells, 1, 0, e);
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
  scheduler_splittasks(sched);

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
  if (e->links != NULL) free(e->links);
  e->size_links = e->sched.nr_tasks * e->links_per_tasks;

  /* Make sure that we have space for more links than last time. */
  if (e->size_links < e->nr_links * engine_rebuild_link_alloc_margin)
    e->size_links = e->nr_links * engine_rebuild_link_alloc_margin;

  /* Allocate the new link list */
  if ((e->links = (struct link *)malloc(sizeof(struct link) * e->size_links)) ==
      NULL)
    error("Failed to allocate cell-task links.");
  e->nr_links = 0;

  tic2 = getticks();

  /* Count the number of tasks associated with each cell and
     store the density tasks in each cell, and make each sort
     depend on the sorts of its progeny. */
  threadpool_map(&e->threadpool, engine_count_and_link_tasks_mapper,
                 sched->tasks, sched->nr_tasks, sizeof(struct task), 0, e);

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
                 sizeof(struct cell), 0, e);

  if (e->verbose)
    message("Setting super-pointers took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Append hierarchical tasks to each cell. */
  threadpool_map(&e->threadpool, engine_make_hierarchical_tasks_mapper, cells,
                 nr_cells, sizeof(struct cell), 0, e);

  tic2 = getticks();

  /* Run through the tasks and make force tasks for each density task.
     Each force task depends on the cell ghosts and unlocks the kick task
     of its super-cell. */
  if (e->policy & engine_policy_hydro)
    threadpool_map(&e->threadpool, engine_make_extra_hydroloop_tasks_mapper,
                   sched->tasks, sched->nr_tasks, sizeof(struct task), 0, e);

  if (e->verbose)
    message("Making extra hydroloop tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Run through the tasks and make stars feedback tasks for each stars density
     task. Each stars feedback task depends on the stars ghosts and unlocks the
     kick task of its super-cell. */
  if (e->policy & engine_policy_stars)
    threadpool_map(&e->threadpool, engine_make_extra_starsloop_tasks_mapper,
                   sched->tasks, sched->nr_tasks, sizeof(struct task), 0, e);

  if (e->verbose)
    message("Making extra starsloop tasks took %.3f %s.",
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
  if (e->policy & engine_policy_feedback)
    error("Cannot run stellar feedback with MPI (yet).");

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
                   sizeof(struct cell_type_pair),
                   /*chunk=*/0, e);

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
                   sizeof(struct cell_type_pair),
                   /*chunk=*/0, e);
    free(recv_cell_type_pairs);

    if (e->verbose)
      message("Creating recv tasks took %.3f %s.",
              clocks_from_ticks(getticks() - tic2), clocks_getunit());
  }

  /* Allocate memory for foreign particles */
  engine_allocate_foreign_particles(e);

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
