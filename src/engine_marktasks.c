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
#include "active.h"
#include "atomic.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "proxy.h"
#include "timers.h"

/**
 * @brief Mark tasks to be un-skipped and set the sort flags accordingly.
 *        Threadpool mapper function.
 *
 * @param map_data pointer to the tasks
 * @param num_elements number of tasks
 * @param extra_data pointer to int that will define if a rebuild is needed.
 */
void engine_marktasks_mapper(void *map_data, int num_elements,
                             void *extra_data) {
  /* Unpack the arguments. */
  struct task *tasks = (struct task *)map_data;
  size_t *rebuild_space = &((size_t *)extra_data)[1];
  struct scheduler *s = (struct scheduler *)(((size_t *)extra_data)[2]);
  struct engine *e = (struct engine *)((size_t *)extra_data)[0];
  const int nodeID = e->nodeID;
  const int with_limiter = e->policy & engine_policy_limiter;

  for (int ind = 0; ind < num_elements; ind++) {

    /* Get basic task information */
    struct task *t = &tasks[ind];
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

    /* Single-cell task? */
    if (t_type == task_type_self || t_type == task_type_sub_self) {

      /* Local pointer. */
      struct cell *ci = t->ci;

      if (ci->nodeID != engine_rank) error("Non-local self task found");

      /* Activate the hydro drift */
      if (t_type == task_type_self && t_subtype == task_subtype_density) {
        if (cell_is_active_hydro(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
          if (with_limiter) cell_activate_limiter(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_density) {
        if (cell_is_active_hydro(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_hydro_tasks(ci, NULL, s);
          if (with_limiter) cell_activate_limiter(ci, s);
        }
      }

      else if (t_type == task_type_self && t_subtype == task_subtype_force) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_force) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t->type == task_type_self &&
               t->subtype == task_subtype_limiter) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t->type == task_type_sub_self &&
               t->subtype == task_subtype_limiter) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

#ifdef EXTRA_HYDRO_LOOP
      else if (t_type == task_type_self && t_subtype == task_subtype_gradient) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_gradient) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }
#endif

      /* Activate the star density */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_stars_density) {
        if (cell_is_active_stars(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
          cell_activate_drift_spart(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_stars_density) {
        if (cell_is_active_stars(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_stars_tasks(ci, NULL, s);
        }
      }

      /* Activate the star feedback */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_stars_feedback) {
        if (cell_is_active_stars(ci, e)) {
          scheduler_activate(s, t);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_stars_feedback) {
        if (cell_is_active_stars(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_stars_tasks(ci, NULL, s);
        }
      }

      /* Activate the gravity drift */
      else if (t_type == task_type_self && t_subtype == task_subtype_grav) {
        if (cell_is_active_gravity(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_grav_tasks(t->ci, NULL, s);
        }
      }

      /* Activate the gravity drift */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_external_grav) {
        if (cell_is_active_gravity(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_drift_gpart(t->ci, s);
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      else {
        error("Invalid task type / sub-type encountered");
      }
#endif
    }

    /* Pair? */
    else if (t_type == task_type_pair || t_type == task_type_sub_pair) {

      /* Local pointers. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;
#ifdef WITH_MPI
      const int ci_nodeID = ci->nodeID;
      const int cj_nodeID = cj->nodeID;
#else
      const int ci_nodeID = nodeID;
      const int cj_nodeID = nodeID;
#endif
      const int ci_active_hydro = cell_is_active_hydro(ci, e);
      const int cj_active_hydro = cell_is_active_hydro(cj, e);

      const int ci_active_gravity = cell_is_active_gravity(ci, e);
      const int cj_active_gravity = cell_is_active_gravity(cj, e);

      const int ci_active_stars = cell_is_active_stars(ci, e);
      const int cj_active_stars = cell_is_active_stars(cj, e);

      /* Only activate tasks that involve a local active cell. */
      if ((t_subtype == task_subtype_density ||
           t_subtype == task_subtype_gradient ||
           t_subtype == task_subtype_limiter ||
           t_subtype == task_subtype_force) &&
          ((ci_active_hydro && ci_nodeID == nodeID) ||
           (cj_active_hydro && cj_nodeID == nodeID))) {

        scheduler_activate(s, t);

        /* Set the correct sorting flags */
        if (t_type == task_type_pair && t_subtype == task_subtype_density) {

          /* Store some values. */
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

          /* Activate the hydro drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

          /* And the limiter */
          if (ci_nodeID == nodeID && with_limiter) cell_activate_limiter(ci, s);
          if (cj_nodeID == nodeID && with_limiter) cell_activate_limiter(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
          cell_activate_hydro_sorts(cj, t->flags, s);

        }

        /* Store current values of dx_max and h_max. */
        else if (t_type == task_type_sub_pair &&
                 t_subtype == task_subtype_density) {
          cell_activate_subcell_hydro_tasks(t->ci, t->cj, s);
        }
      }

      /* Stars density */
      if (t_subtype == task_subtype_stars_density &&
          ((ci_active_stars && ci->nodeID == engine_rank) ||
           (cj_active_stars && cj->nodeID == engine_rank))) {

        // MATTHIEU: The logic here can be improved.
        // If ci is active for stars but not cj, then we can only drift the
        // stars in ci and parts in cj. (and vice-versa). The same logic can be
        // applied in cell_unskip_stars().

        scheduler_activate(s, t);

        /* Set the correct sorting flags */
        if (t_type == task_type_pair) {

          /* Do ci */
          /* Store some values. */
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          atomic_or(&ci->stars.requires_sorts, 1 << t->flags);

          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;
          ci->stars.dx_max_sort_old = ci->stars.dx_max_sort;

          /* Activate the hydro drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_spart(ci, s);

          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(cj, t->flags, s);

          cell_activate_stars_sorts(ci, t->flags, s);

          /* Do cj */
          /* Store some values. */
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          atomic_or(&cj->stars.requires_sorts, 1 << t->flags);

          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
          cj->stars.dx_max_sort_old = cj->stars.dx_max_sort;

          /* Activate the hydro drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);

          if (cj_nodeID == nodeID) cell_activate_drift_spart(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
          cell_activate_stars_sorts(cj, t->flags, s);
        }

        /* Store current values of dx_max and h_max. */
        else if (t_type == task_type_sub_pair) {
          cell_activate_subcell_stars_tasks(t->ci, t->cj, s);
        }
      }

      /* Stars feedback */
      if (t_subtype == task_subtype_stars_feedback &&
          ((ci_active_stars && ci->nodeID == engine_rank) ||
           (cj_active_stars && cj->nodeID == engine_rank))) {

        scheduler_activate(s, t);
      }

      /* Gravity */
      if ((t_subtype == task_subtype_grav) &&
          ((ci_active_gravity && ci_nodeID == nodeID) ||
           (cj_active_gravity && cj_nodeID == nodeID))) {

        scheduler_activate(s, t);

        if (t_type == task_type_pair && t_subtype == task_subtype_grav) {
          /* Activate the gravity drift */
          cell_activate_subcell_grav_tasks(t->ci, t->cj, s);
        }

#ifdef SWIFT_DEBUG_CHECKS
        else if (t_type == task_type_sub_pair &&
                 t_subtype == task_subtype_grav) {
          error("Invalid task sub-type encountered");
        }
#endif
      }

      /* Only interested in density tasks as of here. */
      if (t_subtype == task_subtype_density) {

        /* Too much particle movement? */
        if (cell_need_rebuild_for_hydro_pair(ci, cj)) *rebuild_space = 1;

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (cj_active_hydro) {
            scheduler_activate(s, ci->mpi.hydro.recv_xv);
            if (ci_active_hydro) {
              scheduler_activate(s, ci->mpi.hydro.recv_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate(s, ci->mpi.hydro.recv_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_hydro) scheduler_activate(s, ci->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_hydro) {

            struct link *l =
                scheduler_activate_send(s, cj->mpi.hydro.send_xv, ci_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (cj_active_hydro) {
              scheduler_activate_send(s, cj->mpi.hydro.send_rho, ci_nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, cj->mpi.hydro.send_gradient,
                                      ci_nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_hydro)
            scheduler_activate_send(s, cj->mpi.send_ti, ci_nodeID);

        } else if (cj_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_hydro) {
            scheduler_activate(s, cj->mpi.hydro.recv_xv);
            if (cj_active_hydro) {
              scheduler_activate(s, cj->mpi.hydro.recv_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate(s, cj->mpi.hydro.recv_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_hydro) scheduler_activate(s, cj->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_hydro) {

            struct link *l =
                scheduler_activate_send(s, ci->mpi.hydro.send_xv, cj_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (ci_active_hydro) {

              scheduler_activate_send(s, ci->mpi.hydro.send_rho, cj_nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, ci->mpi.hydro.send_gradient,
                                      cj_nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_hydro)
            scheduler_activate_send(s, ci->mpi.send_ti, cj_nodeID);
        }
#endif
      }

      /* Only interested in stars_density tasks as of here. */
      if (t->subtype == task_subtype_stars_density) {

        /* Too much particle movement? */
        if (cell_need_rebuild_for_stars_pair(ci, cj)) *rebuild_space = 1;

        // LOIC: Need implementing MPI case
      }

      /* Only interested in gravity tasks as of here. */
      if (t_subtype == task_subtype_grav) {

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (cj_active_gravity) scheduler_activate(s, ci->mpi.grav.recv);

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_gravity) scheduler_activate(s, ci->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_gravity) {

            struct link *l =
                scheduler_activate_send(s, cj->mpi.grav.send, ci_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_gravity)
            scheduler_activate_send(s, cj->mpi.send_ti, ci_nodeID);

        } else if (cj_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_gravity) scheduler_activate(s, cj->mpi.grav.recv);

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_gravity) scheduler_activate(s, cj->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_gravity) {

            struct link *l =
                scheduler_activate_send(s, ci->mpi.grav.send, cj_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_gravity)
            scheduler_activate_send(s, ci->mpi.send_ti, cj_nodeID);
        }
#endif
      }
    }

    /* End force ? */
    else if (t_type == task_type_end_force) {

      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Kick ? */
    else if (t_type == task_type_kick1 || t_type == task_type_kick2) {

      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Hydro ghost tasks ? */
    else if (t_type == task_type_ghost || t_type == task_type_extra_ghost ||
             t_type == task_type_ghost_in || t_type == task_type_ghost_out) {
      if (cell_is_active_hydro(t->ci, e)) scheduler_activate(s, t);
    }

    /* logger tasks ? */
    else if (t->type == task_type_logger) {
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e) ||
          cell_is_active_stars(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Gravity stuff ? */
    else if (t_type == task_type_grav_down || t_type == task_type_grav_mesh ||
             t_type == task_type_grav_long_range ||
             t_type == task_type_init_grav ||
             t_type == task_type_init_grav_out ||
             t_type == task_type_drift_gpart_out ||
             t_type == task_type_grav_down_in) {
      if (cell_is_active_gravity(t->ci, e)) scheduler_activate(s, t);
    }

    /* Multipole - Multipole interaction task */
    else if (t_type == task_type_grav_mm) {

      /* Local pointers. */
      const struct cell *ci = t->ci;
      const struct cell *cj = t->cj;
#ifdef WITH_MPI
      const int ci_nodeID = ci->nodeID;
      const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
      const int ci_nodeID = nodeID;
      const int cj_nodeID = nodeID;
#endif
      const int ci_active_gravity = cell_is_active_gravity_mm(ci, e);
      const int cj_active_gravity = cell_is_active_gravity_mm(cj, e);

      if ((ci_active_gravity && ci_nodeID == nodeID) ||
          (cj_active_gravity && cj_nodeID == nodeID))
        scheduler_activate(s, t);
    }

    /* Star ghost tasks ? */
    else if (t_type == task_type_stars_ghost ||
             t_type == task_type_stars_ghost_in ||
             t_type == task_type_stars_ghost_out) {
      if (cell_is_active_stars(t->ci, e)) scheduler_activate(s, t);
    }

    /* Time-step? */
    else if (t_type == task_type_timestep) {
      t->ci->hydro.updated = 0;
      t->ci->grav.updated = 0;
      t->ci->stars.updated = 0;
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Subgrid tasks */
    else if (t_type == task_type_cooling) {
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    } else if (t_type == task_type_star_formation) {
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    } else if (t_type == task_type_fof_self || t_type == task_type_fof_pair) {
      scheduler_activate(s, t);
    }
  }
}

/**
 * @brief Mark tasks to be un-skipped and set the sort flags accordingly.
 *
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */
int engine_marktasks(struct engine *e) {

  struct scheduler *s = &e->sched;
  const ticks tic = getticks();
  int rebuild_space = 0;

  /* Run through the tasks and mark as skip or not. */
  size_t extra_data[3] = {(size_t)e, (size_t)rebuild_space, (size_t)&e->sched};
  threadpool_map(&e->threadpool, engine_marktasks_mapper, s->tasks, s->nr_tasks,
                 sizeof(struct task), 0, extra_data);
  rebuild_space = extra_data[1];

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* All is well... */
  return rebuild_space;
}
