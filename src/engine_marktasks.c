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
  const int with_star_formation = e->policy & engine_policy_star_formation;
#ifdef WITH_MPI
  const int with_feedback = e->policy & engine_policy_feedback;
#endif

  for (int ind = 0; ind < num_elements; ind++) {

    /* Get basic task information */
    struct task *t = &tasks[ind];
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

    /* Single-cell task? */
    if (t_type == task_type_self || t_type == task_type_sub_self) {

      /* Local pointer. */
      struct cell *ci = t->ci;

#ifdef SWIFT_DEBUG_CHECKS
      if (ci->nodeID != nodeID) error("Non-local self task found");
#endif
      const int ci_active_hydro = cell_is_active_hydro(ci, e);
      const int ci_active_gravity = cell_is_active_gravity(ci, e);
      const int ci_active_black_holes = cell_is_active_black_holes(ci, e);
      const int ci_active_stars = cell_is_active_stars(ci, e) ||
                                  (with_star_formation && ci_active_hydro);

      /* Activate the hydro drift */
      if (t_type == task_type_self && t_subtype == task_subtype_density) {
        if (ci_active_hydro) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
          if (with_limiter) cell_activate_limiter(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_density) {
        if (ci_active_hydro) {
          scheduler_activate(s, t);
          cell_activate_subcell_hydro_tasks(ci, NULL, s);
          if (with_limiter) cell_activate_limiter(ci, s);
        }
      }

      else if (t_type == task_type_self && t_subtype == task_subtype_force) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_force) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      else if (t->type == task_type_self &&
               t->subtype == task_subtype_limiter) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      else if (t->type == task_type_sub_self &&
               t->subtype == task_subtype_limiter) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      else if (t_type == task_type_self && t_subtype == task_subtype_gradient) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_gradient) {
        if (ci_active_hydro) scheduler_activate(s, t);
      }

      /* Activate the star density */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_stars_density) {
        if (ci_active_stars) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
          cell_activate_drift_spart(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_stars_density) {
        if (ci_active_stars) {
          scheduler_activate(s, t);
          cell_activate_subcell_stars_tasks(ci, NULL, s, with_star_formation);
        }
      }

      else if (t_type == task_type_self &&
               t_subtype == task_subtype_stars_feedback) {
        if (ci_active_stars) {
          scheduler_activate(s, t);
        }
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_stars_feedback) {
        if (ci_active_stars) scheduler_activate(s, t);
      }

      /* Activate the black hole density */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_bh_density) {
        if (ci_active_black_holes) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
          cell_activate_drift_bpart(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_bh_density) {
        if (ci_active_black_holes) {
          scheduler_activate(s, t);
          cell_activate_subcell_black_holes_tasks(ci, NULL, s);
        }
      }

      else if (t_type == task_type_self &&
               t_subtype == task_subtype_bh_swallow) {
        if (ci_active_black_holes) {
          scheduler_activate(s, t);
        }
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_bh_swallow) {
        if (ci_active_black_holes) scheduler_activate(s, t);
      }

      else if (t_type == task_type_self &&
               t_subtype == task_subtype_do_swallow) {
        if (ci_active_black_holes) {
          scheduler_activate(s, t);
        }
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_do_swallow) {
        if (ci_active_black_holes) scheduler_activate(s, t);
      }

      else if (t_type == task_type_self &&
               t_subtype == task_subtype_bh_feedback) {
        if (ci_active_black_holes) {
          scheduler_activate(s, t);
        }
      }

      else if (t_type == task_type_sub_self &&
               t_subtype == task_subtype_bh_feedback) {
        if (ci_active_black_holes) scheduler_activate(s, t);
      }

      /* Activate the gravity drift */
      else if (t_type == task_type_self && t_subtype == task_subtype_grav) {
        if (ci_active_gravity) {
          scheduler_activate(s, t);
          cell_activate_subcell_grav_tasks(t->ci, NULL, s);
        }
      }

      /* Activate the gravity drift */
      else if (t_type == task_type_self &&
               t_subtype == task_subtype_external_grav) {
        if (ci_active_gravity) {
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

      const int ci_active_black_holes = cell_is_active_black_holes(ci, e);
      const int cj_active_black_holes = cell_is_active_black_holes(cj, e);

      const int ci_active_stars = cell_is_active_stars(ci, e) ||
                                  (with_star_formation && ci_active_hydro);
      const int cj_active_stars = cell_is_active_stars(cj, e) ||
                                  (with_star_formation && cj_active_hydro);

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
      else if ((t_subtype == task_subtype_stars_density) &&
               (ci_active_stars || cj_active_stars) &&
               (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

        scheduler_activate(s, t);

        /* Set the correct sorting flags */
        if (t_type == task_type_pair) {

          /* Do ci */
          if (ci_active_stars) {

            /* stars for ci */
            atomic_or(&ci->stars.requires_sorts, 1 << t->flags);
            ci->stars.dx_max_sort_old = ci->stars.dx_max_sort;

            /* hydro for cj */
            atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
            cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

            /* Activate the drift tasks. */
            if (ci_nodeID == nodeID) cell_activate_drift_spart(ci, s);
            if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

            /* Check the sorts and activate them if needed. */
            cell_activate_hydro_sorts(cj, t->flags, s);
            cell_activate_stars_sorts(ci, t->flags, s);
          }

          /* Do cj */
          if (cj_active_stars) {

            /* hydro for ci */
            atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
            ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;

            /* stars for cj */
            atomic_or(&cj->stars.requires_sorts, 1 << t->flags);
            cj->stars.dx_max_sort_old = cj->stars.dx_max_sort;

            /* Activate the drift tasks. */
            if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
            if (cj_nodeID == nodeID) cell_activate_drift_spart(cj, s);

            /* Check the sorts and activate them if needed. */
            cell_activate_hydro_sorts(ci, t->flags, s);
            cell_activate_stars_sorts(cj, t->flags, s);
          }
        }

        /* Store current values of dx_max and h_max. */
        else if (t_type == task_type_sub_pair &&
                 t_subtype == task_subtype_stars_density) {
          cell_activate_subcell_stars_tasks(ci, cj, s, with_star_formation);
        }
      }

      /* Stars feedback */
      else if (t_subtype == task_subtype_stars_feedback) {

        /* We only want to activate the task if the cell is active and is
           going to update some gas on the *local* node */
        if ((ci_nodeID == nodeID && cj_nodeID == nodeID) &&
            (ci_active_stars || cj_active_stars)) {

          scheduler_activate(s, t);

        } else if ((ci_nodeID == nodeID && cj_nodeID != nodeID) &&
                   (cj_active_stars)) {

          scheduler_activate(s, t);

        } else if ((ci_nodeID != nodeID && cj_nodeID == nodeID) &&
                   (ci_active_stars)) {

          scheduler_activate(s, t);
        }
      }

      /* Black_Holes density */
      else if ((t_subtype == task_subtype_bh_density ||
                t_subtype == task_subtype_bh_swallow ||
                t_subtype == task_subtype_do_swallow ||
                t_subtype == task_subtype_bh_feedback) &&
               (ci_active_black_holes || cj_active_black_holes) &&
               (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

        scheduler_activate(s, t);

        /* Set the correct sorting flags */
        if (t_type == task_type_pair && t_subtype == task_subtype_bh_density) {
          if (ci_nodeID == nodeID) cell_activate_drift_bpart(ci, s);
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);

          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);
          if (cj_nodeID == nodeID) cell_activate_drift_bpart(cj, s);
        }

        /* Store current values of dx_max and h_max. */
        else if (t_type == task_type_sub_pair &&
                 t_subtype == task_subtype_bh_density) {
          cell_activate_subcell_black_holes_tasks(ci, cj, s);
        }
      }

      /* Gravity */
      else if ((t_subtype == task_subtype_grav) &&
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
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_xv);
            if (ci_active_hydro) {
              scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_recv(s, ci->mpi.recv, task_subtype_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_hydro)
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_tend_part);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_hydro) {
            struct link *l = scheduler_activate_send(
                s, cj->mpi.send, task_subtype_xv, ci_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (cj_active_hydro) {
              scheduler_activate_send(s, cj->mpi.send, task_subtype_rho,
                                      ci_nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, cj->mpi.send, task_subtype_gradient,
                                      ci_nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_hydro)
            scheduler_activate_send(s, cj->mpi.send, task_subtype_tend_part,
                                    ci_nodeID);

          /* Propagating new star counts? */
          if (with_star_formation && with_feedback) {
            if (ci_active_hydro && ci->hydro.count > 0) {
              scheduler_activate_recv(s, ci->mpi.recv, task_subtype_sf_counts);
              scheduler_activate_recv(s, ci->mpi.recv, task_subtype_tend_spart);
            }
            if (cj_active_hydro && cj->hydro.count > 0) {
              scheduler_activate_send(s, cj->mpi.send, task_subtype_sf_counts,
                                      ci_nodeID);
              scheduler_activate_send(s, cj->mpi.send, task_subtype_tend_spart,
                                      ci_nodeID);
            }
          }

        } else if (cj_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_hydro) {

            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_xv);
            if (cj_active_hydro) {
              scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_recv(s, cj->mpi.recv, task_subtype_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_hydro)
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_tend_part);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_hydro) {

            struct link *l = scheduler_activate_send(
                s, ci->mpi.send, task_subtype_xv, cj_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (ci_active_hydro) {

              scheduler_activate_send(s, ci->mpi.send, task_subtype_rho,
                                      cj_nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, ci->mpi.send, task_subtype_gradient,
                                      cj_nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_hydro)
            scheduler_activate_send(s, ci->mpi.send, task_subtype_tend_part,
                                    cj_nodeID);

          /* Propagating new star counts? */
          if (with_star_formation && with_feedback) {
            if (cj_active_hydro && cj->hydro.count > 0) {
              scheduler_activate_recv(s, cj->mpi.recv, task_subtype_sf_counts);
              scheduler_activate_recv(s, cj->mpi.recv, task_subtype_tend_spart);
            }
            if (ci_active_hydro && ci->hydro.count > 0) {
              scheduler_activate_send(s, ci->mpi.send, task_subtype_sf_counts,
                                      cj_nodeID);
              scheduler_activate_send(s, ci->mpi.send, task_subtype_tend_spart,
                                      cj_nodeID);
            }
          }
        }
#endif
      }

      /* Only interested in stars_density tasks as of here. */
      else if (t->subtype == task_subtype_stars_density) {

        /* Too much particle movement? */
        if (cell_need_rebuild_for_stars_pair(ci, cj)) *rebuild_space = 1;
        if (cell_need_rebuild_for_stars_pair(cj, ci)) *rebuild_space = 1;

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci_nodeID != nodeID) {

          if (cj_active_stars) {
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_xv);
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);

            /* If the local cell is active, more stuff will be needed. */
            scheduler_activate_send(s, cj->mpi.send, task_subtype_spart,
                                    ci_nodeID);
            cell_activate_drift_spart(cj, s);

            /* If the local cell is active, send its ti_end values. */
            scheduler_activate_send(s, cj->mpi.send, task_subtype_tend_spart,
                                    ci_nodeID);
          }

          if (ci_active_stars) {
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_spart);

            /* If the foreign cell is active, we want its ti_end values. */
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_tend_spart);

            /* Is the foreign cell active and will need stuff from us? */
            scheduler_activate_send(s, cj->mpi.send, task_subtype_xv,
                                    ci_nodeID);
            scheduler_activate_send(s, cj->mpi.send, task_subtype_rho,
                                    ci_nodeID);

            /* Drift the cell which will be sent; note that not all sent
               particles will be drifted, only those that are needed. */
            cell_activate_drift_part(cj, s);
          }

        } else if (cj_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_stars) {
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_xv);
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);

            /* If the local cell is active, more stuff will be needed. */
            scheduler_activate_send(s, ci->mpi.send, task_subtype_spart,
                                    cj_nodeID);
            cell_activate_drift_spart(ci, s);

            /* If the local cell is active, send its ti_end values. */
            scheduler_activate_send(s, ci->mpi.send, task_subtype_tend_spart,
                                    cj_nodeID);
          }

          if (cj_active_stars) {
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_spart);

            /* If the foreign cell is active, we want its ti_end values. */
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_tend_spart);

            /* Is the foreign cell active and will need stuff from us? */
            scheduler_activate_send(s, ci->mpi.send, task_subtype_xv,
                                    cj_nodeID);
            scheduler_activate_send(s, ci->mpi.send, task_subtype_rho,
                                    cj_nodeID);

            /* Drift the cell which will be sent; note that not all sent
               particles will be drifted, only those that are needed. */
            cell_activate_drift_part(ci, s);
          }
        }
#endif
      }

      /* Only interested in black hole density tasks as of here. */
      else if (t->subtype == task_subtype_bh_density) {

        /* Too much particle movement? */
        if (cell_need_rebuild_for_black_holes_pair(ci, cj)) *rebuild_space = 1;
        if (cell_need_rebuild_for_black_holes_pair(cj, ci)) *rebuild_space = 1;

        scheduler_activate(s, ci->hydro.super->black_holes.swallow_ghost[0]);
        scheduler_activate(s, cj->hydro.super->black_holes.swallow_ghost[0]);

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci_nodeID != nodeID) {

          /* Receive the foreign parts to compute BH accretion rates and do the
           * swallowing */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_part_swallow);

          /* Send the local BHs to tag the particles to swallow and to do
           * feedback */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_rho,
                                  ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_feedback,
                                  ci_nodeID);

          /* Drift before you send */
          cell_activate_drift_bpart(cj, s);

          /* Send the new BH time-steps */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_tend_bpart,
                                  ci_nodeID);

          /* Receive the foreign BHs to tag particles to swallow and for
           * feedback */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_bpart_rho);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_bpart_feedback);

          /* Receive the foreign BH time-steps */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_tend_bpart);

          /* Send the local part information */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_rho, ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_part_swallow,
                                  ci_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);

        } else if (cj_nodeID != nodeID) {

          /* Receive the foreign parts to compute BH accretion rates and do the
           * swallowing */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_part_swallow);

          /* Send the local BHs to tag the particles to swallow and to do
           * feedback */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_rho,
                                  cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_feedback,
                                  cj_nodeID);

          /* Drift before you send */
          cell_activate_drift_bpart(ci, s);

          /* Send the new BH time-steps */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_tend_bpart,
                                  cj_nodeID);

          /* Receive the foreign BHs to tag particles to swallow and for
           * feedback */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_bpart_rho);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_bpart_feedback);

          /* Receive the foreign BH time-steps */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_tend_bpart);

          /* Send the local part information */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_rho, cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_part_swallow,
                                  cj_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
        }
#endif
      }

      /* Only interested in gravity tasks as of here. */
      else if (t_subtype == task_subtype_grav) {

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (cj_active_gravity)
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_gpart);

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_gravity)
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_tend_gpart);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_gravity) {

            struct link *l = scheduler_activate_send(
                s, cj->mpi.send, task_subtype_gpart, ci_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_gravity)
            scheduler_activate_send(s, cj->mpi.send, task_subtype_tend_gpart,
                                    ci_nodeID);

        } else if (cj_nodeID != nodeID) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_gravity)
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_gpart);

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_gravity)
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_tend_gpart);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_gravity) {

            struct link *l = scheduler_activate_send(
                s, ci->mpi.send, task_subtype_gpart, cj_nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_gravity)
            scheduler_activate_send(s, ci->mpi.send, task_subtype_tend_gpart,
                                    cj_nodeID);
        }
#endif
      }
    }

    /* End force for hydro ? */
    else if (t_type == task_type_end_hydro_force) {

      if (cell_is_active_hydro(t->ci, e)) scheduler_activate(s, t);
    }

    /* End force for gravity ? */
    else if (t_type == task_type_end_grav_force) {

      if (cell_is_active_gravity(t->ci, e)) scheduler_activate(s, t);
    }

    /* Kick ? */
    else if (t_type == task_type_kick1 || t_type == task_type_kick2) {

      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e) ||
          cell_is_active_stars(t->ci, e) ||
          cell_is_active_black_holes(t->ci, e))
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
    else if (t_type == task_type_stars_ghost) {
      if (cell_is_active_stars(t->ci, e) ||
          (with_star_formation && cell_is_active_hydro(t->ci, e)))
        scheduler_activate(s, t);
    }

    /* Feedback implicit tasks? */
    else if (t_type == task_type_stars_in || t_type == task_type_stars_out) {
      if (cell_is_active_stars(t->ci, e) ||
          (with_star_formation && cell_is_active_hydro(t->ci, e)))
        scheduler_activate(s, t);
    }

    /* Black hole ghost tasks ? */
    else if (t_type == task_type_bh_density_ghost ||
             t_type == task_type_bh_swallow_ghost1 ||
             t_type == task_type_bh_swallow_ghost2) {
      if (cell_is_active_black_holes(t->ci, e)) scheduler_activate(s, t);
    }

    /* Black holes implicit tasks? */
    else if (t_type == task_type_bh_in || t_type == task_type_bh_out) {
      if (cell_is_active_black_holes(t->ci, e)) scheduler_activate(s, t);
    }

    /* Time-step? */
    else if (t_type == task_type_timestep) {
      t->ci->hydro.updated = 0;
      t->ci->grav.updated = 0;
      t->ci->stars.updated = 0;
      t->ci->black_holes.updated = 0;
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e) ||
          cell_is_active_stars(t->ci, e) ||
          cell_is_active_black_holes(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Subgrid tasks: cooling */
    else if (t_type == task_type_cooling) {
      if (cell_is_active_hydro(t->ci, e)) scheduler_activate(s, t);
    }

    /* Subgrid tasks: star formation */
    else if (t_type == task_type_star_formation) {
      if (cell_is_active_hydro(t->ci, e)) {
        scheduler_activate(s, t);
        cell_activate_super_spart_drifts(t->ci, s);
      }
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
