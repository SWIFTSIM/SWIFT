/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2026 Will J. Roper (w.roper@sussex.ac.uk)
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

/* This object's header. */
#include "scheduler.h"

/* Local headers. */
#include "engine.h"
#include "error.h"
#include "sort_part.h"
#include "space.h"
#include "space_getsid.h"
#include "threadpool.h"

/**
 * @brief Split a hydrodynamic task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_hydro(struct task *t, struct scheduler *s) {
  /* Are we considering both stars and hydro when splitting? */
  /* Note this is not very clean as the scheduler should not really
     access the engine... */
  const int with_feedback = (s->space->e->policy & engine_policy_feedback);
  const int with_stars = (s->space->e->policy & engine_policy_stars);
  const int with_sinks = (s->space->e->policy & engine_policy_sinks);
  const int with_black_holes =
      (s->space->e->policy & engine_policy_black_holes);

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Is this a non-empty self-task? */
    const int is_self =
        (t->type == task_type_self) && (t->ci != NULL) &&
        ((t->ci->hydro.count > 0) || (with_stars && t->ci->stars.count > 0) ||
         (with_sinks && t->ci->sinks.count > 0) ||
         (with_black_holes && t->ci->black_holes.count > 0));

    /* Is this a non-empty pair-task? */
    const int is_pair = (t->type == task_type_pair) && (t->ci != NULL) &&
                        (t->cj != NULL) &&
                        ((t->ci->hydro.count > 0) ||
                         (with_feedback && t->ci->stars.count > 0) ||
                         (with_sinks && t->ci->sinks.count > 0) ||
                         (with_black_holes && t->ci->black_holes.count > 0)) &&
                        ((t->cj->hydro.count > 0) ||
                         (with_feedback && t->cj->stars.count > 0) ||
                         (with_sinks && t->cj->sinks.count > 0) ||
                         (with_black_holes && t->cj->black_holes.count > 0));

    /* Empty task? */
    if (!is_self && !is_pair) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split and the task does not violate h ? */
      if (cell_can_split_self_hydro_task(ci)) {
        /* Make a sub? */
        if (scheduler_dosub && (ci->hydro.count < space_subsize_self_hydro) &&
            (ci->stars.count < space_subsize_self_stars)) {

          /* Nothing to do here */

          /* Otherwise, make tasks explicitly. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[k] != NULL &&
                (ci->progeny[k]->hydro.count ||
                 (with_stars && ci->progeny[k]->stars.count))) {
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);
            }
          }

          /* Make a task for each pair of progeny */
          for (int j = 0; j < 8; j++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[j] != NULL &&
                (ci->progeny[j]->hydro.count ||
                 (with_feedback && ci->progeny[j]->stars.count))) {
              for (int k = j + 1; k < 8; k++) {
                /* Do we have a second non-empty progenitor? */
                if (ci->progeny[k] != NULL &&
                    (ci->progeny[k]->hydro.count ||
                     (with_feedback && ci->progeny[k]->stars.count))) {
                  scheduler_splittask_hydro(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
                }
              }
            }
          }
        }
      }
    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Get the sort ID, use space_getsid_and_swap_cells and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      const int sid = space_getsid_and_swap_cells(s->space, &ci, &cj, shift);

#ifdef SWIFT_DEBUG_CHECKS
      if (sid != t->flags)
        error("Got pair task with incorrect flags: sid=%d flags=%lld", sid,
              t->flags);
#endif

      /* Should this task be split-up? */
      if (cell_can_split_pair_hydro_task(ci) &&
          cell_can_split_pair_hydro_task(cj)) {

        const int h_count_i = ci->hydro.count;
        const int h_count_j = cj->hydro.count;

        const int s_count_i = ci->stars.count;
        const int s_count_j = cj->stars.count;

        int do_sub_hydro = 1;
        int do_sub_stars_i = 1;
        int do_sub_stars_j = 1;
        if (h_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_hydro =
              h_count_i * sid_scale[sid] < space_subsize_pair_hydro / h_count_j;
        }
        if (s_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_i =
              s_count_i * sid_scale[sid] < space_subsize_pair_stars / h_count_j;
        }
        if (s_count_j > 0 && h_count_i > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_j =
              s_count_j * sid_scale[sid] < space_subsize_pair_stars / h_count_i;
        }

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            (do_sub_hydro && do_sub_stars_i && do_sub_stars_j) &&
            !sort_is_corner(sid)) {

          /* Nothing to do here! */

          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Loop over the sub-cell pairs for the current sid and add new tasks
           * for them. */
          struct cell_split_pair *csp = &cell_split_pairs[sid];

          t->ci = ci->progeny[csp->pairs[0].pid];
          t->cj = cj->progeny[csp->pairs[0].pjd];
          if (t->ci != NULL) cell_set_flag(t->ci, cell_flag_has_tasks);
          if (t->cj != NULL) cell_set_flag(t->cj, cell_flag_has_tasks);

          t->flags = csp->pairs[0].sid;
          for (int k = 1; k < csp->count; k++) {
            scheduler_splittask_hydro(
                scheduler_addtask(s, task_type_pair, t->subtype,
                                  csp->pairs[k].sid, 0,
                                  ci->progeny[csp->pairs[k].pid],
                                  cj->progeny[csp->pairs[k].pjd]),
                s);
          }
        }

        /* Otherwise, break it up if it is too large? */
      } else if (scheduler_doforcesplit && ci->split && cj->split &&
                 (ci->hydro.count > space_maxsize / cj->hydro.count)) {
        // message( "force splitting pair with %i and %i parts." ,
        // ci->hydro.count , cj->hydro.count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->hydro.count)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL && cj->progeny[k]->hydro.count) {
                struct task *tl =
                    scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k]);
                scheduler_splittask_hydro(tl, s);
                tl->flags = space_getsid_and_swap_cells(s->space, &t->ci,
                                                        &t->cj, shift);
              }
      }
    } /* pair interaction? */
  } /* iterate over the current task. */
}

/**
 * @brief Split a gravity task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_gravity(struct task *t, struct scheduler *s) {
  const struct space *sp = s->space;
  struct engine *e = sp->e;

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_pair && t->cj == NULL)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      const struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should we split this task? */
      if (cell_can_split_self_gravity_task(ci)) {
        if (scheduler_dosub && ci->grav.count < space_subsize_self_grav) {
          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              scheduler_splittask_gravity(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);

          /* Make a task for each pair of progeny */
          if (t->subtype != task_subtype_external_grav) {
            for (int j = 0; j < 8; j++) {
              if (ci->progeny[j] == NULL) continue;
              struct cell *cpi = ci->progeny[j];
              for (int k = j + 1; k < 8; k++) {
                if (ci->progeny[k] == NULL) continue;
                struct cell *cpj = ci->progeny[k];

                /* If running with the mesh this pair may be beyond the mesh
                 * criterion meaning we won't need a task here. */
                if (cell_can_use_mesh(e, cpi, cpj)) {
                  continue;
                }

                scheduler_splittask_gravity(
                    scheduler_addtask(s, task_type_pair, t->subtype,
                                      sub_sid_flag[j][k], 0, cpi, cpj),
                    s);
              } /* cpj loop */
            } /* cpi loop */
          } /* Self-gravity only */
        } /* Make tasks explicitly */
      } /* Cell is split */
    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should this task be split-up? */
      if (cell_can_split_pair_gravity_task(ci) &&
          cell_can_split_pair_gravity_task(cj)) {
        const long long gcount_i = ci->grav.count;
        const long long gcount_j = cj->grav.count;

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            gcount_i * gcount_j < ((long long)space_subsize_pair_grav)) {
          /* Otherwise, split it. */
        } else {
          /* Turn the task into a M-M task that will take care of all the
           * progeny pairs */
          t->type = task_type_grav_mm;
          t->subtype = task_subtype_none;
          t->flags = 0;

          /* Make a task for every other pair of progeny */
          for (int i = 0; i < 8; i++) {
            if (ci->progeny[i] != NULL) {
              struct cell *cpi = ci->progeny[i];
              for (int j = 0; j < 8; j++) {
                if (cj->progeny[j] != NULL) {
                  struct cell *cpj = cj->progeny[j];

                  /* If running with the mesh this pair may be beyond the mesh
                   * criterion meaning we won't need a task here. */
                  if (cell_can_use_mesh(e, cpi, cpj)) {
                    continue;
                  }

                  /* Can we use a M-M interaction here? */
                  if (cell_can_use_pair_mm(cpi, cpj, e, sp,
                                           /*use_rebuild_data=*/1,
                                           /*is_tree_walk=*/1)) {

                    /* Flag this pair as being treated by the M-M task.
                     * We use the 64 bits in the task->flags field to store
                     * this information. The corresponding task will unpack
                     * the information and operate according to the choices
                     * made here. */
                    const int flag = i * 8 + j;
                    t->flags |= (1ULL << flag);

                  } else {
                    /* Ok, we actually have to create a task */
                    scheduler_splittask_gravity(
                        scheduler_addtask(s, task_type_pair, task_subtype_grav,
                                          0, 0, cpi, cpj),
                        s);
                  }
                }
              }
            }
          }

          /* Can none of the progenies use M-M calculations? */
          if (t->flags == 0) {
            t->type = task_type_none;
            t->subtype = task_subtype_none;
            t->ci = NULL;
            t->cj = NULL;
            t->skip = 1;
          }

        } /* Split the pair */
      }
    } /* pair interaction? */
  } /* iterate over the current task. */
}

/**
 * @brief Split a FOF task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_fof(struct task *t, struct scheduler *s) {

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {

    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_fof_pair && t->cj == NULL) ||
        t->ci->grav.count == 0 || (t->cj != NULL && t->cj->grav.count == 0)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_fof_self) {

      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split? */
      if (cell_can_split_self_fof_task(ci)) {

        /* Take a step back (we're going to recycle the current task)... */
        redo = 1;

        /* Add the self tasks. */
        int first_child = 0;
        while (ci->progeny[first_child] == NULL) first_child++;
        t->ci = ci->progeny[first_child];
        for (int k = first_child + 1; k < 8; k++)
          if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
            scheduler_splittask_fof(
                scheduler_addtask(s, task_type_fof_self, t->subtype, 0, 0,
                                  ci->progeny[k], NULL),
                s);

        /* Make a task for each pair of progeny */
        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->grav.count)
            for (int k = j + 1; k < 8; k++)
              if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
                scheduler_splittask_fof(
                    scheduler_addtask(s, task_type_fof_pair, t->subtype, 0, 0,
                                      ci->progeny[j], ci->progeny[k]),
                    s);
      } /* Cell is split */

    } /* Self interaction */

  } /* iterate over the current task. */
}

/**
 * @brief Split an SIDM task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_sidm(struct task *t, struct scheduler *s) {
  /* Are we considering both stars and hydro when splitting? */
  /* Note this is not very clean as the scheduler should not really
     access the engine... */

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Is this a non-empty self-task? */
    const int is_self = (t->type == task_type_self) && (t->ci != NULL) &&
                        (t->ci->sidm.count > 0);

    /* Is this a non-empty pair-task? */
    const int is_pair = (t->type == task_type_pair) && (t->ci != NULL) &&
                        (t->cj != NULL) && (t->ci->sidm.count > 0) &&
                        (t->cj->sidm.count > 0);

    /* Empty task? */
    if (!is_self && !is_pair) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split and the task does not violate h ? */
      if (cell_can_split_self_sidm_task(ci)) {
        /* Make a sub? */
        if (scheduler_dosub && ci->sidm.count < space_subsize_self_sidm) {

          /* Nothing to do here */

          /* Otherwise, make tasks explicitly. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[k] != NULL && (ci->progeny[k]->sidm.count)) {
              scheduler_splittask_sidm(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);
            }
          }

          /* Make a task for each pair of progeny */
          for (int j = 0; j < 8; j++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[j] != NULL && ci->progeny[j]->sidm.count) {
              for (int k = j + 1; k < 8; k++) {
                /* Do we have a second non-empty progenitor? */
                if (ci->progeny[k] != NULL && ci->progeny[k]->sidm.count) {
                  scheduler_splittask_sidm(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
                }
              }
            }
          }
        }
      }
    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Get the sort ID, use space_getsid_and_swap_cells and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      const int sid = space_getsid_and_swap_cells(s->space, &ci, &cj, shift);

#ifdef SWIFT_DEBUG_CHECKS
      if (sid != t->flags)
        error("Got pair task with incorrect flags: sid=%d flags=%lld", sid,
              t->flags);
#endif

      /* Should this task be split-up? */
      if (cell_can_split_pair_sidm_task(ci) &&
          cell_can_split_pair_sidm_task(cj)) {

        const int h_count_i = ci->sidm.count;
        const int h_count_j = cj->sidm.count;

        int do_sub_sidm = 1;
        if (h_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_sidm =
              h_count_i * sid_scale[sid] < space_subsize_pair_sidm / h_count_j;
        }

        /* Replace by a single sub-task? */
        if (scheduler_dosub && do_sub_sidm
            // !sort_is_corner(sid) //TODO?
        ) {

          /* Nothing to do here! */

          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Loop over the sub-cell pairs for the current sid and add new tasks
           * for them. */
          struct cell_split_pair *csp = &cell_split_pairs[sid];

          t->ci = ci->progeny[csp->pairs[0].pid];
          t->cj = cj->progeny[csp->pairs[0].pjd];
          if (t->ci != NULL) cell_set_flag(t->ci, cell_flag_has_tasks);
          if (t->cj != NULL) cell_set_flag(t->cj, cell_flag_has_tasks);

          t->flags = csp->pairs[0].sid;
          for (int k = 1; k < csp->count; k++) {
            scheduler_splittask_sidm(
                scheduler_addtask(s, task_type_pair, t->subtype,
                                  csp->pairs[k].sid, 0,
                                  ci->progeny[csp->pairs[k].pid],
                                  cj->progeny[csp->pairs[k].pjd]),
                s);
          }
        }

        /* Otherwise, break it up if it is too large? */
      } else if (scheduler_doforcesplit && ci->split && cj->split &&
                 (ci->sidm.count > space_maxsize / cj->sidm.count)) {
        // message( "force splitting pair with %i and %i parts." ,
        // ci->hydro.count , cj->hydro.count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->sidm.count)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL && cj->progeny[k]->sidm.count) {
                struct task *tl =
                    scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k]);
                scheduler_splittask_sidm(tl, s);
                tl->flags = space_getsid_and_swap_cells(s->space, &t->ci,
                                                        &t->cj, shift);
              }
      }
    } /* pair interaction? */
  } /* iterate over the current task. */
}

/**
 * @brief Mapper function to split FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_fof_mapper(void *map_data, int num_elements,
                                     void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->type == task_type_fof_self || t->type == task_type_fof_pair) {
      scheduler_splittask_fof(t, s);
    }
  }
}

/**
 * @brief Mapper function to split non-FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_mapper(void *map_data, int num_elements,
                                 void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->subtype == task_subtype_density) {
      scheduler_splittask_hydro(t, s);
    } else if (t->subtype == task_subtype_sidm_density) {
      scheduler_splittask_sidm(t, s);      
    } else if (t->subtype == task_subtype_external_grav) {
      scheduler_splittask_gravity(t, s);
    } else if (t->subtype == task_subtype_grav) {
      scheduler_splittask_gravity(t, s);
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Unexpected task sub-type %s/%s", taskID_names[t->type],
            subtaskID_names[t->subtype]);
#endif
    }
  }
}

/**
 * @brief Splits all the tasks in the scheduler that are too large.
 *
 * @param s The #scheduler.
 * @param fof_tasks Are we splitting the FOF tasks (1)? Or the regular tasks
 * (0)?
 * @param verbose Are we talkative?
 */
void scheduler_splittasks(struct scheduler *s, const int fof_tasks,
                          const int verbose) {

  if (verbose) {
    message("space_subsize_self_hydro= %d", space_subsize_self_hydro);
    message("space_subsize_pair_hydro= %d", space_subsize_pair_hydro);
    message("space_subsize_self_stars= %d", space_subsize_self_stars);
    message("space_subsize_pair_stars= %d", space_subsize_pair_stars);
    message("space_subsize_self_grav= %d", space_subsize_self_grav);
    message("space_subsize_pair_grav= %d", space_subsize_pair_grav);
  }

  if (fof_tasks) {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_fof_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);

  } else {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);
  }
}
