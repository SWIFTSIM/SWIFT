/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "scheduler.h"

/* Local headers. */
#include "atomic.h"
#include "cycle.h"
#include "engine.h"
#include "error.h"
#include "intrinsics.h"
#include "kernel_hydro.h"
#include "queue.h"
#include "sort_part.h"
#include "space.h"
#include "task.h"
#include "timers.h"

/**
 * @brief Re-set the list of active tasks.
 */
void scheduler_clear_active(struct scheduler *s) { s->active_count = 0; }

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param s The #scheduler.
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */
void scheduler_addunlock(struct scheduler *s, struct task *ta,
                         struct task *tb) {
#ifdef SWIFT_DEBUG_CHECKS
  if (ta == NULL) error("Unlocking task is NULL.");
  if (tb == NULL) error("Unlocked task is NULL.");
#endif

  /* Get an index at which to store this unlock. */
  const int ind = atomic_inc(&s->nr_unlocks);

  /* Does the buffer need to be grown? */
  if (ind == s->size_unlocks) {
    /* Allocate the new buffer. */
    struct task **unlocks_new;
    int *unlock_ind_new;
    const int size_unlocks_new = s->size_unlocks * 2;
    if ((unlocks_new = (struct task **)malloc(sizeof(struct task *) *
                                              size_unlocks_new)) == NULL ||
        (unlock_ind_new = (int *)malloc(sizeof(int) * size_unlocks_new)) ==
            NULL)
      error("Failed to re-allocate unlocks.");

    /* Wait for all writes to the old buffer to complete. */
    while (s->completed_unlock_writes < ind)
      ;

    /* Copy the buffers. */
    memcpy(unlocks_new, s->unlocks, sizeof(struct task *) * ind);
    memcpy(unlock_ind_new, s->unlock_ind, sizeof(int) * ind);
    free(s->unlocks);
    free(s->unlock_ind);
    s->unlocks = unlocks_new;
    s->unlock_ind = unlock_ind_new;

    /* Publish the new buffer size. */
    s->size_unlocks = size_unlocks_new;
  }

  /* Wait for there to actually be space at my index. */
  while (ind > s->size_unlocks)
    ;

  /* Write the unlock to the scheduler. */
  s->unlocks[ind] = tb;
  s->unlock_ind[ind] = ta - s->tasks;
  atomic_inc(&s->completed_unlock_writes);
}

/**
 * @brief Split a hydrodynamic task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_hydro(struct task *t, struct scheduler *s) {

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {

    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_pair && t->cj == NULL) ||
        t->ci->count == 0 || (t->cj != NULL && t->cj->count == 0)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
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
      if (cell_can_split_self_task(ci)) {

        /* Make a sub? */
        if (scheduler_dosub && ci->count < space_subsize_self) {

          /* convert to a self-subtask. */
          t->type = task_type_sub_self;

          /* Otherwise, make tasks explicitly. */
        } else {

          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;
          t->ci = ci->progeny[first_child];
          for (int k = first_child + 1; k < 8; k++)
            if (ci->progeny[k] != NULL && ci->progeny[k]->count)
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);

          /* Make a task for each pair of progeny */
          for (int j = 0; j < 8; j++)
            if (ci->progeny[j] != NULL && ci->progeny[j]->count)
              for (int k = j + 1; k < 8; k++)
                if (ci->progeny[k] != NULL && ci->progeny[k]->count)
                  scheduler_splittask_hydro(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
        }
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

      /* Get the sort ID, use space_getsid and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      const int sid = space_getsid(s->space, &ci, &cj, shift);

      /* Should this task be split-up? */
      if (cell_can_split_pair_task(ci) && cell_can_split_pair_task(cj)) {

        /* Replace by a single sub-task? */
        if (scheduler_dosub && /* Use division to avoid integer overflow. */
            ci->count * sid_scale[sid] < space_subsize_pair / cj->count &&
            !sort_is_corner(sid)) {

          /* Make this task a sub task. */
          t->type = task_type_sub_pair;

          /* Otherwise, split it. */
        } else {

          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* For each different sorting type... */
          switch (sid) {

            case 0: /* (  1 ,  1 ,  1 ) */
              t->ci = ci->progeny[7];
              t->cj = cj->progeny[0];
              t->flags = 0;
              break;

            case 1: /* (  1 ,  1 ,  0 ) */
              t->ci = ci->progeny[6];
              t->cj = cj->progeny[0];
              t->flags = 1;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[7], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[6], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              break;

            case 2: /* (  1 ,  1 , -1 ) */
              t->ci = ci->progeny[6];
              t->cj = cj->progeny[1];
              t->flags = 2;
              break;

            case 3: /* (  1 ,  0 ,  1 ) */
              t->ci = ci->progeny[5];
              t->cj = cj->progeny[0];
              t->flags = 3;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[7], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[5], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              break;

            case 4: /* (  1 ,  0 ,  0 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[0];
              t->flags = 4;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[5], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[6], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[4], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[5], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[6], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[7], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[4], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[5], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[6], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[7], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[4], cj->progeny[3]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[5], cj->progeny[3]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[6], cj->progeny[3]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[7], cj->progeny[3]),
                  s);
              break;

            case 5: /* (  1 ,  0 , -1 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[1];
              t->flags = 5;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[6], cj->progeny[3]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[4], cj->progeny[3]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[6], cj->progeny[1]),
                  s);
              break;

            case 6: /* (  1 , -1 ,  1 ) */
              t->ci = ci->progeny[5];
              t->cj = cj->progeny[2];
              t->flags = 6;
              break;

            case 7: /* (  1 , -1 ,  0 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[3];
              t->flags = 6;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[5], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[4], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[5], cj->progeny[3]),
                  s);
              break;

            case 8: /* (  1 , -1 , -1 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[3];
              t->flags = 8;
              break;

            case 9: /* (  0 ,  1 ,  1 ) */
              t->ci = ci->progeny[3];
              t->cj = cj->progeny[0];
              t->flags = 9;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[7], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[3], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              break;

            case 10: /* (  0 ,  1 ,  0 ) */
              t->ci = ci->progeny[2];
              t->cj = cj->progeny[0];
              t->flags = 10;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[3], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[6], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[2], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[3], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[6], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[7], cj->progeny[1]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[2], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[3], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[6], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[7], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[2], cj->progeny[5]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[3], cj->progeny[5]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[6], cj->progeny[5]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[7], cj->progeny[5]),
                  s);
              break;

            case 11: /* (  0 ,  1 , -1 ) */
              t->ci = ci->progeny[2];
              t->cj = cj->progeny[1];
              t->flags = 11;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[6], cj->progeny[5]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[2], cj->progeny[5]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[6], cj->progeny[1]),
                  s);
              break;

            case 12: /* (  0 ,  0 ,  1 ) */
              t->ci = ci->progeny[1];
              t->cj = cj->progeny[0];
              t->flags = 12;
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[3], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[5], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[7], cj->progeny[0]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[1], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[3], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[5], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[7], cj->progeny[2]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[1], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[3], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[5], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[7], cj->progeny[4]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[1], cj->progeny[6]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[3], cj->progeny[6]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[5], cj->progeny[6]),
                  s);
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[7], cj->progeny[6]),
                  s);
              break;
          } /* switch(sid) */
        }

        /* Otherwise, break it up if it is too large? */
      } else if (scheduler_doforcesplit && ci->split && cj->split &&
                 (ci->count > space_maxsize / cj->count)) {

        // message( "force splitting pair with %i and %i parts." , ci->count ,
        // cj->count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->count)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL && cj->progeny[k]->count) {
                struct task *tl =
                    scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k]);
                scheduler_splittask_hydro(tl, s);
                tl->flags = space_getsid(s->space, &t->ci, &t->cj, shift);
              }
      }
    } /* pair interaction? */
  }   /* iterate over the current task. */
}

/**
 * @brief Split a gravity task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_gravity(struct task *t, struct scheduler *s) {

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {

    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_pair && t->cj == NULL)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
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

      /* Should we split this task? */
      if (ci->split && ci->gcount > space_subsize_self_grav) {

        /* Take a step back (we're going to recycle the current task)... */
        redo = 1;

        /* Add the self tasks. */
        int first_child = 0;
        while (ci->progeny[first_child] == NULL) first_child++;
        t->ci = ci->progeny[first_child];
        for (int k = first_child + 1; k < 8; k++)
          if (ci->progeny[k] != NULL)
            scheduler_splittask_gravity(
                scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                  ci->progeny[k], NULL),
                s);

        /* Make a task for each pair of progeny */
        if (t->subtype != task_subtype_external_grav) {
          for (int j = 0; j < 8; j++)
            if (ci->progeny[j] != NULL)
              for (int k = j + 1; k < 8; k++)
                if (ci->progeny[k] != NULL)
                  scheduler_splittask_gravity(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
        }
      }

      /* Otherwise, make sure the self task has a drift task */
      else {

        lock_lock(&ci->lock);

        if (ci->drift_gpart == NULL)
          ci->drift_gpart = scheduler_addtask(
              s, task_type_drift_gpart, task_subtype_none, 0, 0, ci, NULL);
        lock_unlock_blind(&ci->lock);
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

      /* Should this task be split-up? */
      if (0 && ci->split && cj->split) {

        // MATTHIEU: nothing here for now

      } else {

        /* Create the drift for ci. */
        lock_lock(&ci->lock);
        if (ci->drift_gpart == NULL && ci->nodeID == engine_rank)
          ci->drift_gpart = scheduler_addtask(
              s, task_type_drift_gpart, task_subtype_none, 0, 0, ci, NULL);
        lock_unlock_blind(&ci->lock);

        /* Create the drift for cj. */
        lock_lock(&cj->lock);
        if (cj->drift_gpart == NULL && cj->nodeID == engine_rank)
          cj->drift_gpart = scheduler_addtask(
              s, task_type_drift_gpart, task_subtype_none, 0, 0, cj, NULL);
        lock_unlock_blind(&cj->lock);
      }
    } /* pair interaction? */
  }   /* iterate over the current task. */
}

/**
 * @brief Mapper function to split tasks that may be too large.
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
    } else if (t->subtype == task_subtype_external_grav) {
      scheduler_splittask_gravity(t, s);
    } else if (t->subtype == task_subtype_grav) {
      scheduler_splittask_gravity(t, s);
    } else if (t->type == task_type_grav_top_level ||
               t->type == task_type_grav_ghost) {
      // MATTHIEU: for the future
    } else {
      error("Unexpected task sub-type");
    }
  }
}

/**
 * @brief Splits all the tasks in the scheduler that are too large.
 *
 * @param s The #scheduler.
 */
void scheduler_splittasks(struct scheduler *s) {

  /* Call the mapper on each current task. */
  threadpool_map(s->threadpool, scheduler_splittasks_mapper, s->tasks,
                 s->nr_tasks, sizeof(struct task), 0, s);
}

/**
 * @brief Add a #task to the #scheduler.
 *
 * @param s The #scheduler we are working in.
 * @param type The type of the task.
 * @param subtype The sub-type of the task.
 * @param flags The flags of the task.
 * @param implicit If true, only use this task to unlock dependencies, i.e.
 *        this task is never enqueued.
 * @param ci The first cell to interact.
 * @param cj The second cell to interact.
 */
struct task *scheduler_addtask(struct scheduler *s, enum task_types type,
                               enum task_subtypes subtype, int flags,
                               int implicit, struct cell *ci, struct cell *cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == NULL && cj != NULL)
    error("Added a task with ci==NULL and cj!=NULL type=%s/%s",
          taskID_names[type], subtaskID_names[subtype]);
#endif

  /* Get the next free task. */
  const int ind = atomic_inc(&s->tasks_next);

  /* Overflow? */
  if (ind >= s->size) error("Task list overflow.");

  /* Get a pointer to the new task. */
  struct task *t = &s->tasks[ind];

  /* Copy the data. */
  t->type = type;
  t->subtype = subtype;
  t->flags = flags;
  t->wait = 0;
  t->ci = ci;
  t->cj = cj;
  t->skip = 1; /* Mark tasks as skip by default. */
  t->implicit = implicit;
  t->weight = 0;
  t->rank = 0;
  t->nr_unlock_tasks = 0;
#ifdef SWIFT_DEBUG_TASKS
  t->rid = -1;
  t->tic = 0;
  t->toc = 0;
#endif

  /* Add an index for it. */
  // lock_lock( &s->lock );
  s->tasks_ind[atomic_inc(&s->nr_tasks)] = ind;
  // lock_unlock_blind( &s->lock );

  /* Return a pointer to the new task. */
  return t;
}

/**
 * @brief Set the unlock pointers in each task.
 *
 * @param s The #scheduler.
 */
void scheduler_set_unlocks(struct scheduler *s) {

  /* Store the counts for each task. */
  short int *counts;
  if ((counts = (short int *)malloc(sizeof(short int) * s->nr_tasks)) == NULL)
    error("Failed to allocate temporary counts array.");
  bzero(counts, sizeof(short int) * s->nr_tasks);
  for (int k = 0; k < s->nr_unlocks; k++) {
    counts[s->unlock_ind[k]] += 1;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we are not overflowing */
    if (counts[s->unlock_ind[k]] < 0)
      error("Task unlocking more than %d other tasks!",
            (1 << (8 * sizeof(short int) - 1)) - 1);

#endif
  }

  /* Compute the offset for each unlock block. */
  int *offsets;
  if ((offsets = (int *)malloc(sizeof(int) * (s->nr_tasks + 1))) == NULL)
    error("Failed to allocate temporary offsets array.");
  offsets[0] = 0;
  for (int k = 0; k < s->nr_tasks; k++) {
    offsets[k + 1] = offsets[k] + counts[k];

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we are not overflowing */
    if (offsets[k + 1] < 0) error("Task unlock offset array overflowing");
#endif
  }

  /* Create and fill a temporary array with the sorted unlocks. */
  struct task **unlocks;
  if ((unlocks = (struct task **)malloc(sizeof(struct task *) *
                                        s->size_unlocks)) == NULL)
    error("Failed to allocate temporary unlocks array.");
  for (int k = 0; k < s->nr_unlocks; k++) {
    const int ind = s->unlock_ind[k];
    unlocks[offsets[ind]] = s->unlocks[k];
    offsets[ind] += 1;
  }

  /* Swap the unlocks. */
  free(s->unlocks);
  s->unlocks = unlocks;

  /* Re-set the offsets. */
  offsets[0] = 0;
  for (int k = 1; k < s->nr_tasks; k++)
    offsets[k] = offsets[k - 1] + counts[k - 1];

  /* Set the unlocks in the tasks. */
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &s->tasks[k];
    t->nr_unlock_tasks = counts[k];
    t->unlock_tasks = &s->unlocks[offsets[k]];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that there are no duplicate unlocks. */
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &s->tasks[k];
    for (int i = 0; i < t->nr_unlock_tasks; i++) {
      for (int j = i + 1; j < t->nr_unlock_tasks; j++) {
        if (t->unlock_tasks[i] == t->unlock_tasks[j])
          error("duplicate unlock! t->type=%s/%s unlocking type=%s/%s",
                taskID_names[t->type], subtaskID_names[t->subtype],
                taskID_names[t->unlock_tasks[i]->type],
                subtaskID_names[t->unlock_tasks[i]->subtype]);
      }
    }
  }
#endif

  /* Clean up. */
  free(counts);
  free(offsets);
}

/**
 * @brief Sort the tasks in topological order over all queues.
 *
 * @param s The #scheduler.
 */
void scheduler_ranktasks(struct scheduler *s) {

  struct task *tasks = s->tasks;
  int *tid = s->tasks_ind;
  const int nr_tasks = s->nr_tasks;

  /* Run through the tasks and get all the waits right. */
  for (int i = 0; i < nr_tasks; i++) {
    struct task *t = &tasks[i];

    // Increment the waits of the dependances
    for (int k = 0; k < t->nr_unlock_tasks; k++) {
      t->unlock_tasks[k]->wait++;
    }
  }

  /* Load the tids of tasks with no waits. */
  int left = 0;
  for (int k = 0; k < nr_tasks; k++)
    if (tasks[k].wait == 0) {
      tid[left] = k;
      left += 1;
    }

  /* Main loop. */
  for (int j = 0, rank = 0; j < nr_tasks; rank++) {

    /* Did we get anything? */
    if (j == left) error("Unsatisfiable task dependencies detected.");

    /* Unlock the next layer of tasks. */
    const int left_old = left;
    for (; j < left_old; j++) {
      struct task *t = &tasks[tid[j]];
      t->rank = rank;
      /* message( "task %i of type %s has rank %i." , i ,
          (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ?
         "pair" : "sort" , rank ); */
      for (int k = 0; k < t->nr_unlock_tasks; k++) {
        struct task *u = t->unlock_tasks[k];
        if (--u->wait == 0) {
          tid[left] = u - tasks;
          left += 1;
        }
      }
    }

    /* Move back to the old left (like Sanders!). */
    j = left_old;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the tasks were ranked correctly. */
  for (int k = 1; k < s->nr_tasks; k++)
    if (tasks[tid[k - 1]].rank > tasks[tid[k]].rank)
      error("Task ranking failed.");
#endif
}

/**
 * @brief (Re)allocate the task arrays.
 *
 * @param s The #scheduler.
 * @param size The maximum number of tasks in the #scheduler.
 */
void scheduler_reset(struct scheduler *s, int size) {

  /* Do we need to re-allocate? */
  if (size > s->size) {

    /* Free existing task lists if necessary. */
    scheduler_free_tasks(s);

    /* Allocate the new lists. */
    if (posix_memalign((void *)&s->tasks, task_align,
                       size * sizeof(struct task)) != 0)
      error("Failed to allocate task array.");

    if ((s->tasks_ind = (int *)malloc(sizeof(int) * size)) == NULL)
      error("Failed to allocate task lists.");

    if ((s->tid_active = (int *)malloc(sizeof(int) * size)) == NULL)
      error("Failed to allocate aactive task lists.");
  }

  /* Reset the counters. */
  s->size = size;
  s->nr_tasks = 0;
  s->tasks_next = 0;
  s->waiting = 0;
  s->nr_unlocks = 0;
  s->completed_unlock_writes = 0;
  s->active_count = 0;

  /* Set the task pointers in the queues. */
  for (int k = 0; k < s->nr_queues; k++) s->queues[k].tasks = s->tasks;
}

/**
 * @brief Compute the task weights
 *
 * @param s The #scheduler.
 * @param verbose Are we talkative?
 */
void scheduler_reweight(struct scheduler *s, int verbose) {

  const int nr_tasks = s->nr_tasks;
  int *tid = s->tasks_ind;
  struct task *tasks = s->tasks;
  const int nodeID = s->nodeID;
  const float wscale = 0.001;
  const ticks tic = getticks();

  /* Run through the tasks backwards and set their weights. */
  for (int k = nr_tasks - 1; k >= 0; k--) {
    struct task *t = &tasks[tid[k]];
    t->weight = 0;
    for (int j = 0; j < t->nr_unlock_tasks; j++)
      if (t->unlock_tasks[j]->weight > t->weight)
        t->weight = t->unlock_tasks[j]->weight;
    int cost = 0;
    switch (t->type) {
      case task_type_sort:
        cost = wscale * intrinsics_popcount(t->flags) * t->ci->count *
               (sizeof(int) * 8 - intrinsics_clz(t->ci->count));
        break;
      case task_type_self:
        cost = 1 * wscale * t->ci->count * t->ci->count;
        break;
      case task_type_pair:
        if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID)
          cost = 3 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
        else
          cost = 2 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
        break;
      case task_type_sub_pair:
        if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID) {
          if (t->flags < 0)
            cost = 3 * wscale * t->ci->count * t->cj->count;
          else
            cost =
                3 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
        } else {
          if (t->flags < 0)
            cost = 2 * wscale * t->ci->count * t->cj->count;
          else
            cost =
                2 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
        }
        break;
      case task_type_sub_self:
        cost = 1 * wscale * t->ci->count * t->ci->count;
        break;
      case task_type_ghost:
        if (t->ci == t->ci->super) cost = wscale * t->ci->count;
        break;
      case task_type_drift_part:
        cost = wscale * t->ci->count;
        break;
      case task_type_drift_gpart:
        cost = wscale * t->ci->gcount;
        break;
      case task_type_kick1:
        cost = wscale * t->ci->count;
        break;
      case task_type_kick2:
        cost = wscale * t->ci->count;
        break;
      case task_type_timestep:
        cost = wscale * t->ci->count;
        break;
      default:
        cost = 0;
        break;
    }
#if defined(WITH_MPI) && defined(HAVE_METIS)
    t->cost = cost;
#endif
    t->weight += cost;
  }

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* int min = tasks[0].weight, max = tasks[0].weight;
  for ( int k = 1 ; k < nr_tasks ; k++ )
      if ( tasks[k].weight < min )
          min = tasks[k].weight;
      else if ( tasks[k].weight > max )
          max = tasks[k].weight;
  message( "task weights are in [ %i , %i ]." , min , max ); */
}

/**
 * @brief #threadpool_map function which runs through the task
 *        graph and re-computes the task wait counters.
 */
void scheduler_rewait_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  struct scheduler *s = (struct scheduler *)extra_data;
  const int *tid = (int *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &s->tasks[tid[ind]];

    /* Ignore skipped tasks. */
    if (t->skip) continue;

    /* Increment the task's own wait counter for the enqueueing. */
    atomic_inc(&t->wait);

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that we don't have more waits that what can be stored. */
    if (t->wait < 0)
      error("Task unlocked by more than %d tasks!",
            (1 << (8 * sizeof(t->wait) - 1)) - 1);
#endif

    /* Sets the waits of the dependances */
    for (int k = 0; k < t->nr_unlock_tasks; k++) {
      struct task *u = t->unlock_tasks[k];
      atomic_inc(&u->wait);
    }
  }
}

void scheduler_enqueue_mapper(void *map_data, int num_elements,
                              void *extra_data) {
  struct scheduler *s = (struct scheduler *)extra_data;
  const int *tid = (int *)map_data;
  struct task *tasks = s->tasks;
  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[tid[ind]];
    if (atomic_dec(&t->wait) == 1 && !t->skip) {
      scheduler_enqueue(s, t);
    }
  }
  pthread_cond_broadcast(&s->sleep_cond);
}

/**
 * @brief Start the scheduler, i.e. fill the queues with ready tasks.
 *
 * @param s The #scheduler.
 */
void scheduler_start(struct scheduler *s) {

/* Reset all task debugging timers */
#ifdef SWIFT_DEBUG_TASKS
  for (int i = 0; i < s->nr_tasks; ++i) {
    s->tasks[i].tic = 0;
    s->tasks[i].toc = 0;
    s->tasks[i].rid = -1;
  }
#endif

  /* Re-wait the tasks. */
  if (s->active_count > 1000) {
    threadpool_map(s->threadpool, scheduler_rewait_mapper, s->tid_active,
                   s->active_count, sizeof(int), 0, s);
  } else {
    scheduler_rewait_mapper(s->tid_active, s->active_count, s);
  }

/* Check we have not missed an active task */
#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_current = s->space->e->ti_current;

  if (ti_current > 0) {

    for (int k = 0; k < s->nr_tasks; k++) {

      struct task *t = &s->tasks[k];
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      if (t->type == task_type_none) continue;

      /* Don't check MPI stuff */
      if (t->type == task_type_send || t->type == task_type_recv) continue;

      /* Don't check the FFT task */
      if (t->type == task_type_grav_top_level ||
          t->type == task_type_grav_ghost)
        continue;

      if (ci == NULL && cj == NULL) {

        error("Task not associated with cells!");

      } else if (cj == NULL) { /* self */

        if (ci->ti_end_min == ti_current && t->skip &&
            t->type != task_type_sort && t->type != task_type_drift_part &&
            t->type != task_type_drift_gpart)
          error(
              "Task (type='%s/%s') should not have been skipped "
              "ti_current=%lld "
              "c->ti_end_min=%lld",
              taskID_names[t->type], subtaskID_names[t->subtype], ti_current,
              ci->ti_end_min);

        /* Special treatment for sort tasks */
        /* if (ci->ti_end_min == ti_current && t->skip &&
            t->type == task_type_sort && t->flags == 0)
          error(
              "Task (type='%s/%s') should not have been skipped "
              "ti_current=%lld "
              "c->ti_end_min=%lld t->flags=%d",
              taskID_names[t->type], subtaskID_names[t->subtype], ti_current,
              ci->ti_end_min, t->flags); */

      } else { /* pair */

        if (t->skip) {

          /* Check that the pair is active if the local cell is active */
          if ((ci->ti_end_min == ti_current && ci->nodeID == engine_rank) ||
              (cj->ti_end_min == ti_current && cj->nodeID == engine_rank))
            error(
                "Task (type='%s/%s') should not have been skipped "
                "ti_current=%lld "
                "ci->ti_end_min=%lld cj->ti_end_min=%lld",
                taskID_names[t->type], subtaskID_names[t->subtype], ti_current,
                ci->ti_end_min, cj->ti_end_min);
        }
      }
    }
  }
#endif

  /* Loop over the tasks and enqueue whoever is ready. */
  if (s->active_count > 1000) {
    threadpool_map(s->threadpool, scheduler_enqueue_mapper, s->tid_active,
                   s->active_count, sizeof(int), 0, s);
  } else {
    scheduler_enqueue_mapper(s->tid_active, s->active_count, s);
  }

  /* Clear the list of active tasks. */
  s->active_count = 0;

  /* To be safe, fire of one last sleep_cond in a safe way. */
  pthread_mutex_lock(&s->sleep_mutex);
  pthread_cond_broadcast(&s->sleep_cond);
  pthread_mutex_unlock(&s->sleep_mutex);
}

/**
 * @brief Put a task on one of the queues.
 *
 * @param s The #scheduler.
 * @param t The #task.
 */
void scheduler_enqueue(struct scheduler *s, struct task *t) {

  /* The target queue for this task. */
  int qid = -1;

  /* Ignore skipped tasks */
  if (t->skip) return;

  /* If this is an implicit task, just pretend it's done. */
  if (t->implicit) {
#ifdef SWIFT_DEBUG_CHECKS
    t->ti_run = s->space->e->ti_current;
#endif
    t->skip = 1;
    for (int j = 0; j < t->nr_unlock_tasks; j++) {
      struct task *t2 = t->unlock_tasks[j];
      if (atomic_dec(&t2->wait) == 1) scheduler_enqueue(s, t2);
    }
  }

  /* Otherwise, look for a suitable queue. */
  else {
#ifdef WITH_MPI
    int err = MPI_SUCCESS;
#endif

    /* Find the previous owner for each task type, and do
       any pre-processing needed. */
    switch (t->type) {
      case task_type_self:
      case task_type_sub_self:
      case task_type_sort:
      case task_type_ghost:
      case task_type_kick1:
      case task_type_kick2:
      case task_type_drift_part:
      case task_type_drift_gpart:
      case task_type_timestep:
        qid = t->ci->super->owner;
        break;
      case task_type_pair:
      case task_type_sub_pair:
        qid = t->ci->super->owner;
        if (qid < 0 ||
            s->queues[qid].count > s->queues[t->cj->super->owner].count)
          qid = t->cj->super->owner;
        break;
      case task_type_recv:
#ifdef WITH_MPI
        if (t->subtype == task_subtype_tend) {
          t->buff = malloc(sizeof(integertime_t) * t->ci->pcell_size);
          err = MPI_Irecv(t->buff, t->ci->pcell_size * sizeof(integertime_t),
                          MPI_BYTE, t->ci->nodeID, t->flags, MPI_COMM_WORLD,
                          &t->req);
        } else if (t->subtype == task_subtype_xv ||
                   t->subtype == task_subtype_rho ||
                   t->subtype == task_subtype_gradient) {
          err = MPI_Irecv(t->ci->parts, t->ci->count, part_mpi_type,
                          t->ci->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
          // message( "receiving %i parts with tag=%i from %i to %i." ,
          //     t->ci->count , t->flags , t->ci->nodeID , s->nodeID );
          // fflush(stdout);
        } else if (t->subtype == task_subtype_gpart) {
          err = MPI_Irecv(t->ci->gparts, t->ci->gcount, gpart_mpi_type,
                          t->ci->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else if (t->subtype == task_subtype_spart) {
          err = MPI_Irecv(t->ci->sparts, t->ci->scount, spart_mpi_type,
                          t->ci->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else if (t->subtype == task_subtype_multipole) {
          err = MPI_Irecv(t->ci->multipole, 1, multipole_mpi_type,
                          t->ci->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else {
          error("Unknown communication sub-type");
        }
        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit irecv for particle data.");
        }
        qid = 1 % s->nr_queues;
#else
        error("SWIFT was not compiled with MPI support.");
#endif
        break;
      case task_type_send:
#ifdef WITH_MPI
        if (t->subtype == task_subtype_tend) {
          t->buff = malloc(sizeof(integertime_t) * t->ci->pcell_size);
          cell_pack_ti_ends(t->ci, t->buff);
          err = MPI_Isend(t->buff, t->ci->pcell_size * sizeof(integertime_t),
                          MPI_BYTE, t->cj->nodeID, t->flags, MPI_COMM_WORLD,
                          &t->req);
        } else if (t->subtype == task_subtype_xv ||
                   t->subtype == task_subtype_rho ||
                   t->subtype == task_subtype_gradient) {
          err = MPI_Isend(t->ci->parts, t->ci->count, part_mpi_type,
                          t->cj->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
          // message( "sending %i parts with tag=%i from %i to %i." ,
          //     t->ci->count , t->flags , s->nodeID , t->cj->nodeID );
          // fflush(stdout);
        } else if (t->subtype == task_subtype_gpart) {
          err = MPI_Isend(t->ci->gparts, t->ci->gcount, gpart_mpi_type,
                          t->cj->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else if (t->subtype == task_subtype_spart) {
          err = MPI_Isend(t->ci->sparts, t->ci->scount, spart_mpi_type,
                          t->cj->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else if (t->subtype == task_subtype_multipole) {
          err = MPI_Isend(t->ci->multipole, 1, multipole_mpi_type,
                          t->cj->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        } else {
          error("Unknown communication sub-type");
        }
        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit isend for particle data.");
        }
        qid = 0;
#else
        error("SWIFT was not compiled with MPI support.");
#endif
        break;
      default:
        qid = -1;
    }

    if (qid >= s->nr_queues) error("Bad computed qid.");

    /* If no previous owner, pick a random queue. */
    if (qid < 0) qid = rand() % s->nr_queues;

    /* Increase the waiting counter. */
    atomic_inc(&s->waiting);

    /* Insert the task into that queue. */
    queue_insert(&s->queues[qid], t);
  }
}

/**
 * @brief Take care of a tasks dependencies.
 *
 * @param s The #scheduler.
 * @param t The finished #task.
 *
 * @return A pointer to the next task, if a suitable one has
 *         been identified.
 */
struct task *scheduler_done(struct scheduler *s, struct task *t) {

  /* Release whatever locks this task held. */
  if (!t->implicit) task_unlock(t);

  /* Loop through the dependencies and add them to a queue if
     they are ready. */
  for (int k = 0; k < t->nr_unlock_tasks; k++) {
    struct task *t2 = t->unlock_tasks[k];
    if (t2->skip) continue;

    const int res = atomic_dec(&t2->wait);
    if (res < 1) {
      error("Negative wait!");
    } else if (res == 1) {
      scheduler_enqueue(s, t2);
    }
  }

  /* Task definitely done, signal any sleeping runners. */
  if (!t->implicit) {
#ifdef SWIFT_DEBUG_TASKS
    t->toc = getticks();
#endif
    pthread_mutex_lock(&s->sleep_mutex);
    atomic_dec(&s->waiting);
    pthread_cond_broadcast(&s->sleep_cond);
    pthread_mutex_unlock(&s->sleep_mutex);
  }

  /* Mark the task as skip. */
  t->skip = 1;

  /* Return the next best task. Note that we currently do not
     implement anything that does this, as getting it to respect
     priorities is too tricky and currently unnecessary. */
  return NULL;
}

/**
 * @brief Resolve a single dependency by hand.
 *
 * @param s The #scheduler.
 * @param t The dependent #task.
 *
 * @return A pointer to the next task, if a suitable one has
 *         been identified.
 */
struct task *scheduler_unlock(struct scheduler *s, struct task *t) {

  /* Loop through the dependencies and add them to a queue if
     they are ready. */
  for (int k = 0; k < t->nr_unlock_tasks; k++) {
    struct task *t2 = t->unlock_tasks[k];
    const int res = atomic_dec(&t2->wait);
    if (res < 1) {
      error("Negative wait!");
    } else if (res == 1) {
      scheduler_enqueue(s, t2);
    }
  }

  /* Task definitely done. */
  if (!t->implicit) {
#ifdef SWIFT_DEBUG_TASKS
    t->toc = getticks();
#endif
    pthread_mutex_lock(&s->sleep_mutex);
    atomic_dec(&s->waiting);
    pthread_cond_broadcast(&s->sleep_cond);
    pthread_mutex_unlock(&s->sleep_mutex);
  }

  /* Return the next best task. Note that we currently do not
     implement anything that does this, as getting it to respect
     priorities is too tricky and currently unnecessary. */
  return NULL;
}

/**
 * @brief Get a task, preferably from the given queue.
 *
 * @param s The #scheduler.
 * @param qid The ID of the preferred #queue.
 * @param prev the previous task that was run.
 *
 * @return A pointer to a #task or @c NULL if there are no available tasks.
 */
struct task *scheduler_gettask(struct scheduler *s, int qid,
                               const struct task *prev) {

  struct task *res = NULL;
  const int nr_queues = s->nr_queues;
  unsigned int seed = qid;

  /* Check qid. */
  if (qid >= nr_queues || qid < 0) error("Bad queue ID.");

  /* Loop as long as there are tasks... */
  while (s->waiting > 0 && res == NULL) {

    /* Try more than once before sleeping. */
    for (int tries = 0; res == NULL && s->waiting && tries < scheduler_maxtries;
         tries++) {

      /* Try to get a task from the suggested queue. */
      if (s->queues[qid].count > 0 || s->queues[qid].count_incoming > 0) {
        TIMER_TIC
        res = queue_gettask(&s->queues[qid], prev, 0);
        TIMER_TOC(timer_qget);
        if (res != NULL) break;
      }

      /* If unsuccessful, try stealing from the other queues. */
      if (s->flags & scheduler_flag_steal) {
        int count = 0, qids[nr_queues];
        for (int k = 0; k < nr_queues; k++)
          if (s->queues[k].count > 0 || s->queues[k].count_incoming > 0) {
            qids[count++] = k;
          }
        for (int k = 0; k < scheduler_maxsteal && count > 0; k++) {
          const int ind = rand_r(&seed) % count;
          TIMER_TIC
          res = queue_gettask(&s->queues[qids[ind]], prev, 0);
          TIMER_TOC(timer_qsteal);
          if (res != NULL)
            break;
          else
            qids[ind] = qids[--count];
        }
        if (res != NULL) break;
      }
    }

/* If we failed, take a short nap. */
#ifdef WITH_MPI
    if (res == NULL && qid > 1)
#else
    if (res == NULL)
#endif
    {
      pthread_mutex_lock(&s->sleep_mutex);
      res = queue_gettask(&s->queues[qid], prev, 1);
      if (res == NULL && s->waiting > 0) {
        pthread_cond_wait(&s->sleep_cond, &s->sleep_mutex);
      }
      pthread_mutex_unlock(&s->sleep_mutex);
    }
  }

#ifdef SWIFT_DEBUG_TASKS
  /* Start the timer on this task, if we got one. */
  if (res != NULL) {
    res->tic = getticks();
    res->rid = qid;
  }
#endif

  /* No milk today. */
  return res;
}

/**
 * @brief Initialize the #scheduler.
 *
 * @param s The #scheduler.
 * @param space The #space we are working with
 * @param nr_tasks The number of tasks to allocate initially.
 * @param nr_queues The number of queues in this scheduler.
 * @param flags The #scheduler flags.
 * @param nodeID The MPI rank
 * @param tp Parallel processing threadpool.
 */
void scheduler_init(struct scheduler *s, struct space *space, int nr_tasks,
                    int nr_queues, unsigned int flags, int nodeID,
                    struct threadpool *tp) {

  /* Init the lock. */
  lock_init(&s->lock);

  /* Allocate the queues. */
  if (posix_memalign((void **)&s->queues, queue_struct_align,
                     sizeof(struct queue) * nr_queues) != 0)
    error("Failed to allocate queues.");

  /* Initialize each queue. */
  for (int k = 0; k < nr_queues; k++) queue_init(&s->queues[k], NULL);

  /* Init the sleep mutex and cond. */
  if (pthread_cond_init(&s->sleep_cond, NULL) != 0 ||
      pthread_mutex_init(&s->sleep_mutex, NULL) != 0)
    error("Failed to initialize sleep barrier.");

  /* Init the unlocks. */
  if ((s->unlocks = (struct task **)malloc(
           sizeof(struct task *) * scheduler_init_nr_unlocks)) == NULL ||
      (s->unlock_ind =
           (int *)malloc(sizeof(int) * scheduler_init_nr_unlocks)) == NULL)
    error("Failed to allocate unlocks.");
  s->nr_unlocks = 0;
  s->size_unlocks = scheduler_init_nr_unlocks;

  /* Set the scheduler variables. */
  s->nr_queues = nr_queues;
  s->flags = flags;
  s->space = space;
  s->nodeID = nodeID;
  s->threadpool = tp;

  /* Init the tasks array. */
  s->size = 0;
  s->tasks = NULL;
  s->tasks_ind = NULL;
  scheduler_reset(s, nr_tasks);
}

/**
 * @brief Prints the list of tasks to a file
 *
 * @param s The #scheduler
 * @param fileName Name of the file to write to
 */
void scheduler_print_tasks(const struct scheduler *s, const char *fileName) {

  const int nr_tasks = s->nr_tasks, *tid = s->tasks_ind;
  struct task *t, *tasks = s->tasks;

  FILE *file = fopen(fileName, "w");

  fprintf(file, "# Rank  Name  Subname  unlocks  waits\n");

  for (int k = nr_tasks - 1; k >= 0; k--) {
    t = &tasks[tid[k]];
    if (!((1 << t->type)) || t->skip) continue;
    fprintf(file, "%d %s %s %d %d\n", k, taskID_names[t->type],
            subtaskID_names[t->subtype], t->nr_unlock_tasks, t->wait);
  }

  fclose(file);
}

/**
 * @brief Frees up the memory allocated for this #scheduler
 */
void scheduler_clean(struct scheduler *s) {

  scheduler_free_tasks(s);
  free(s->unlocks);
  free(s->unlock_ind);
  for (int i = 0; i < s->nr_queues; ++i) queue_clean(&s->queues[i]);
  free(s->queues);
}

/**
 * @brief Free the task arrays allocated by this #scheduler.
 */
void scheduler_free_tasks(struct scheduler *s) {

  if (s->tasks != NULL) {
    free(s->tasks);
    s->tasks = NULL;
  }
  if (s->tasks_ind != NULL) {
    free(s->tasks_ind);
    s->tasks_ind = NULL;
  }
  if (s->tid_active != NULL) {
    free(s->tid_active);
    s->tid_active = NULL;
  }
  s->size = 0;
}
