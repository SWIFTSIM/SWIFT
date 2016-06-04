/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "const.h"
#include "cycle.h"
#include "error.h"
#include "intrinsics.h"
#include "kernel_hydro.h"
#include "timers.h"

/**
 * @brief Add an unlock_task to the given task.
 *
 * @param s The #scheduler.
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */

void scheduler_addunlock(struct scheduler *s, struct task *ta,
                         struct task *tb) {

  /* Lock the scheduler since re-allocating the unlocks is not
     thread-safe. */
  if (lock_lock(&s->lock) != 0) error("Unable to lock scheduler.");

  /* Does the buffer need to be grown? */
  if (s->nr_unlocks == s->size_unlocks) {
    struct task **unlocks_new;
    int *unlock_ind_new;
    s->size_unlocks *= 2;
    if ((unlocks_new = (struct task **)malloc(sizeof(struct task *) *
                                              s->size_unlocks)) == NULL ||
        (unlock_ind_new = (int *)malloc(sizeof(int) * s->size_unlocks)) == NULL)
      error("Failed to re-allocate unlocks.");
    memcpy(unlocks_new, s->unlocks, sizeof(struct task *) * s->nr_unlocks);
    memcpy(unlock_ind_new, s->unlock_ind, sizeof(int) * s->nr_unlocks);
    free(s->unlocks);
    free(s->unlock_ind);
    s->unlocks = unlocks_new;
    s->unlock_ind = unlock_ind_new;
  }

  /* Write the unlock to the scheduler. */
  const int ind = atomic_inc(&s->nr_unlocks);
  s->unlocks[ind] = tb;
  s->unlock_ind[ind] = ta - s->tasks;

  /* Release the scheduler. */
  if (lock_unlock(&s->lock) != 0) error("Unable to unlock scheduler.");
}

/**
 * @brief Split tasks that may be too large.
 *
 * @param s The #scheduler we are working in.
 */

void scheduler_splittasks(struct scheduler *s) {

  const int pts[7][8] = {
      {-1, 12, 10, 9, 4, 3, 1, 0},     {-1, -1, 11, 10, 5, 4, 2, 1},
      {-1, -1, -1, 12, 7, 6, 4, 3},    {-1, -1, -1, -1, 8, 7, 5, 4},
      {-1, -1, -1, -1, -1, 12, 10, 9}, {-1, -1, -1, -1, -1, -1, 11, 10},
      {-1, -1, -1, -1, -1, -1, -1, 12}};
  const float sid_scale[13] = {0.1897, 0.4025, 0.1897, 0.4025, 0.5788,
                               0.4025, 0.1897, 0.4025, 0.1897, 0.4025,
                               0.5788, 0.4025, 0.5788};

  /* Loop through the tasks... */
  int tid = 0, redo = 0;
  struct task *t_old = NULL;
  while (1) {

    /* Get a pointer on the task. */
    struct task *t = t_old;
    if (redo) {
      redo = 0;
    } else {
      const int ind = atomic_inc(&tid);
      if (ind < s->nr_tasks)
        t_old = t = &s->tasks[s->tasks_ind[ind]];
      else
        break;
    }

    /* Skip sorting tasks. */
    if (t->type == task_type_part_sort) continue;

    if (t->type == task_type_gpart_sort) continue;

    /* Empty task? */
    if (t->ci == NULL || (t->type == task_type_pair && t->cj == NULL)) {
      t->type = task_type_none;
      t->skip = 1;
      continue;
    }

    /* Non-local kick task? */
    if ((t->type == task_type_kick) && t->ci->nodeID != s->nodeID) {
      t->type = task_type_none;
      t->skip = 1;
      continue;
    }

    /* Non-local drift task? */
    if ((t->type == task_type_drift) && t->ci->nodeID != s->nodeID) {
      t->type = task_type_none;
      t->skip = 1;
      continue;
    }

    /* Non-local init task? */
    if ((t->type == task_type_init) && t->ci->nodeID != s->nodeID) {
      t->type = task_type_none;
      t->skip = 1;
      continue;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {

      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        continue;
      }

      /* Is this cell even split? */
      if (ci->split) {

        /* Make a sub? */
        if (scheduler_dosub && ci->count < space_subsize / ci->count) {

          /* convert to a self-subtask. */
          t->type = task_type_sub;

        }

        /* Otherwise, make tasks explicitly. */
        else {

          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self task. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;
          t->ci = ci->progeny[first_child];
          for (int k = first_child + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                ci->progeny[k], NULL, 0);

          /* Make a task for each pair of progeny. */
          for (int j = 0; j < 8; j++)
            if (ci->progeny[j] != NULL)
              for (int k = j + 1; k < 8; k++)
                if (ci->progeny[k] != NULL)
                  scheduler_addtask(s, task_type_pair, t->subtype, pts[j][k], 0,
                                    ci->progeny[j], ci->progeny[k], 0);
        }
      }

    }

    /* Pair interaction? */
    else if (t->type == task_type_pair) {

      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;
      const double hi = ci->dmin;
      const double hj = cj->dmin;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        continue;
      }

      /* Get the sort ID, use space_getsid and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      int sid = space_getsid(s->space, &ci, &cj, shift);

      /* Should this task be split-up? */
      if (ci->split && cj->split &&
          ci->h_max * kernel_gamma * space_stretch < hi / 2 &&
          cj->h_max * kernel_gamma * space_stretch < hj / 2) {

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            ci->count * sid_scale[sid] < space_subsize / cj->count &&
            sid != 0 && sid != 2 && sid != 6 && sid != 8) {

          /* Make this task a sub task. */
          t->type = task_type_sub;

        }

        /* Otherwise, split it. */
        else {

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
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[7], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[6], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              break;

            case 2: /* (  1 ,  1 , -1 ) */
              t->ci = ci->progeny[6];
              t->cj = cj->progeny[1];
              t->flags = 2;
              t->tight = 1;
              break;

            case 3: /* (  1 ,  0 ,  1 ) */
              t->ci = ci->progeny[5];
              t->cj = cj->progeny[0];
              t->flags = 3;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[7], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[5], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              break;

            case 4: /* (  1 ,  0 ,  0 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[0];
              t->flags = 4;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[5], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[6], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[4], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[5], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[6], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[7], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[4], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[5], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[6], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[7], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[4], cj->progeny[3], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[5], cj->progeny[3], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[6], cj->progeny[3], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 4, 0,
                                    ci->progeny[7], cj->progeny[3], 1);
              break;

            case 5: /* (  1 ,  0 , -1 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[1];
              t->flags = 5;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[6], cj->progeny[3], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[4], cj->progeny[3], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[6], cj->progeny[1], 1);
              break;

            case 6: /* (  1 , -1 ,  1 ) */
              t->ci = ci->progeny[5];
              t->cj = cj->progeny[2];
              t->flags = 6;
              t->tight = 1;
              break;

            case 7: /* (  1 , -1 ,  0 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[3];
              t->flags = 6;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[5], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[4], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[5], cj->progeny[3], 1);
              break;

            case 8: /* (  1 , -1 , -1 ) */
              t->ci = ci->progeny[4];
              t->cj = cj->progeny[3];
              t->flags = 8;
              t->tight = 1;
              break;

            case 9: /* (  0 ,  1 ,  1 ) */
              t->ci = ci->progeny[3];
              t->cj = cj->progeny[0];
              t->flags = 9;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[7], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[3], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              break;

            case 10: /* (  0 ,  1 ,  0 ) */
              t->ci = ci->progeny[2];
              t->cj = cj->progeny[0];
              t->flags = 10;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[3], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[6], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[2], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[3], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[6], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 7, 0,
                                    ci->progeny[7], cj->progeny[1], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[2], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[3], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[6], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[7], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[2], cj->progeny[5], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 1, 0,
                                    ci->progeny[3], cj->progeny[5], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[6], cj->progeny[5], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 10, 0,
                                    ci->progeny[7], cj->progeny[5], 1);
              break;

            case 11: /* (  0 ,  1 , -1 ) */
              t->ci = ci->progeny[2];
              t->cj = cj->progeny[1];
              t->flags = 11;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[6], cj->progeny[5], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[2], cj->progeny[5], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[6], cj->progeny[1], 1);
              break;

            case 12: /* (  0 ,  0 ,  1 ) */
              t->ci = ci->progeny[1];
              t->cj = cj->progeny[0];
              t->flags = 12;
              t->tight = 1;
              t = scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[3], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[5], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 2, 0,
                                    ci->progeny[7], cj->progeny[0], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[1], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[3], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 8, 0,
                                    ci->progeny[5], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 5, 0,
                                    ci->progeny[7], cj->progeny[2], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[1], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 6, 0,
                                    ci->progeny[3], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[5], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 11, 0,
                                    ci->progeny[7], cj->progeny[4], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                    ci->progeny[1], cj->progeny[6], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 3, 0,
                                    ci->progeny[3], cj->progeny[6], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 9, 0,
                                    ci->progeny[5], cj->progeny[6], 1);
              t = scheduler_addtask(s, task_type_pair, t->subtype, 12, 0,
                                    ci->progeny[7], cj->progeny[6], 1);
              break;
          }
        }

      } /* split this task? */

      /* Otherwise, break it up if it is too large? */
      else if (scheduler_doforcesplit && ci->split && cj->split &&
               (ci->count > space_maxsize / cj->count)) {

        // message( "force splitting pair with %i and %i parts." , ci->count ,
        // cj->count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL) {
                t = scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k], 0);
                t->flags = space_getsid(s->space, &t->ci, &t->cj, shift);
              }

      }

      /* Otherwise, if not spilt, stitch-up the sorting. */
      else {

        /* Create the sort for ci. */
        // lock_lock( &ci->lock );
        if (ci->sorts == NULL)
          ci->sorts =
              scheduler_addtask(s, task_type_sort, 0, 1 << sid, 0, ci, NULL, 0);
        else
          ci->sorts->flags |= (1 << sid);
        // lock_unlock_blind( &ci->lock );
        scheduler_addunlock(s, ci->sorts, t);

        /* Create the sort for cj. */
        // lock_lock( &cj->lock );
        if (cj->sorts == NULL)
          cj->sorts =
              scheduler_addtask(s, task_type_sort, 0, 1 << sid, 0, cj, NULL, 0);
        else
          cj->sorts->flags |= (1 << sid);
        // lock_unlock_blind( &cj->lock );
        scheduler_addunlock(s, cj->sorts, t);
      }

    } /* pair interaction? */

    /* Gravity interaction? */
    else if (t->type == task_type_grav_mm) {

      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Self-interaction? */
      if (cj == NULL) {

        /* Ignore this task if the cell has no gparts. */
        if (ci->gcount == 0) t->type = task_type_none;

        /* If the cell is split, recurse. */
        else if (ci->split) {

          /* Make a single sub-task? */
          if (scheduler_dosub && ci->gcount < space_subsize / ci->gcount) {

            t->type = task_type_sub;
            t->subtype = task_subtype_grav;

          }

          /* Otherwise, just split the task. */
          else {

            /* Split this task into tasks on its progeny. */
            t->type = task_type_none;
            for (int j = 0; j < 8; j++)
              if (ci->progeny[j] != NULL && ci->progeny[j]->gcount > 0) {
                if (t->type == task_type_none) {
                  t->type = task_type_grav_mm;
                  t->ci = ci->progeny[j];
                  t->cj = NULL;
                } else
                  t = scheduler_addtask(s, task_type_grav_mm, task_subtype_none,
                                        0, 0, ci->progeny[j], NULL, 0);
                for (int k = j + 1; k < 8; k++)
                  if (ci->progeny[k] != NULL && ci->progeny[k]->gcount > 0) {
                    if (t->type == task_type_none) {
                      t->type = task_type_grav_mm;
                      t->ci = ci->progeny[j];
                      t->cj = ci->progeny[k];
                    } else
                      t = scheduler_addtask(s, task_type_grav_mm,
                                            task_subtype_none, 0, 0,
                                            ci->progeny[j], ci->progeny[k], 0);
                  }
              }
            redo = (t->type != task_type_none);
          }

        }

        /* Otherwise, just make a pp task out of it. */
        else
          t->type = task_type_grav_pp;

      }

      /* Nope, pair. */
      else {

        /* Make a sub-task? */
        if (scheduler_dosub && ci->gcount < space_subsize / cj->gcount) {

          t->type = task_type_sub;
          t->subtype = task_subtype_grav;

        }

        /* Otherwise, split the task. */
        else {

          /* Get the opening angle theta. */
          float dx[3], theta;
          for (int k = 0; k < 3; k++) {
            dx[k] = fabs(ci->loc[k] - cj->loc[k]);
            if (s->space->periodic && dx[k] > 0.5 * s->space->dim[k])
              dx[k] = -dx[k] + s->space->dim[k];
            if (dx[k] > 0.0f) dx[k] -= ci->h[k];
          }
          theta =
              (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]) /
              (ci->h[0] * ci->h[0] + ci->h[1] * ci->h[1] + ci->h[2] * ci->h[2]);

          /* Ignore this task if the cell has no gparts. */
          if (ci->gcount == 0 || cj->gcount == 0) t->type = task_type_none;

          /* Split the interaction? */
          else if (theta < const_theta_max * const_theta_max) {

            /* Are both ci and cj split? */
            if (ci->split && cj->split) {

              /* Split this task into tasks on its progeny. */
              t->type = task_type_none;
              for (int j = 0; j < 8; j++)
                if (ci->progeny[j] != NULL && ci->progeny[j]->gcount > 0) {
                  for (int k = 0; k < 8; k++)
                    if (cj->progeny[k] != NULL && cj->progeny[k]->gcount > 0) {
                      if (t->type == task_type_none) {
                        t->type = task_type_grav_mm;
                        t->ci = ci->progeny[j];
                        t->cj = cj->progeny[k];
                      } else
                        t = scheduler_addtask(
                            s, task_type_grav_mm, task_subtype_none, 0, 0,
                            ci->progeny[j], cj->progeny[k], 0);
                    }
                }
              redo = (t->type != task_type_none);

            }

            /* Otherwise, make a pp task out of it. */
            else
              t->type = task_type_grav_pp;
          }
        }

      } /* gravity pair interaction? */

    } /* gravity interaction? */

  } /* loop over all tasks. */
}

/**
 * @brief Add a #task to the #scheduler.
 *
 * @param s The #scheduler we are working in.
 * @param type The type of the task.
 * @param subtype The sub-type of the task.
 * @param flags The flags of the task.
 * @param wait
 * @param ci The first cell to interact.
 * @param cj The second cell to interact.
 * @param tight
 */

struct task *scheduler_addtask(struct scheduler *s, int type, int subtype,
                               int flags, int wait, struct cell *ci,
                               struct cell *cj, int tight) {

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
  t->wait = wait;
  t->ci = ci;
  t->cj = cj;
  t->skip = 0;
  t->tight = tight;
  t->implicit = 0;
  t->weight = 0;
  t->rank = 0;
  t->tic = 0;
  t->toc = 0;
  t->nr_unlock_tasks = 0;
  t->rid = -1;
  t->last_rid = -1;

  /* Init the lock. */
  lock_init(&t->lock);

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
  int *counts;
  if ((counts = (int *)malloc(sizeof(int) * s->nr_tasks)) == NULL)
    error("Failed to allocate temporary counts array.");
  bzero(counts, sizeof(int) * s->nr_tasks);
  for (int k = 0; k < s->nr_unlocks; k++) counts[s->unlock_ind[k]] += 1;

  /* Compute the offset for each unlock block. */
  int *offsets;
  if ((offsets = (int *)malloc(sizeof(int) * (s->nr_tasks + 1))) == NULL)
    error("Failed to allocate temporary offsets array.");
  offsets[0] = 0;
  for (int k = 0; k < s->nr_tasks; k++) offsets[k + 1] = offsets[k] + counts[k];

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
  for (int k = 0; k < s->nr_tasks; k++)
    for (int j = offsets[k]; j < offsets[k + 1]; j++) s->unlock_ind[j] = k;

  /* Set the unlocks in the tasks. */
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &s->tasks[k];
    t->nr_unlock_tasks = counts[k];
    t->unlock_tasks = &s->unlocks[offsets[k]];
    for (int j = offsets[k]; j < offsets[k + 1]; j++) s->unlock_ind[j] = k;
  }

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
  for (int k = 0; k < nr_tasks; k++) {
    tid[k] = k;
    for (int j = 0; j < tasks[k].nr_unlock_tasks; j++)
      tasks[k].unlock_tasks[j]->wait += 1;
  }

  /* Main loop. */
  for (int j = 0, rank = 0, left = 0; left < nr_tasks; rank++) {

    /* Load the tids of tasks with no waits. */
    for (int k = left; k < nr_tasks; k++)
      if (tasks[tid[k]].wait == 0) {
        int temp = tid[j];
        tid[j] = tid[k];
        tid[k] = temp;
        j += 1;
      }

    /* Did we get anything? */
    if (j == left) error("Unsatisfiable task dependencies detected.");

    /* Unlock the next layer of tasks. */
    for (int i = left; i < j; i++) {
      struct task *t = &tasks[tid[i]];
      t->rank = rank;
      tid[i] = t - tasks;
      if (tid[i] >= nr_tasks) error("Task index overshoot.");
      /* message( "task %i of type %s has rank %i." , i ,
          (t->type == task_type_self) ? "self" : (t->type == task_type_pair) ?
         "pair" : "sort" , rank ); */
      for (int k = 0; k < t->nr_unlock_tasks; k++)
        t->unlock_tasks[k]->wait -= 1;
    }

    /* The new left (no, not tony). */
    left = j;
  }

  /* Verify that the tasks were ranked correctly. */
  /* for ( k = 1 ; k < s->nr_tasks ; k++ )
      if ( tasks[ tid[k-1] ].rank > tasks[ tid[k-1] ].rank )
          error( "Task ranking failed." ); */
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
    if (s->tasks != NULL) free(s->tasks);
    if (s->tasks_ind != NULL) free(s->tasks_ind);

    /* Allocate the new lists. */
    if ((s->tasks = (struct task *)malloc(sizeof(struct task) * size)) ==
            NULL ||
        (s->tasks_ind = (int *)malloc(sizeof(int) * size)) == NULL)
      error("Failed to allocate task lists.");
  }

  /* Reset the task data. */
  bzero(s->tasks, sizeof(struct task) * size);

  /* Reset the counters. */
  s->size = size;
  s->nr_tasks = 0;
  s->tasks_next = 0;
  s->waiting = 0;
  s->mask = 0;
  s->submask = 0;
  s->nr_unlocks = 0;

  /* Set the task pointers in the queues. */
  for (int k = 0; k < s->nr_queues; k++) s->queues[k].tasks = s->tasks;
}

/**
 * @brief Compute the task weights
 *
 * @param s The #scheduler.
 */

void scheduler_reweight(struct scheduler *s) {

  const int nr_tasks = s->nr_tasks;
  int *tid = s->tasks_ind;
  struct task *tasks = s->tasks;
  const int nodeID = s->nodeID;
  const float sid_scale[13] = {0.1897, 0.4025, 0.1897, 0.4025, 0.5788,
                               0.4025, 0.1897, 0.4025, 0.1897, 0.4025,
                               0.5788, 0.4025, 0.5788};
  const float wscale = 0.001;
  // ticks tic;

  /* Run through the tasks backwards and set their waits and
     weights. */
  // tic = getticks();
  for (int k = nr_tasks - 1; k >= 0; k--) {
    struct task *t = &tasks[tid[k]];
    t->weight = 0;
    for (int j = 0; j < t->nr_unlock_tasks; j++)
      if (t->unlock_tasks[j]->weight > t->weight)
        t->weight = t->unlock_tasks[j]->weight;
    if (!t->implicit && t->tic > 0)
      t->weight += wscale * (t->toc - t->tic);
    else
      switch (t->type) {
        case task_type_sort:
          t->weight += wscale * intrinsics_popcount(t->flags) * t->ci->count *
                       (sizeof(int) * 8 - intrinsics_clz(t->ci->count));
          break;
        case task_type_self:
          t->weight += 1 * t->ci->count * t->ci->count;
          break;
        case task_type_pair:
          if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID)
            t->weight +=
                3 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
          else
            t->weight +=
                2 * wscale * t->ci->count * t->cj->count * sid_scale[t->flags];
          break;
        case task_type_sub:
          if (t->cj != NULL) {
            if (t->ci->nodeID != nodeID || t->cj->nodeID != nodeID) {
              if (t->flags < 0)
                t->weight += 3 * wscale * t->ci->count * t->cj->count;
              else
                t->weight += 3 * wscale * t->ci->count * t->cj->count *
                             sid_scale[t->flags];
            } else {
              if (t->flags < 0)
                t->weight += 2 * wscale * t->ci->count * t->cj->count;
              else
                t->weight += 2 * wscale * t->ci->count * t->cj->count *
                             sid_scale[t->flags];
            }
          } else
            t->weight += 1 * wscale * t->ci->count * t->ci->count;
          break;
        case task_type_ghost:
          if (t->ci == t->ci->super) t->weight += wscale * t->ci->count;
          break;
        case task_type_kick:
          t->weight += wscale * t->ci->count;
          break;
        case task_type_drift:
          t->weight += wscale * t->ci->count;
          break;
        case task_type_init:
          t->weight += wscale * t->ci->count;
          break;
        default:
          break;
      }
    if (t->type == task_type_send) t->weight = INT_MAX / 8;
    if (t->type == task_type_recv) t->weight *= 1.41;
  }
  // message( "weighting tasks took %.3f %s." ,
  // clocks_from_ticks( getticks() - tic ), clocks_getunit());

  /* int min = tasks[0].weight, max = tasks[0].weight;
  for ( k = 1 ; k < nr_tasks ; k++ )
      if ( tasks[k].weight < min )
          min = tasks[k].weight;
      else if ( tasks[k].weight > max )
          max = tasks[k].weight;
  message( "task weights are in [ %i , %i ]." , min , max ); */
}

/**
 * @brief Start the scheduler, i.e. fill the queues with ready tasks.
 *
 * @param s The #scheduler.
 * @param mask The task types to enqueue.
 * @param submask The sub-task types to enqueue.
 */

void scheduler_start(struct scheduler *s, unsigned int mask,
                     unsigned int submask) {

  const int nr_tasks = s->nr_tasks;
  int *tid = s->tasks_ind;
  struct task *tasks = s->tasks;
  // ticks tic;

  /* Store the masks */
  s->mask = mask | (1 << task_type_rewait);
  s->submask = submask | (1 << task_subtype_none);

  /* Clear all the waits and rids. */
  // ticks tic = getticks();
  for (int k = 0; k < s->nr_tasks; k++) {
    s->tasks[k].wait = 1;
    s->tasks[k].rid = -1;
  }
  // message( "waiting tasks took %.3f %s." ,
  // clocks_from_ticks(getticks() - tic), clocks_getunit() );

  /* Enqueue a set of extraenous tasks to set the task waits. */
  struct task *rewait_tasks = &s->tasks[s->nr_tasks];
  const int num_rewait_tasks = s->nr_queues > s->size - s->nr_tasks
                                   ? s->size - s->nr_tasks
                                   : s->nr_queues;

  /* Remember that engine_launch may fiddle with this value. */
  const int waiting_old = s->waiting;

  /* We are going to use the task structure in a modified way to pass
     information to the task. Don't do this at home !
     - ci and cj will give the range of tasks to which the waits will be applied
     - the flags will be used to transfer the mask
     - the rank will be used to transfer the submask
     - the rest is unused.
  */
  for (int k = 0; k < num_rewait_tasks; k++) {
    rewait_tasks[k].type = task_type_rewait;
    rewait_tasks[k].ci = (struct cell *)&s->tasks[k * nr_tasks / s->nr_queues];
    rewait_tasks[k].cj =
        (struct cell *)&s->tasks[(k + 1) * nr_tasks / s->nr_queues];
    rewait_tasks[k].flags = s->mask;
    rewait_tasks[k].rank = s->submask;
    rewait_tasks[k].skip = 0;
    rewait_tasks[k].wait = 0;
    rewait_tasks[k].rid = -1;
    rewait_tasks[k].weight = 1;
    rewait_tasks[k].implicit = 0;
    rewait_tasks[k].nr_unlock_tasks = 0;
    scheduler_enqueue(s, &rewait_tasks[k]);
    pthread_cond_broadcast(&s->sleep_cond);
  }

  /* Wait for the rewait tasks to have executed. */
  pthread_mutex_lock(&s->sleep_mutex);
  pthread_cond_broadcast(&s->sleep_cond);
  while (s->waiting > waiting_old) {
    pthread_cond_wait(&s->sleep_cond, &s->sleep_mutex);
  }
  pthread_mutex_unlock(&s->sleep_mutex);
  /* message("waiting tasks took %.3f %s.",
     clocks_from_ticks(getticks() - tic), clocks_getunit());*/

  s->mask = mask;
  s->submask = submask | (1 << task_subtype_none);

  /* Loop over the tasks and enqueue whoever is ready. */
  // tic = getticks();
  for (int k = 0; k < s->nr_tasks; k++) {
    struct task *t = &tasks[tid[k]];
    if (atomic_dec(&t->wait) == 1 && ((1 << t->type) & s->mask) &&
        ((1 << t->subtype) & s->submask) && !t->skip) {
      scheduler_enqueue(s, t);
      pthread_cond_broadcast(&s->sleep_cond);
    }
  }

  /* To be safe, fire of one last sleep_cond in a safe way. */
  pthread_mutex_lock(&s->sleep_mutex);
  pthread_cond_broadcast(&s->sleep_cond);
  pthread_mutex_unlock(&s->sleep_mutex);

  // message( "enqueueing tasks took %.3f %s." ,
  // clocks_from_ticks( getticks() - tic ), clocks_getunit());
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

  /* Fail if this task has already been enqueued before. */
  if (t->rid >= 0) error("Task has already been enqueued.");

  /* Ignore skipped tasks and tasks not in the masks. */
  if (t->skip || (1 << t->type) & ~(s->mask) ||
      (1 << t->subtype) & ~(s->submask)) {
    return;
  }

  /* If this is an implicit task, just pretend it's done. */
  if (t->implicit) {
    for (int j = 0; j < t->nr_unlock_tasks; j++) {
      struct task *t2 = t->unlock_tasks[j];

      if (atomic_dec(&t2->wait) == 1) scheduler_enqueue(s, t2);
    }
  }

  /* Otherwise, look for a suitable queue. */
  else {
#ifdef WITH_MPI
    int err;
#endif

    /* Find the previous owner for each task type, and do
       any pre-processing needed. */
    switch (t->type) {
      case task_type_self:
      case task_type_sort:
      case task_type_ghost:
      case task_type_kick:
      case task_type_drift:
      case task_type_init:
        qid = t->ci->super->owner;
        break;
      case task_type_pair:
      case task_type_sub:
        qid = t->ci->super->owner;
        if (t->cj != NULL &&
            (qid < 0 ||
             s->queues[qid].count > s->queues[t->cj->super->owner].count))
          qid = t->cj->super->owner;
        break;
      case task_type_recv:
#ifdef WITH_MPI
        err = MPI_Irecv(t->ci->parts, t->ci->count, part_mpi_type,
                        t->ci->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit irecv for particle data.");
        }
        // message( "receiving %i parts with tag=%i from %i to %i." ,
        //     t->ci->count , t->flags , t->ci->nodeID , s->nodeID );
        // fflush(stdout);
        qid = 1 % s->nr_queues;
#else
        error("SWIFT was not compiled with MPI support.");
#endif
        break;
      case task_type_send:
#ifdef WITH_MPI
        err = MPI_Isend(t->ci->parts, t->ci->count, part_mpi_type,
                        t->cj->nodeID, t->flags, MPI_COMM_WORLD, &t->req);
        if (err != MPI_SUCCESS) {
          mpi_error(err, "Failed to emit isend for particle data.");
        }
        // message( "sending %i parts with tag=%i from %i to %i." ,
        //     t->ci->count , t->flags , s->nodeID , t->cj->nodeID );
        // fflush(stdout);
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

    const int res = atomic_dec(&t2->wait);
    if (res < 1) {
      error("Negative wait!");
    } else if (res == 1) {
      scheduler_enqueue(s, t2);
    }
  }

  /* Task definitely done, signal any sleeping runners. */
  if (!t->implicit) {
    t->toc = getticks();
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
    t->toc = getticks();
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
      if (s->queues[qid].count > 0) {
        TIMER_TIC
        res = queue_gettask(&s->queues[qid], prev, 0);
        TIMER_TOC(timer_qget);
        if (res != NULL) break;
      }

      /* If unsuccessful, try stealing from the other queues. */
      if (s->flags & scheduler_flag_steal) {
        int count = 0, qids[nr_queues];
        for (int k = 0; k < nr_queues; k++)
          if (s->queues[k].count > 0) qids[count++] = k;
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
    if (res == NULL && qid > 1) {
#else
    if (res == NULL) {
#endif
      pthread_mutex_lock(&s->sleep_mutex);
      res = queue_gettask(&s->queues[qid], prev, 1);
      if (res == NULL && s->waiting > 0) {
        pthread_cond_wait(&s->sleep_cond, &s->sleep_mutex);
      }
      pthread_mutex_unlock(&s->sleep_mutex);
    }
  }

  /* Start the timer on this task, if we got one. */
  if (res != NULL) {
    res->tic = getticks();
    res->rid = qid;
  }

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
 */

void scheduler_init(struct scheduler *s, struct space *space, int nr_tasks,
                    int nr_queues, unsigned int flags, int nodeID) {

  /* Init the lock. */
  lock_init(&s->lock);

  /* Allocate the queues. */
  if ((s->queues = (struct queue *)malloc(sizeof(struct queue) * nr_queues)) ==
      NULL)
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
 * @brief Sets the waits of the dependants of a range of task
 *
 * @param t_begin Beginning of the #task range
 * @param t_end End of the #task range
 * @param mask The scheduler task mask
 * @param submask The scheduler subtask mask
 */
void scheduler_do_rewait(struct task *t_begin, struct task *t_end,
                         unsigned int mask, unsigned int submask) {
  for (struct task *t2 = t_begin; t2 != t_end; t2++) {

    if (t2->skip) continue;

    /* Skip tasks not in the mask */
    if (!((1 << t2->type) & mask) || !((1 << t2->subtype) & submask)) continue;

    /* Skip sort tasks that have already been performed */
    if (t2->type == task_type_sort && t2->flags == 0) continue;

    /* Sets the waits of the dependances */
    for (int k = 0; k < t2->nr_unlock_tasks; k++) {
      struct task *t3 = t2->unlock_tasks[k];
      atomic_inc(&t3->wait);
    }
  }
}
