/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <float.h>
#include <limits.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "task.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"
#include "lock.h"

/* Task type names. */
const char *taskID_names[task_type_count] = {
    "none",       "sort",      "self",          "pair",      "sub_self",
    "sub_pair",   "init",      "ghost",         "drift",     "kick",
    "kick_fixdt", "send",      "recv",          "grav_pp",   "grav_mm",
    "grav_up",    "grav_down", "grav_external", "part_sort", "gpart_sort",
    "split_cell", "rewait"};

const char *subtaskID_names[task_type_count] = {"none", "density", "force",
                                                "grav"};

/**
 * @brief Computes the overlap between the parts array of two given cells.
 */

size_t task_cell_overlap(const struct cell *ci, const struct cell *cj) {
  if (ci == NULL || cj == NULL) return 0;
  if (ci->parts <= cj->parts &&
      ci->parts + ci->count >= cj->parts + cj->count) {
    return cj->count;
  } else if (cj->parts <= ci->parts &&
             cj->parts + cj->count >= ci->parts + ci->count) {
    return ci->count;
  }
  return 0;
}

/**
 * @brief Compute the Jaccard similarity of the data used by two
 *        different tasks.
 *
 * @param ta The first #task.
 * @param tb The second #task.
 */

float task_overlap(const struct task *ta, const struct task *tb) {
  /* First check if any of the two tasks are of a type that don't
     use cells. */
  if (ta == NULL || tb == NULL || ta->type == task_type_none ||
      ta->type == task_type_part_sort || ta->type == task_type_gpart_sort ||
      ta->type == task_type_split_cell || ta->type == task_type_rewait ||
      tb->type == task_type_none || tb->type == task_type_part_sort ||
      tb->type == task_type_gpart_sort || tb->type == task_type_split_cell ||
      tb->type == task_type_rewait)
    return 0.0f;

  /* Compute the union of the cell data. */
  size_t size_union = 0;
  if (ta->ci != NULL) size_union += ta->ci->count;
  if (ta->cj != NULL) size_union += ta->cj->count;
  if (tb->ci != NULL) size_union += tb->ci->count;
  if (tb->cj != NULL) size_union += tb->cj->count;

  /* Compute the intersection of the cell data. */
  const size_t size_intersect =
      task_cell_overlap(ta->ci, tb->ci) + task_cell_overlap(ta->ci, tb->cj) +
      task_cell_overlap(ta->cj, tb->ci) + task_cell_overlap(ta->cj, tb->cj);

  return ((float)size_intersect) / (size_union - size_intersect);
}

/**
 * @brief Unlock the cell held by this task.
 *
 * @param t The #task.
 */

void task_unlock(struct task *t) {

  /* Act based on task type. */
  switch (t->type) {
    case task_type_self:
    case task_type_sub_self:
    case task_type_sort:
      cell_unlocktree(t->ci);
      break;
    case task_type_pair:
    case task_type_sub_pair:
      cell_unlocktree(t->ci);
      cell_unlocktree(t->cj);
      break;
    case task_type_grav_pp:
    case task_type_grav_mm:
    case task_type_grav_down:
      cell_gunlocktree(t->ci);
      if (t->cj != NULL) cell_gunlocktree(t->cj);
      break;
    default:
      break;
  }
}

/**
 * @brief Try to lock the cells associated with this task.
 *
 * @param t the #task.
 */

int task_lock(struct task *t) {

  const int type = t->type;
  const int subtype = t->subtype;
  struct cell *ci = t->ci, *cj = t->cj;
#ifdef WITH_MPI
  int res = 0, err = 0;
  MPI_Status stat;
#endif

  switch (type) {

    /* Communication task? */
    case task_type_recv:
    case task_type_send:
#ifdef WITH_MPI
      /* Check the status of the MPI request. */
      if ((err = MPI_Test(&t->req, &res, &stat)) != MPI_SUCCESS) {
        char buff[MPI_MAX_ERROR_STRING];
        int len;
        MPI_Error_string(err, buff, &len);
        error("Failed to test request on send/recv task (tag=%i, %s).",
              t->flags, buff);
      }
      return res;
#else
      error("SWIFT was not compiled with MPI support.");
#endif
      break;

    case task_type_sort:
      if (cell_locktree(ci) != 0) return 0;
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
        if (cell_glocktree(ci) != 0) return 0;
      } else {
        if (cell_locktree(ci) != 0) return 0;
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
        if (ci->ghold || cj->ghold) return 0;
        if (cell_glocktree(ci) != 0) return 0;
        if (cell_glocktree(cj) != 0) {
          cell_gunlocktree(ci);
          return 0;
        }
      } else {
        if (ci->hold || cj->hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
        if (cell_locktree(cj) != 0) {
          cell_unlocktree(ci);
          return 0;
        }
      }
      break;

    default:
      break;
  }

  /* If we made it this far, we've got a lock. */
  return 1;
}

/**
 * @brief Remove all unlocks to tasks that are of the given type.
 *
 * @param t The #task.
 * @param type The task type ID to remove.
 */

void task_cleanunlock(struct task *t, int type) {

  int k;

  lock_lock(&t->lock);

  for (k = 0; k < t->nr_unlock_tasks; k++)
    if (t->unlock_tasks[k]->type == type) {
      t->nr_unlock_tasks -= 1;
      t->unlock_tasks[k] = t->unlock_tasks[t->nr_unlock_tasks];
    }

  lock_unlock_blind(&t->lock);
}

/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 */

void task_rmunlock(struct task *ta, struct task *tb) {

  int k;

  lock_lock(&ta->lock);

  for (k = 0; k < ta->nr_unlock_tasks; k++)
    if (ta->unlock_tasks[k] == tb) {
      ta->nr_unlock_tasks -= 1;
      ta->unlock_tasks[k] = ta->unlock_tasks[ta->nr_unlock_tasks];
      lock_unlock_blind(&ta->lock);
      return;
    }
  error("Task not found.");
}

/**
 * @brief Remove an unlock_task from the given task.
 *
 * @param ta The unlocking #task.
 * @param tb The #task that will be unlocked.
 *
 * Differs from #task_rmunlock in that it will not fail if
 * the task @c tb is not in the unlocks of @c ta.
 */

void task_rmunlock_blind(struct task *ta, struct task *tb) {

  int k;

  lock_lock(&ta->lock);

  for (k = 0; k < ta->nr_unlock_tasks; k++)
    if (ta->unlock_tasks[k] == tb) {
      ta->nr_unlock_tasks -= 1;
      ta->unlock_tasks[k] = ta->unlock_tasks[ta->nr_unlock_tasks];
      break;
    }

  lock_unlock_blind(&ta->lock);
}

/**
 * @brief Prints the list of tasks contained in a given mask
 *
 * @param mask The mask to analyse
 */
void task_print_mask(unsigned int mask) {

  printf("task_print_mask: The tasks to run are [");
  for (int k = 1; k < task_type_count; k++)
    printf(" %s=%s", taskID_names[k], (mask & (1 << k)) ? "yes" : "no");
  printf(" ]\n");
}

/**
 * @brief Prints the list of subtasks contained in a given submask
 *
 * @param submask The submask to analyse
 */
void task_print_submask(unsigned int submask) {

  printf("task_print_submask: The subtasks to run are [");
  for (int k = 1; k < task_subtype_count; k++)
    printf(" %s=%s", subtaskID_names[k], (submask & (1 << k)) ? "yes" : "no");
  printf(" ]\n");
}
