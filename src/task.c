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
#include "inline.h"
#include "lock.h"

/* Task type names. */
const char *taskID_names[task_type_count] = {
    "none",       "sort",           "self",
    "pair",       "sub_self",       "sub_pair",
    "init_grav",  "ghost",          "extra_ghost",
    "drift_part", "drift_gpart",    "kick1",
    "kick2",      "timestep",       "send",
    "recv",       "grav_top_level", "grav_long_range",
    "grav_ghost", "grav_mm",        "grav_down",
    "cooling",    "sourceterms"};

/* Sub-task type names. */
const char *subtaskID_names[task_subtype_count] = {
    "none", "density", "gradient", "force", "grav",      "external_grav",
    "tend", "xv",      "rho",      "gpart", "multipole", "spart"};

/**
 * @brief Computes the overlap between the parts array of two given cells.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline)) INLINE static size_t task_cell_overlap_part(
    const struct cell *restrict ci, const struct cell *restrict cj) {

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
 * @brief Computes the overlap between the gparts array of two given cells.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline)) INLINE static size_t task_cell_overlap_gpart(
    const struct cell *restrict ci, const struct cell *restrict cj) {

  if (ci == NULL || cj == NULL) return 0;

  if (ci->gparts <= cj->gparts &&
      ci->gparts + ci->gcount >= cj->gparts + cj->gcount) {
    return cj->gcount;
  } else if (cj->gparts <= ci->gparts &&
             cj->gparts + cj->gcount >= ci->gparts + ci->gcount) {
    return ci->gcount;
  }

  return 0;
}

/**
 * @brief Returns the #task_actions for a given task.
 *
 * @param t The #task.
 */
__attribute__((always_inline)) INLINE static enum task_actions task_acts_on(
    const struct task *t) {

  switch (t->type) {

    case task_type_none:
      return task_action_none;
      break;

    case task_type_drift_part:
    case task_type_sort:
    case task_type_ghost:
    case task_type_extra_ghost:
    case task_type_cooling:
    case task_type_sourceterms:
      return task_action_part;
      break;

    case task_type_self:
    case task_type_pair:
    case task_type_sub_self:
    case task_type_sub_pair:
      switch (t->subtype) {

        case task_subtype_density:
        case task_subtype_gradient:
        case task_subtype_force:
          return task_action_part;
          break;

        case task_subtype_grav:
        case task_subtype_external_grav:
          return task_action_gpart;
          break;

        default:
          error("Unknow task_action for task");
          return task_action_none;
          break;
      }
      break;

    case task_type_kick1:
    case task_type_kick2:
    case task_type_timestep:
    case task_type_send:
    case task_type_recv:
      if (t->ci->count > 0 && t->ci->gcount > 0)
        return task_action_all;
      else if (t->ci->count > 0)
        return task_action_part;
      else if (t->ci->gcount > 0)
        return task_action_gpart;
      else
        error("Task without particles");
      break;

    case task_type_init_grav:
    case task_type_grav_top_level:
    case task_type_grav_long_range:
    case task_type_grav_mm:
      return task_action_multipole;
      break;

    case task_type_drift_gpart:
    case task_type_grav_down:
      return task_action_gpart;
      break;

    default:
      error("Unknown task_action for task");
      return task_action_none;
      break;
  }

  /* Silence compiler warnings */
  error("Unknown task_action for task");
  return task_action_none;
}

/**
 * @brief Compute the Jaccard similarity of the data used by two
 *        different tasks.
 *
 * @param ta The first #task.
 * @param tb The second #task.
 */
float task_overlap(const struct task *restrict ta,
                   const struct task *restrict tb) {

  if (ta == NULL || tb == NULL) return 0.f;

  const enum task_actions ta_act = task_acts_on(ta);
  const enum task_actions tb_act = task_acts_on(tb);

  /* First check if any of the two tasks are of a type that don't
     use cells. */
  if (ta_act == task_action_none || tb_act == task_action_none) return 0.f;

  const int ta_part = (ta_act == task_action_part || ta_act == task_action_all);
  const int ta_gpart =
      (ta_act == task_action_gpart || ta_act == task_action_all);
  const int tb_part = (tb_act == task_action_part || tb_act == task_action_all);
  const int tb_gpart =
      (tb_act == task_action_gpart || tb_act == task_action_all);

  /* In the case where both tasks act on parts */
  if (ta_part && tb_part) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->count;
    if (ta->cj != NULL) size_union += ta->cj->count;
    if (tb->ci != NULL) size_union += tb->ci->count;
    if (tb->cj != NULL) size_union += tb->cj->count;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_part(ta->ci, tb->ci) +
                                  task_cell_overlap_part(ta->ci, tb->cj) +
                                  task_cell_overlap_part(ta->cj, tb->ci) +
                                  task_cell_overlap_part(ta->cj, tb->cj);

    return ((float)size_intersect) / (size_union - size_intersect);
  }

  /* In the case where both tasks act on gparts */
  else if (ta_gpart && tb_gpart) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->gcount;
    if (ta->cj != NULL) size_union += ta->cj->gcount;
    if (tb->ci != NULL) size_union += tb->ci->gcount;
    if (tb->cj != NULL) size_union += tb->cj->gcount;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_gpart(ta->ci, tb->ci) +
                                  task_cell_overlap_gpart(ta->ci, tb->cj) +
                                  task_cell_overlap_gpart(ta->cj, tb->ci) +
                                  task_cell_overlap_gpart(ta->cj, tb->cj);

    return ((float)size_intersect) / (size_union - size_intersect);
  }

  /* Else, no overlap */
  return 0.f;
}

/**
 * @brief Unlock the cell held by this task.
 *
 * @param t The #task.
 */
void task_unlock(struct task *t) {

  const enum task_types type = t->type;
  const enum task_subtypes subtype = t->subtype;
  struct cell *ci = t->ci, *cj = t->cj;

  /* Act based on task type. */
  switch (type) {

    case task_type_kick1:
    case task_type_kick2:
    case task_type_timestep:
      cell_unlocktree(ci);
      cell_gunlocktree(ci);
      break;

    case task_type_drift_part:
    case task_type_sort:
      cell_unlocktree(ci);
      break;

    case task_type_drift_gpart:
      cell_gunlocktree(ci);
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
        cell_gunlocktree(ci);
        cell_munlocktree(ci);
      } else {
        cell_unlocktree(ci);
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
        cell_gunlocktree(ci);
        cell_gunlocktree(cj);
        cell_munlocktree(ci);
        cell_munlocktree(cj);
      } else {
        cell_unlocktree(ci);
        cell_unlocktree(cj);
      }
      break;

    case task_type_grav_down:
      cell_gunlocktree(ci);
      cell_munlocktree(ci);
      break;

    case task_type_grav_long_range:
    case task_type_grav_mm:
      cell_munlocktree(ci);
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

  const enum task_types type = t->type;
  const enum task_subtypes subtype = t->subtype;
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

    case task_type_kick1:
    case task_type_kick2:
    case task_type_timestep:
      if (ci->hold || ci->ghold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      if (cell_glocktree(ci) != 0) {
        cell_unlocktree(ci);
        return 0;
      }
      break;

    case task_type_drift_part:
    case task_type_sort:
      if (ci->hold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      break;

    case task_type_drift_gpart:
      if (ci->ghold) return 0;
      if (cell_glocktree(ci) != 0) return 0;
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
        /* Lock the gparts and the m-pole */
        if (ci->ghold || ci->mhold) return 0;
        if (cell_glocktree(ci) != 0)
          return 0;
        else if (cell_mlocktree(ci) != 0) {
          cell_gunlocktree(ci);
          return 0;
        }
      } else {
        if (cell_locktree(ci) != 0) return 0;
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
        /* Lock the gparts and the m-pole in both cells */
        if (ci->ghold || cj->ghold) return 0;
        if (cell_glocktree(ci) != 0) return 0;
        if (cell_glocktree(cj) != 0) {
          cell_gunlocktree(ci);
          return 0;
        } else if (cell_mlocktree(ci) != 0) {
          cell_gunlocktree(ci);
          cell_gunlocktree(cj);
          return 0;
        } else if (cell_mlocktree(cj) != 0) {
          cell_gunlocktree(ci);
          cell_gunlocktree(cj);
          cell_munlocktree(ci);
          return 0;
        }
      } else {
        /* Lock the parts in both cells */
        if (ci->hold || cj->hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
        if (cell_locktree(cj) != 0) {
          cell_unlocktree(ci);
          return 0;
        }
      }
      break;

    case task_type_grav_down:
      /* Lock the gparts and the m-poles */
      if (ci->ghold || ci->mhold) return 0;
      if (cell_glocktree(ci) != 0)
        return 0;
      else if (cell_mlocktree(ci) != 0) {
        cell_gunlocktree(ci);
        return 0;
      }
      break;

    case task_type_grav_long_range:
    case task_type_grav_mm:
      /* Lock the m-poles */
      if (ci->mhold) return 0;
      if (cell_mlocktree(ci) != 0) return 0;
      break;

    default:
      break;
  }

  /* If we made it this far, we've got a lock. */
  return 1;
}

/**
 * @brief Print basic information about a task.
 *
 * @param t The #task.
 */
void task_print(const struct task *t) {

  message("Type:'%s' sub_type:'%s' wait=%d nr_unlocks=%d skip=%d",
          taskID_names[t->type], subtaskID_names[t->subtype], t->wait,
          t->nr_unlock_tasks, t->skip);
}
