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
#include "engine.h"
#include "error.h"
#include "inline.h"
#include "lock.h"

/* Task type names. */
const char *taskID_names[task_type_count] = {"none",
                                             "sort",
                                             "self",
                                             "pair",
                                             "sub_self",
                                             "sub_pair",
                                             "init_grav",
                                             "init_grav_out",
                                             "ghost_in",
                                             "ghost",
                                             "ghost_out",
                                             "extra_ghost",
                                             "drift_part",
                                             "drift_spart",
                                             "drift_gpart",
                                             "drift_gpart_out",
                                             "end_hydro_force",
                                             "kick1",
                                             "kick2",
                                             "timestep",
                                             "timestep_limiter",
                                             "send",
                                             "recv",
                                             "grav_long_range",
                                             "grav_mm",
                                             "grav_down_in",
                                             "grav_down",
                                             "grav_mesh",
                                             "grav_end_force",
                                             "cooling",
                                             "star_formation",
                                             "logger",
                                             "stars_in",
                                             "stars_out",
                                             "stars_ghost_in",
                                             "stars_ghost",
                                             "stars_ghost_out",
                                             "stars_sort",
                                             "fof_self",
                                             "fof_pair"};

/* Sub-task type names. */
const char *subtaskID_names[task_subtype_count] = {
    "none",    "density",       "gradient",      "force",
    "limiter", "grav",          "external_grav", "tend",
    "xv",      "rho",           "gpart",         "multipole",
    "spart",   "stars_density", "stars_feedback"};

#ifdef WITH_MPI
/* MPI communicators for the subtypes. */
MPI_Comm subtaskMPI_comms[task_subtype_count];
#endif

/**
 * @brief Computes the overlap between the parts array of two given cells.
 *
 * @param TYPE is the type of parts (e.g. #part, #gpart, #spart)
 * @param ARRAY is the array of this specific type.
 * @param COUNT is the number of elements in the array.
 */
#define TASK_CELL_OVERLAP(TYPE, ARRAY, COUNT)                               \
  __attribute__((always_inline))                                            \
      INLINE static size_t task_cell_overlap_##TYPE(                        \
          const struct cell *restrict ci, const struct cell *restrict cj) { \
                                                                            \
    if (ci == NULL || cj == NULL) return 0;                                 \
                                                                            \
    if (ci->ARRAY <= cj->ARRAY &&                                           \
        ci->ARRAY + ci->COUNT >= cj->ARRAY + cj->COUNT) {                   \
      return cj->COUNT;                                                     \
    } else if (cj->ARRAY <= ci->ARRAY &&                                    \
               cj->ARRAY + cj->COUNT >= ci->ARRAY + ci->COUNT) {            \
      return ci->COUNT;                                                     \
    }                                                                       \
                                                                            \
    return 0;                                                               \
  }

TASK_CELL_OVERLAP(part, hydro.parts, hydro.count);
TASK_CELL_OVERLAP(gpart, grav.parts, grav.count);
TASK_CELL_OVERLAP(spart, stars.parts, stars.count);

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
    case task_type_timestep_limiter:
    case task_type_cooling:
    case task_type_end_hydro_force:
      return task_action_part;
      break;

    case task_type_star_formation:
      return task_action_all;

    case task_type_drift_spart:
    case task_type_stars_ghost:
    case task_type_stars_sort:
      return task_action_spart;
      break;

    case task_type_self:
    case task_type_pair:
    case task_type_sub_self:
    case task_type_sub_pair:
      switch (t->subtype) {

        case task_subtype_density:
        case task_subtype_gradient:
        case task_subtype_force:
        case task_subtype_limiter:
          return task_action_part;
          break;

        case task_subtype_stars_density:
        case task_subtype_stars_feedback:
          return task_action_all;
          break;

        case task_subtype_grav:
        case task_subtype_external_grav:
          return task_action_gpart;
          break;

        default:
#ifdef SWIFT_DEBUG_CHECKS
          error("Unknown task_action for task %s/%s", taskID_names[t->type],
                subtaskID_names[t->subtype]);
#endif
          return task_action_none;
          break;
      }
      break;

    case task_type_kick1:
    case task_type_kick2:
    case task_type_logger:
    case task_type_fof_self:
    case task_type_fof_pair:
    case task_type_timestep:
    case task_type_send:
    case task_type_recv:
      if (t->ci->hydro.count > 0 && t->ci->grav.count > 0)
        return task_action_all;
      else if (t->ci->hydro.count > 0)
        return task_action_part;
      else if (t->ci->grav.count > 0)
        return task_action_gpart;
      else {
#ifdef SWIFT_DEBUG_CHECKS
        error("Task without particles");
#endif
      }
      break;

    case task_type_init_grav:
    case task_type_grav_mm:
    case task_type_grav_long_range:
      return task_action_multipole;
      break;

    case task_type_drift_gpart:
    case task_type_grav_down:
    case task_type_end_grav_force:
    case task_type_grav_mesh:
      return task_action_gpart;
      break;

    default:
#ifdef SWIFT_DEBUG_CHECKS
      error("Unknown task_action for task %s/%s", taskID_names[t->type],
            subtaskID_names[t->subtype]);
#endif
      return task_action_none;
      break;
  }

#ifdef SWIFT_DEBUG_CHECKS
  error("Unknown task_action for task %s/%s", taskID_names[t->type],
        subtaskID_names[t->subtype]);
#endif
  /* Silence compiler warnings. We should never get here. */
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
  const int ta_spart =
      (ta_act == task_action_spart || ta_act == task_action_all);
  const int tb_part = (tb_act == task_action_part || tb_act == task_action_all);
  const int tb_gpart =
      (tb_act == task_action_gpart || tb_act == task_action_all);
  const int tb_spart =
      (tb_act == task_action_spart || tb_act == task_action_all);

  /* In the case where both tasks act on parts */
  if (ta_part && tb_part) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->hydro.count;
    if (ta->cj != NULL) size_union += ta->cj->hydro.count;
    if (tb->ci != NULL) size_union += tb->ci->hydro.count;
    if (tb->cj != NULL) size_union += tb->cj->hydro.count;

    if (size_union == 0) return 0.f;

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
    if (ta->ci != NULL) size_union += ta->ci->grav.count;
    if (ta->cj != NULL) size_union += ta->cj->grav.count;
    if (tb->ci != NULL) size_union += tb->ci->grav.count;
    if (tb->cj != NULL) size_union += tb->cj->grav.count;

    if (size_union == 0) return 0.f;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_gpart(ta->ci, tb->ci) +
                                  task_cell_overlap_gpart(ta->ci, tb->cj) +
                                  task_cell_overlap_gpart(ta->cj, tb->ci) +
                                  task_cell_overlap_gpart(ta->cj, tb->cj);

    return ((float)size_intersect) / (size_union - size_intersect);
  }

  /* In the case where both tasks act on sparts */
  else if (ta_spart && tb_spart) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->stars.count;
    if (ta->cj != NULL) size_union += ta->cj->stars.count;
    if (tb->ci != NULL) size_union += tb->ci->stars.count;
    if (tb->cj != NULL) size_union += tb->cj->stars.count;

    if (size_union == 0) return 0.f;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_spart(ta->ci, tb->ci) +
                                  task_cell_overlap_spart(ta->ci, tb->cj) +
                                  task_cell_overlap_spart(ta->cj, tb->ci) +
                                  task_cell_overlap_spart(ta->cj, tb->cj);

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
    case task_type_logger:
    case task_type_timestep:
      cell_unlocktree(ci);
      cell_gunlocktree(ci);
      break;

    case task_type_drift_part:
    case task_type_sort:
    case task_type_ghost:
    case task_type_end_hydro_force:
    case task_type_timestep_limiter:
      cell_unlocktree(ci);
      break;

    case task_type_drift_gpart:
    case task_type_grav_mesh:
    case task_type_end_grav_force:
      cell_gunlocktree(ci);
      break;

    case task_type_stars_sort:
      cell_sunlocktree(ci);
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
        cell_gunlocktree(ci);
        cell_munlocktree(ci);
      } else if (subtype == task_subtype_stars_density) {
        cell_sunlocktree(ci);
      } else if (subtype == task_subtype_stars_feedback) {
        cell_sunlocktree(ci);
        cell_unlocktree(ci);
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
      } else if (subtype == task_subtype_stars_density) {
        cell_sunlocktree(ci);
        cell_sunlocktree(cj);
      } else if (subtype == task_subtype_stars_feedback) {
        cell_sunlocktree(ci);
        cell_sunlocktree(cj);
        cell_unlocktree(ci);
        cell_unlocktree(cj);
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
      cell_munlocktree(ci);
      break;

    case task_type_grav_mm:
      cell_munlocktree(ci);
      cell_munlocktree(cj);
      break;

    case task_type_star_formation:
      cell_unlocktree(ci);
      cell_sunlocktree(ci);
      cell_gunlocktree(ci);
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
        error(
            "Failed to test request on send/recv task (type=%s/%s tag=%lld, "
            "%s).",
            taskID_names[t->type], subtaskID_names[t->subtype], t->flags, buff);
      }
      return res;
#else
      error("SWIFT was not compiled with MPI support.");
#endif
      break;

    case task_type_kick1:
    case task_type_kick2:
    case task_type_logger:
    case task_type_timestep:
      if (ci->hydro.hold || ci->grav.phold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      if (cell_glocktree(ci) != 0) {
        cell_unlocktree(ci);
        return 0;
      }
      break;

    case task_type_drift_part:
    case task_type_sort:
    case task_type_ghost:
    case task_type_end_hydro_force:
    case task_type_timestep_limiter:
      if (ci->hydro.hold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      break;

    case task_type_stars_sort:
      if (ci->stars.hold) return 0;
      if (cell_slocktree(ci) != 0) return 0;
      break;

    case task_type_drift_gpart:
    case task_type_end_grav_force:
    case task_type_grav_mesh:
      if (ci->grav.phold) return 0;
      if (cell_glocktree(ci) != 0) return 0;
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
        /* Lock the gparts and the m-pole */
        if (ci->grav.phold || ci->grav.mhold) return 0;
        if (cell_glocktree(ci) != 0)
          return 0;
        else if (cell_mlocktree(ci) != 0) {
          cell_gunlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_stars_density) {
        if (ci->stars.hold) return 0;
        if (cell_slocktree(ci) != 0) return 0;
      } else if (subtype == task_subtype_stars_feedback) {
        if (ci->stars.hold) return 0;
        if (ci->hydro.hold) return 0;
        if (cell_slocktree(ci) != 0) return 0;
        if (cell_locktree(ci) != 0) {
          cell_sunlocktree(ci);
          return 0;
        }
      } else { /* subtype == hydro */
        if (ci->hydro.hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
        /* Lock the gparts and the m-pole in both cells */
        if (ci->grav.phold || cj->grav.phold) return 0;
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
      } else if (subtype == task_subtype_stars_density) {
        if (ci->stars.hold || cj->stars.hold) return 0;
        if (cell_slocktree(ci) != 0) return 0;
        if (cell_slocktree(cj) != 0) {
          cell_sunlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_stars_feedback) {
        /* Lock the stars and the gas particles in both cells */
        if (ci->stars.hold || cj->stars.hold) return 0;
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (cell_slocktree(ci) != 0) return 0;
        if (cell_slocktree(cj) != 0) {
          cell_sunlocktree(ci);
          return 0;
        }
        if (cell_locktree(ci) != 0) {
          cell_sunlocktree(ci);
          cell_sunlocktree(cj);
          return 0;
        }
        if (cell_locktree(cj) != 0) {
          cell_sunlocktree(ci);
          cell_sunlocktree(cj);
          cell_unlocktree(ci);
          return 0;
        }
      } else { /* subtype == hydro */
        /* Lock the parts in both cells */
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
        if (cell_locktree(cj) != 0) {
          cell_unlocktree(ci);
          return 0;
        }
      }
      break;

    case task_type_grav_down:
      /* Lock the gparts and the m-poles */
      if (ci->grav.phold || ci->grav.mhold) return 0;
      if (cell_glocktree(ci) != 0)
        return 0;
      else if (cell_mlocktree(ci) != 0) {
        cell_gunlocktree(ci);
        return 0;
      }
      break;

    case task_type_grav_long_range:
      /* Lock the m-poles */
      if (ci->grav.mhold) return 0;
      if (cell_mlocktree(ci) != 0) return 0;
      break;

    case task_type_grav_mm:
      /* Lock both m-poles */
      if (ci->grav.mhold || cj->grav.mhold) return 0;
      if (cell_mlocktree(ci) != 0) return 0;
      if (cell_mlocktree(cj) != 0) {
        cell_munlocktree(ci);
        return 0;
      }
      break;

    case task_type_star_formation:
      /* Lock the gas, gravity and star particles */
      if (ci->hydro.hold || ci->stars.hold || ci->grav.phold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      if (cell_slocktree(ci) != 0) {
        cell_unlocktree(ci);
        return 0;
      }
      if (cell_glocktree(ci) != 0) {
        cell_unlocktree(ci);
        cell_sunlocktree(ci);
        return 0;
      }

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

/**
 * @brief Get the group name of a task.
 *
 * This is used to group tasks with similar actions in the task dependency
 * graph.
 *
 * @param type The #task type.
 * @param subtype The #task subtype.
 * @param cluster (return) The group name (should be allocated)
 */
void task_get_group_name(int type, int subtype, char *cluster) {

  if (type == task_type_grav_long_range || type == task_type_grav_mm ||
      type == task_type_grav_mesh) {

    strcpy(cluster, "Gravity");
    return;
  }

  switch (subtype) {
    case task_subtype_density:
      strcpy(cluster, "Density");
      break;
    case task_subtype_gradient:
      if (type == task_type_send || type == task_type_recv) {
        strcpy(cluster, "None");
      } else {
        strcpy(cluster, "Gradient");
      }
      break;
    case task_subtype_force:
      strcpy(cluster, "Force");
      break;
    case task_subtype_grav:
      strcpy(cluster, "Gravity");
      break;
    case task_subtype_limiter:
      strcpy(cluster, "Timestep_limiter");
      break;
    case task_subtype_stars_density:
      strcpy(cluster, "StarsDensity");
      break;
    case task_subtype_stars_feedback:
      strcpy(cluster, "StarsFeedback");
      break;
    default:
      strcpy(cluster, "None");
      break;
  }
}

/**
 * @brief Generate the full name of a #task.
 *
 * @param type The #task type.
 * @param subtype The #task type.
 * @param name (return) The formatted string
 */
void task_get_full_name(int type, int subtype, char *name) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Check input */
  if (type >= task_type_count) error("Unknown task type %i", type);

  if (subtype >= task_subtype_count)
    error("Unknown task subtype %i with type %s", subtype, taskID_names[type]);
#endif

  /* Full task name */
  if (subtype == task_subtype_none)
    sprintf(name, "%s", taskID_names[type]);
  else
    sprintf(name, "%s_%s", taskID_names[type], subtaskID_names[subtype]);
}

#ifdef WITH_MPI
/**
 * @brief Create global communicators for each of the subtasks.
 */
void task_create_mpi_comms(void) {
  for (int i = 0; i < task_subtype_count; i++) {
    MPI_Comm_dup(MPI_COMM_WORLD, &subtaskMPI_comms[i]);
  }
}
#endif

/**
 * @brief dump all the tasks of all the known engines into a file for
 * postprocessing.
 *
 * Dumps the information to a file "thread_info-stepn.dat" where n is the
 * given step value, or "thread_info_MPI-stepn.dat", if we are running
 * under MPI. Note if running under MPIU all the ranks are dumped into this
 * one file, which has an additional field to identify the rank.
 *
 * @param e the #engine
 * @param step the current step.
 */
void task_dump_all(struct engine *e, int step) {

#ifdef SWIFT_DEBUG_TASKS

  /* Need this to convert ticks to seconds. */
  unsigned long long cpufreq = clocks_get_cpufreq();

#ifdef WITH_MPI
  /* Make sure output file is empty, only on one rank. */
  char dumpfile[35];
  snprintf(dumpfile, sizeof(dumpfile), "thread_info_MPI-step%d.dat", step);
  FILE *file_thread;
  if (engine_rank == 0) {
    file_thread = fopen(dumpfile, "w");
    fclose(file_thread);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < e->nr_nodes; i++) {

    /* Rank 0 decides the index of the writing node, this happens
     * one-by-one. */
    int kk = i;
    MPI_Bcast(&kk, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (i == engine_rank) {

      /* Open file and position at end. */
      file_thread = fopen(dumpfile, "a");

      /* Add some information to help with the plots and conversion of ticks to
       * seconds. */
      fprintf(file_thread, " %03d 0 0 0 0 %lld %lld %lld %lld %lld 0 0 %lld\n",
              engine_rank, (long long int)e->tic_step,
              (long long int)e->toc_step, e->updates, e->g_updates,
              e->s_updates, cpufreq);
      int count = 0;
      for (int l = 0; l < e->sched.nr_tasks; l++) {
        if (!e->sched.tasks[l].implicit && e->sched.tasks[l].toc != 0) {
          fprintf(
              file_thread, " %03i %i %i %i %i %lli %lli %i %i %i %i %lli %i\n",
              engine_rank, e->sched.tasks[l].rid, e->sched.tasks[l].type,
              e->sched.tasks[l].subtype, (e->sched.tasks[l].cj == NULL),
              (long long int)e->sched.tasks[l].tic,
              (long long int)e->sched.tasks[l].toc,
              (e->sched.tasks[l].ci != NULL) ? e->sched.tasks[l].ci->hydro.count
                                             : 0,
              (e->sched.tasks[l].cj != NULL) ? e->sched.tasks[l].cj->hydro.count
                                             : 0,
              (e->sched.tasks[l].ci != NULL) ? e->sched.tasks[l].ci->grav.count
                                             : 0,
              (e->sched.tasks[l].cj != NULL) ? e->sched.tasks[l].cj->grav.count
                                             : 0,
              e->sched.tasks[l].flags, e->sched.tasks[l].sid);
        }
        count++;
      }
      fclose(file_thread);
    }

    /* And we wait for all to synchronize. */
    MPI_Barrier(MPI_COMM_WORLD);
  }

#else
  /* Non-MPI, so just a single engine's worth of tasks to dump. */
  char dumpfile[32];
  snprintf(dumpfile, sizeof(dumpfile), "thread_info-step%d.dat", step);
  FILE *file_thread;
  file_thread = fopen(dumpfile, "w");

  /* Add some information to help with the plots and conversion of ticks to
   * seconds. */
  fprintf(file_thread, " %d %d %d %d %lld %lld %lld %lld %lld %d %lld\n", -2,
          -1, -1, 1, (unsigned long long)e->tic_step,
          (unsigned long long)e->toc_step, e->updates, e->g_updates,
          e->s_updates, 0, cpufreq);
  for (int l = 0; l < e->sched.nr_tasks; l++) {
    if (!e->sched.tasks[l].implicit && e->sched.tasks[l].toc != 0) {
      fprintf(
          file_thread, " %i %i %i %i %lli %lli %i %i %i %i %i\n",
          e->sched.tasks[l].rid, e->sched.tasks[l].type,
          e->sched.tasks[l].subtype, (e->sched.tasks[l].cj == NULL),
          (unsigned long long)e->sched.tasks[l].tic,
          (unsigned long long)e->sched.tasks[l].toc,
          (e->sched.tasks[l].ci == NULL) ? 0
                                         : e->sched.tasks[l].ci->hydro.count,
          (e->sched.tasks[l].cj == NULL) ? 0
                                         : e->sched.tasks[l].cj->hydro.count,
          (e->sched.tasks[l].ci == NULL) ? 0 : e->sched.tasks[l].ci->grav.count,
          (e->sched.tasks[l].cj == NULL) ? 0 : e->sched.tasks[l].cj->grav.count,
          e->sched.tasks[l].sid);
    }
  }
  fclose(file_thread);
#endif  // WITH_MPI
#endif  // SWIFT_DEBUG_TASKS
}

/**
 * @brief Generate simple statistics about the times used by the tasks of
 *        all the engines and write these into two format, a human readable
 *        version for debugging and one intented for inclusion as the fixed
 *        costs for repartitioning.
 *
 * Note that when running under MPI all the tasks can be summed into this single
 * file. In the fuller, human readable file, the statistics included are the
 * number of task of each type/subtype followed by the minimum, maximum, mean
 * and total time, in millisec and then the fixed costs value.
 *
 * If header is set, only the fixed costs value is written into the output
 * file in a format that is suitable for inclusion in SWIFT (as
 * partition_fixed_costs.h).
 *
 * @param dumpfile name of the file for the output.
 * @param e the #engine
 * @param header whether to write a header include file.
 * @param allranks do the statistics over all ranks, if not just the current
 *                 one, only used if header is false.
 */
void task_dump_stats(const char *dumpfile, struct engine *e, int header,
                     int allranks) {

  /* Need arrays for sum, min and max across all types and subtypes. */
  double sum[task_type_count][task_subtype_count];
  double min[task_type_count][task_subtype_count];
  double max[task_type_count][task_subtype_count];
  int count[task_type_count][task_subtype_count];

  for (int j = 0; j < task_type_count; j++) {
    for (int k = 0; k < task_subtype_count; k++) {
      sum[j][k] = 0.0;
      count[j][k] = 0;
      min[j][k] = DBL_MAX;
      max[j][k] = 0.0;
    }
  }

  double total[1] = {0.0};
  for (int l = 0; l < e->sched.nr_tasks; l++) {
    int type = e->sched.tasks[l].type;

    /* Skip implicit tasks, tasks that didn't run and MPI send/recv as these
     * are not interesting (or meaningfully measured). */
    if (!e->sched.tasks[l].implicit && e->sched.tasks[l].toc != 0 &&
        type != task_type_send && type != task_type_recv) {
      int subtype = e->sched.tasks[l].subtype;

      double dt = e->sched.tasks[l].toc - e->sched.tasks[l].tic;
      sum[type][subtype] += dt;
      count[type][subtype] += 1;
      if (dt < min[type][subtype]) {
        min[type][subtype] = dt;
      }
      if (dt > max[type][subtype]) {
        max[type][subtype] = dt;
      }
      total[0] += dt;
    }
  }

#ifdef WITH_MPI
  if (allranks || header) {
    /* Get these from all ranks for output from rank 0. Could wrap these into a
     * single operation. */
    size_t size = task_type_count * task_subtype_count;
    int res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : sum), sum, size,
                         MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task sums");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : count), count, size,
                     MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task counts");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : min), min, size,
                     MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task minima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : max), max, size,
                     MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task maxima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : total), total, 1,
                     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task total time");
  }

  if (!allranks || (engine_rank == 0 && (allranks || header))) {
#endif

    FILE *dfile = fopen(dumpfile, "w");
    if (header) {
      fprintf(dfile, "/* use as src/partition_fixed_costs.h */\n");
      fprintf(dfile, "#define HAVE_FIXED_COSTS 1\n");
    } else {
      fprintf(dfile, "# task ntasks min max sum mean percent fixed_cost\n");
    }

    for (int j = 0; j < task_type_count; j++) {
      const char *taskID = taskID_names[j];
      for (int k = 0; k < task_subtype_count; k++) {
        if (sum[j][k] > 0.0) {
          double mean = sum[j][k] / (double)count[j][k];
          double perc = 100.0 * sum[j][k] / total[0];

          /* Fixed cost is in .1ns as we want to compare between runs in
           * some absolute units. */
          int fixed_cost = (int)(clocks_from_ticks(mean) * 10000.f);
          if (header) {
            fprintf(dfile, "repartition_costs[%d][%d] = %10d; /* %s/%s */\n", j,
                    k, fixed_cost, taskID, subtaskID_names[k]);
          } else {
            fprintf(dfile,
                    "%15s/%-10s %10d %14.4f %14.4f %14.4f %14.4f %14.4f %10d\n",
                    taskID, subtaskID_names[k], count[j][k],
                    clocks_from_ticks(min[j][k]), clocks_from_ticks(max[j][k]),
                    clocks_from_ticks(sum[j][k]), clocks_from_ticks(mean), perc,
                    fixed_cost);
          }
        }
      }
    }
    fclose(dfile);
#ifdef WITH_MPI
  }
#endif
}
