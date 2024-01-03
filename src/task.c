/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

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
#include "mpiuse.h"

/* Task type names. */
const char *taskID_names[task_type_count] = {
    "none",
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
    "drift_sink",
    "drift_bpart",
    "drift_gpart",
    "drift_gpart_out",
    "hydro_end_force",
    "kick1",
    "kick2",
    "timestep",
    "timestep_limiter",
    "timestep_sync",
    "collect",
    "send",
    "recv",
    "pack",
    "unpack",
    "grav_long_range",
    "grav_mm",
    "grav_down_in",
    "grav_down",
    "grav_end_force",
    "cooling",
    "cooling_in",
    "cooling_out",
    "star_formation",
    "star_formation_in",
    "star_formation_out",
    "star_formation_sink",
    "csds",
    "stars_in",
    "stars_out",
    "stars_ghost_in",
    "stars_density_ghost",
    "stars_ghost_out",
    "stars_prep_ghost1",
    "hydro_prep_ghost1",
    "stars_prep_ghost2",
    "stars_sort",
    "stars_resort",
    "bh_in",
    "bh_out",
    "bh_density_ghost",
    "bh_swallow_ghost1",
    "bh_swallow_ghost2",
    "bh_swallow_ghost3",
    "fof_self",
    "fof_pair",
    "fof_attach_self",
    "fof_attach_pair",
    "neutrino_weight",
    "sink_in",
    "sink_ghost1",
    "sink_ghost2",
    "sink_out",
    "rt_in",
    "rt_out",
    "sink_formation",
    "rt_ghost1",
    "rt_ghost2",
    "rt_transport_out",
    "rt_tchem",
    "rt_advance_cell_time",
    "rt_sorts",
    "rt_collect_times",
};

/* Sub-task type names. */
const char *subtaskID_names[task_subtype_count] = {
    "none",
    "density",
    "gradient",
    "force",
    "limiter",
    "grav",
    "external_grav",
    "tend",
    "xv",
    "rho",
    "part_swallow",
    "bpart_merger",
    "gpart",
    "spart_density",
    "part_prep1",
    "spart_prep2",
    "stars_density",
    "stars_prep1",
    "stars_prep2",
    "stars_feedback",
    "sf_counts",
    "bpart_rho",
    "bpart_feedback",
    "bh_density",
    "bh_swallow",
    "do_gas_swallow",
    "do_bh_swallow",
    "bh_feedback",
    "sink_do_sink_swallow",
    "sink_swallow",
    "sink_do_gas_swallow",
    "rt_gradient",
    "rt_transport",
};

const char *task_category_names[task_category_count] = {
    "drift",       "sorts",        "resort",
    "hydro",       "gravity",      "feedback",
    "black holes", "cooling",      "star formation",
    "limiter",     "sync",         "time integration",
    "mpi",         "pack",         "fof",
    "others",      "neutrino",     "sink",
    "RT",          "CSDS",         "RT tchem",
    "Hydro ghost", "Hydro density"};

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
#define TASK_CELL_OVERLAP(TYPE, ARRAY, COUNT)                           \
  __attribute__((always_inline))                                        \
  INLINE static size_t task_cell_overlap_##TYPE(                        \
      const struct cell *restrict ci, const struct cell *restrict cj) { \
                                                                        \
    if (ci == NULL || cj == NULL) return 0;                             \
                                                                        \
    if (ci->ARRAY <= cj->ARRAY &&                                       \
        ci->ARRAY + ci->COUNT >= cj->ARRAY + cj->COUNT) {               \
      return cj->COUNT;                                                 \
    } else if (cj->ARRAY <= ci->ARRAY &&                                \
               cj->ARRAY + cj->COUNT >= ci->ARRAY + ci->COUNT) {        \
      return ci->COUNT;                                                 \
    }                                                                   \
                                                                        \
    return 0;                                                           \
  }

TASK_CELL_OVERLAP(part, hydro.parts, hydro.count);
TASK_CELL_OVERLAP(gpart, grav.parts, grav.count);
TASK_CELL_OVERLAP(spart, stars.parts, stars.count);
TASK_CELL_OVERLAP(sink, sinks.parts, sinks.count);
TASK_CELL_OVERLAP(bpart, black_holes.parts, black_holes.count);

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
    case task_type_end_hydro_force:
      return task_action_part;
      break;

    case task_type_star_formation:
    case task_type_star_formation_sink:
    case task_type_sink_formation:
      return task_action_all;

    case task_type_drift_spart:
    case task_type_stars_ghost:
    case task_type_stars_sort:
    case task_type_stars_resort:
      return task_action_spart;
      break;

    case task_type_drift_sink:
      return task_action_sink;
      break;

    case task_type_drift_bpart:
    case task_type_bh_density_ghost:
    case task_type_bh_swallow_ghost3:
      return task_action_bpart;
      break;

    case task_type_rt_ghost1:
    case task_type_rt_ghost2:
    case task_type_rt_tchem:
    case task_type_rt_sort:
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
        case task_subtype_limiter:
          return task_action_part;
          break;

        case task_subtype_stars_density:
        case task_subtype_stars_feedback:
          return task_action_all;
          break;

        case task_subtype_bh_density:
        case task_subtype_bh_feedback:
        case task_subtype_bh_swallow:
        case task_subtype_do_gas_swallow:
          return task_action_all;
          break;

        case task_subtype_do_bh_swallow:
          return task_action_bpart;
          break;

        case task_subtype_sink_do_gas_swallow:
        case task_subtype_sink_do_sink_swallow:
        case task_subtype_sink_swallow:
          return task_action_all;

        case task_subtype_rt_transport:
        case task_subtype_rt_gradient:
          return task_action_part;
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
    case task_type_csds:
    case task_type_fof_self:
    case task_type_fof_pair:
    case task_type_fof_attach_self:
    case task_type_fof_attach_pair:
    case task_type_timestep:
    case task_type_timestep_limiter:
    case task_type_timestep_sync:
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
  const int ta_sink = (ta_act == task_action_sink || ta_act == task_action_all);
  const int ta_bpart =
      (ta_act == task_action_bpart || ta_act == task_action_all);
  const int tb_part = (tb_act == task_action_part || tb_act == task_action_all);
  const int tb_gpart =
      (tb_act == task_action_gpart || tb_act == task_action_all);
  const int tb_spart =
      (tb_act == task_action_spart || tb_act == task_action_all);
  const int tb_sink = (tb_act == task_action_sink || tb_act == task_action_all);
  const int tb_bpart =
      (tb_act == task_action_bpart || tb_act == task_action_all);

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

  /* In the case where both tasks act on sink */
  else if (ta_sink && tb_sink) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->sinks.count;
    if (ta->cj != NULL) size_union += ta->cj->sinks.count;
    if (tb->ci != NULL) size_union += tb->ci->sinks.count;
    if (tb->cj != NULL) size_union += tb->cj->sinks.count;

    if (size_union == 0) return 0.f;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_spart(ta->ci, tb->ci) +
                                  task_cell_overlap_sink(ta->ci, tb->cj) +
                                  task_cell_overlap_sink(ta->cj, tb->ci) +
                                  task_cell_overlap_sink(ta->cj, tb->cj);

    return ((float)size_intersect) / (size_union - size_intersect);
  }

  /* In the case where both tasks act on bparts */
  else if (ta_bpart && tb_bpart) {

    /* Compute the union of the cell data. */
    size_t size_union = 0;
    if (ta->ci != NULL) size_union += ta->ci->black_holes.count;
    if (ta->cj != NULL) size_union += ta->cj->black_holes.count;
    if (tb->ci != NULL) size_union += tb->ci->black_holes.count;
    if (tb->cj != NULL) size_union += tb->cj->black_holes.count;

    if (size_union == 0) return 0.f;

    /* Compute the intersection of the cell data. */
    const size_t size_intersect = task_cell_overlap_bpart(ta->ci, tb->ci) +
                                  task_cell_overlap_bpart(ta->ci, tb->cj) +
                                  task_cell_overlap_bpart(ta->cj, tb->ci) +
                                  task_cell_overlap_bpart(ta->cj, tb->cj);

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
    case task_type_csds:
    case task_type_timestep:
      cell_unlocktree(ci);
      cell_gunlocktree(ci);
      break;

    case task_type_drift_part:
    case task_type_sort:
    case task_type_ghost:
    case task_type_extra_ghost:
    case task_type_end_hydro_force:
    case task_type_timestep_limiter:
    case task_type_timestep_sync:
    case task_type_rt_ghost1:
    case task_type_rt_ghost2:
    case task_type_rt_tchem:
    case task_type_rt_sort:
    case task_type_rt_advance_cell_time:
      cell_unlocktree(ci);
      break;

    case task_type_drift_gpart:
    case task_type_end_grav_force:
      cell_gunlocktree(ci);
      break;

    case task_type_drift_sink:
      cell_sink_unlocktree(ci);
      break;

    case task_type_stars_sort:
    case task_type_stars_resort:
      cell_sunlocktree(ci);
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        cell_gunlocktree(ci);
        cell_munlocktree(ci);
#endif
      } else if (subtype == task_subtype_sink_swallow) {
        cell_sink_unlocktree(ci);
        cell_unlocktree(ci);
      } else if (subtype == task_subtype_sink_do_sink_swallow) {
        cell_sink_unlocktree(ci);
        cell_gunlocktree(ci);
      } else if (subtype == task_subtype_sink_do_gas_swallow) {
        cell_unlocktree(ci);
        cell_sink_unlocktree(ci);
        cell_gunlocktree(ci);
      } else if ((subtype == task_subtype_stars_density) ||
                 (subtype == task_subtype_stars_prep1) ||
                 (subtype == task_subtype_stars_prep2) ||
                 (subtype == task_subtype_stars_feedback)) {
        cell_sunlocktree(ci);
        cell_unlocktree(ci);
      } else if ((subtype == task_subtype_bh_density) ||
                 (subtype == task_subtype_bh_feedback) ||
                 (subtype == task_subtype_bh_swallow) ||
                 (subtype == task_subtype_do_gas_swallow)) {
        cell_bunlocktree(ci);
        cell_unlocktree(ci);
      } else if (subtype == task_subtype_do_bh_swallow) {
        cell_bunlocktree(ci);
      } else if (subtype == task_subtype_limiter) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        cell_unlocktree(ci);
#endif
      } else { /* hydro */
        cell_unlocktree(ci);
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        cell_gunlocktree(ci);
        cell_gunlocktree(cj);
        cell_munlocktree(ci);
        cell_munlocktree(cj);
#endif
      } else if (subtype == task_subtype_sink_swallow) {
        cell_sink_unlocktree(ci);
        cell_sink_unlocktree(cj);
        cell_unlocktree(ci);
        cell_unlocktree(cj);
      } else if (subtype == task_subtype_sink_do_sink_swallow) {
        cell_sink_unlocktree(ci);
        cell_sink_unlocktree(cj);
        cell_gunlocktree(ci);
        cell_gunlocktree(cj);
      } else if (subtype == task_subtype_sink_do_gas_swallow) {
        cell_sink_unlocktree(ci);
        cell_sink_unlocktree(cj);
        cell_unlocktree(ci);
        cell_unlocktree(cj);
        cell_gunlocktree(ci);
        cell_gunlocktree(cj);
      } else if ((subtype == task_subtype_stars_density) ||
                 (subtype == task_subtype_stars_prep1) ||
                 (subtype == task_subtype_stars_prep2) ||
                 (subtype == task_subtype_stars_feedback)) {
        cell_sunlocktree(ci);
        cell_sunlocktree(cj);
        cell_unlocktree(ci);
        cell_unlocktree(cj);
      } else if ((subtype == task_subtype_bh_density) ||
                 (subtype == task_subtype_bh_feedback) ||
                 (subtype == task_subtype_bh_swallow) ||
                 (subtype == task_subtype_do_gas_swallow)) {
        cell_bunlocktree(ci);
        cell_bunlocktree(cj);
        cell_unlocktree(ci);
        cell_unlocktree(cj);
      } else if (subtype == task_subtype_do_bh_swallow) {
        cell_bunlocktree(ci);
        cell_bunlocktree(cj);
      } else if (subtype == task_subtype_limiter) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        cell_unlocktree(ci);
        cell_unlocktree(cj);
#endif
      } else { /* hydro */
        cell_unlocktree(ci);
        cell_unlocktree(cj);
      }
      break;

    case task_type_grav_down:
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      cell_gunlocktree(ci);
      cell_munlocktree(ci);
#endif
      break;

    case task_type_grav_long_range:
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      cell_munlocktree(ci);
#endif
      break;

    case task_type_grav_mm:
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      cell_munlocktree(ci);
      cell_munlocktree(cj);
#endif
      break;

    case task_type_fof_self:
    case task_type_fof_attach_self:
      cell_gunlocktree(ci);
      break;

    case task_type_fof_pair:
    case task_type_fof_attach_pair:
      cell_gunlocktree(ci);
      cell_gunlocktree(cj);
      break;

    case task_type_star_formation:
      cell_unlocktree(ci);
      cell_sunlocktree(ci);
      cell_gunlocktree(ci);
      break;

    case task_type_star_formation_sink:
      cell_sink_unlocktree(ci);
      cell_sunlocktree(ci);
      cell_gunlocktree(ci);
      break;

    case task_type_sink_formation:
      cell_unlocktree(ci);
      cell_sink_unlocktree(ci);
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

      /* And log deactivation, if logging enabled. */
      if (res) {
        mpiuse_log_allocation(t->type, t->subtype, &t->req, 0, 0, 0, 0);
      }

      return res;
#else
      error("SWIFT was not compiled with MPI support.");
#endif
      break;

    case task_type_kick1:
    case task_type_kick2:
    case task_type_csds:
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
    case task_type_extra_ghost:
    case task_type_end_hydro_force:
    case task_type_timestep_limiter:
    case task_type_timestep_sync:
    case task_type_rt_ghost1:
    case task_type_rt_ghost2:
    case task_type_rt_tchem:
    case task_type_rt_sort:
    case task_type_rt_advance_cell_time:
      if (ci->hydro.hold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      break;

    case task_type_stars_sort:
    case task_type_stars_resort:
      if (ci->stars.hold) return 0;
      if (cell_slocktree(ci) != 0) return 0;
      break;

    case task_type_drift_gpart:
    case task_type_end_grav_force:
      if (ci->grav.phold) return 0;
      if (cell_glocktree(ci) != 0) return 0;
      break;

    case task_type_drift_sink:
      if (ci->sinks.hold) return 0;
      if (cell_sink_locktree(ci) != 0) return 0;
      break;

    case task_type_self:
    case task_type_sub_self:
      if (subtype == task_subtype_grav) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        /* Lock the gparts and the m-pole */
        if (ci->grav.phold || ci->grav.mhold) return 0;
        if (cell_glocktree(ci) != 0)
          return 0;
        else if (cell_mlocktree(ci) != 0) {
          cell_gunlocktree(ci);
          return 0;
        }
#endif
      } else if (subtype == task_subtype_sink_do_sink_swallow) {
        if (ci->sinks.hold) return 0;
        if (ci->grav.phold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_glocktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_sink_swallow) {
        if (ci->sinks.hold) return 0;
        if (ci->hydro.hold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_locktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_sink_do_gas_swallow) {
        if (ci->sinks.hold) return 0;
        if (ci->grav.phold) return 0;
        if (ci->hydro.hold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_locktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
        if (cell_glocktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          cell_unlocktree(ci);
          return 0;
        }
      } else if ((subtype == task_subtype_stars_density) ||
                 (subtype == task_subtype_stars_prep1) ||
                 (subtype == task_subtype_stars_prep2) ||
                 (subtype == task_subtype_stars_feedback)) {
        if (ci->stars.hold) return 0;
        if (ci->hydro.hold) return 0;
        if (cell_slocktree(ci) != 0) return 0;
        if (cell_locktree(ci) != 0) {
          cell_sunlocktree(ci);
          return 0;
        }
      } else if ((subtype == task_subtype_bh_density) ||
                 (subtype == task_subtype_bh_feedback) ||
                 (subtype == task_subtype_bh_swallow) ||
                 (subtype == task_subtype_do_gas_swallow)) {
        if (ci->black_holes.hold) return 0;
        if (ci->hydro.hold) return 0;
        if (cell_blocktree(ci) != 0) return 0;
        if (cell_locktree(ci) != 0) {
          cell_bunlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_do_bh_swallow) {
        if (ci->black_holes.hold) return 0;
        if (cell_blocktree(ci) != 0) return 0;
      } else if (subtype == task_subtype_limiter) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        if (ci->hydro.hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
#endif
      } else { /* subtype == hydro */
        if (ci->hydro.hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
      }
      break;

    case task_type_pair:
    case task_type_sub_pair:
      if (subtype == task_subtype_grav) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
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
#endif
      } else if (subtype == task_subtype_sink_swallow) {
        /* Lock the sinks and the gas particles in both cells */
        if (ci->sinks.hold || cj->sinks.hold) return 0;
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_sink_locktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
        if (cell_locktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          return 0;
        }
        if (cell_locktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          cell_unlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_sink_do_gas_swallow) {
        /* Lock the sinks and the gas particles in both cells */
        if (ci->sinks.hold || cj->sinks.hold) return 0;
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (ci->grav.phold || cj->grav.phold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_sink_locktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
        if (cell_locktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          return 0;
        }
        if (cell_locktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          cell_unlocktree(ci);
          return 0;
        }
        if (cell_glocktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          cell_unlocktree(ci);
          cell_unlocktree(cj);
          return 0;
        }
        if (cell_glocktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          cell_unlocktree(ci);
          cell_unlocktree(cj);
          cell_gunlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_sink_do_sink_swallow) {
        /* Lock the sink and the dm particles in both cells */
        if (ci->sinks.hold || cj->sinks.hold) return 0;
        if (ci->grav.phold || cj->grav.phold) return 0;
        if (cell_sink_locktree(ci) != 0) return 0;
        if (cell_sink_locktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          return 0;
        }
        if (cell_glocktree(ci) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          return 0;
        }
        if (cell_glocktree(cj) != 0) {
          cell_sink_unlocktree(ci);
          cell_sink_unlocktree(cj);
          cell_gunlocktree(ci);
          return 0;
        }
      } else if ((subtype == task_subtype_stars_density) ||
                 (subtype == task_subtype_stars_prep1) ||
                 (subtype == task_subtype_stars_prep2) ||
                 (subtype == task_subtype_stars_feedback)) {
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
      } else if ((subtype == task_subtype_bh_density) ||
                 (subtype == task_subtype_bh_feedback) ||
                 (subtype == task_subtype_bh_swallow) ||
                 (subtype == task_subtype_do_gas_swallow)) {
        /* Lock the BHs and the gas particles in both cells */
        if (ci->black_holes.hold || cj->black_holes.hold) return 0;
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (cell_blocktree(ci) != 0) return 0;
        if (cell_blocktree(cj) != 0) {
          cell_bunlocktree(ci);
          return 0;
        }
        if (cell_locktree(ci) != 0) {
          cell_bunlocktree(ci);
          cell_bunlocktree(cj);
          return 0;
        }
        if (cell_locktree(cj) != 0) {
          cell_bunlocktree(ci);
          cell_bunlocktree(cj);
          cell_unlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_do_bh_swallow) {
        if (ci->black_holes.hold || cj->black_holes.hold) return 0;
        if (cell_blocktree(ci) != 0) return 0;
        if (cell_blocktree(cj) != 0) {
          cell_bunlocktree(ci);
          return 0;
        }
      } else if (subtype == task_subtype_limiter) {
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
        if (ci->hydro.hold || cj->hydro.hold) return 0;
        if (cell_locktree(ci) != 0) return 0;
        if (cell_locktree(cj) != 0) {
          cell_unlocktree(ci);
          return 0;
        }
#endif
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
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      /* Lock the gparts and the m-poles */
      if (ci->grav.phold || ci->grav.mhold) return 0;
      if (cell_glocktree(ci) != 0)
        return 0;
      else if (cell_mlocktree(ci) != 0) {
        cell_gunlocktree(ci);
        return 0;
      }
#endif
      break;

    case task_type_grav_long_range:
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      /* Lock the m-poles */
      if (ci->grav.mhold) return 0;
      if (cell_mlocktree(ci) != 0) return 0;
#endif
      break;

    case task_type_grav_mm:
#ifdef SWIFT_TASKS_WITHOUT_ATOMICS
      /* Lock both m-poles */
      if (ci->grav.mhold || cj->grav.mhold) return 0;
      if (cell_mlocktree(ci) != 0) return 0;
      if (cell_mlocktree(cj) != 0) {
        cell_munlocktree(ci);
        return 0;
      }
#endif
      break;

    case task_type_fof_self:
    case task_type_fof_attach_self:
      /* Lock the gpart as this this what we act on */
      if (ci->grav.phold) return 0;
      if (cell_glocktree(ci) != 0) return 0;
      break;

    case task_type_fof_pair:
    case task_type_fof_attach_pair:
      /* Lock the gpart as this this what we act on */
      if (ci->grav.phold || cj->grav.phold) return 0;
      if (cell_glocktree(ci) != 0) return 0;
      if (cell_glocktree(cj) != 0) {
        cell_gunlocktree(ci);
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
      break;

    case task_type_star_formation_sink:
      /* Lock the gas, gravity and star particles */
      if (ci->sinks.hold || ci->stars.hold || ci->grav.phold) return 0;
      if (cell_sink_locktree(ci) != 0) return 0;
      if (cell_slocktree(ci) != 0) {
        cell_sink_unlocktree(ci);
        return 0;
      }
      if (cell_glocktree(ci) != 0) {
        cell_sink_unlocktree(ci);
        cell_sunlocktree(ci);
        return 0;
      }
      break;

    case task_type_sink_formation:
      /* Lock the gas, gravity and star particles */
      if (ci->hydro.hold || ci->sinks.hold || ci->grav.phold) return 0;
      if (cell_locktree(ci) != 0) return 0;
      if (cell_sink_locktree(ci) != 0) {
        cell_unlocktree(ci);
        return 0;
      }
      if (cell_glocktree(ci) != 0) {
        cell_unlocktree(ci);
        cell_sink_unlocktree(ci);
        return 0;
      }
      break;

    default:
      break;
  }

  /* If we made it this far, we've got a lock. */
  return 1;
}

/**
 * @brief Returns a pointer to the unique task unlocked by this task.
 *
 * The task MUST have only dependence!
 *
 * @param The #task.
 */
struct task *task_get_unique_dependent(const struct task *t) {

#ifdef SWIFT_DEBUG_CHECKS
  if (t->nr_unlock_tasks != 1)
    error("Task is unlocking more than one dependence!");
#endif

  return t->unlock_tasks[0];
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

  if (type == task_type_grav_long_range || type == task_type_grav_mm) {

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
      if (type == task_type_send || type == task_type_recv) {
        strcpy(cluster, "None");
      } else {
        strcpy(cluster, "Timestep_limiter");
      }
      break;
    case task_subtype_stars_density:
      strcpy(cluster, "StarsDensity");
      break;
    case task_subtype_stars_prep1:
      strcpy(cluster, "StarsKickPrep1");
      break;
    case task_subtype_stars_prep2:
      strcpy(cluster, "StarsKickPrep2");
      break;
    case task_subtype_stars_feedback:
      strcpy(cluster, "StarsFeedback");
      break;
    case task_subtype_bh_density:
      strcpy(cluster, "BHDensity");
      break;
    case task_subtype_bh_swallow:
      strcpy(cluster, "BHSwallow");
      break;
    case task_subtype_do_gas_swallow:
      strcpy(cluster, "DoGasSwallow");
      break;
    case task_subtype_do_bh_swallow:
      strcpy(cluster, "DoBHSwallow");
      break;
    case task_subtype_bh_feedback:
      strcpy(cluster, "BHFeedback");
      break;
    case task_subtype_rt_gradient:
      if (type == task_type_send || type == task_type_recv) {
        strcpy(cluster, "None");
      } else {
        strcpy(cluster, "RTgradient");
      }
      break;
    case task_subtype_rt_transport:
      if (type == task_type_send || type == task_type_recv) {
        strcpy(cluster, "None");
      } else {
        strcpy(cluster, "RTtransport");
      }
      break;
    case task_subtype_sink_swallow:
      strcpy(cluster, "SinkFormation");
      break;
    case task_subtype_sink_do_sink_swallow:
      strcpy(cluster, "SinkMerger");
      break;
    case task_subtype_sink_do_gas_swallow:
      strcpy(cluster, "SinkAccretion");
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

void task_create_name_files(const char *file_prefix) {
  char file_name[200];
  sprintf(file_name, "%s_task_types.txt", file_prefix);
  FILE *file = fopen(file_name, "w");
  if (file == NULL) error("Could not create file '%s'.", file_name);
  fprintf(file, "# type\tname\n");
  for (int type = 0; type < task_type_count; type++) {
    fprintf(file, "%i\t%s\n", type, taskID_names[type]);
  }
  fclose(file);
  sprintf(file_name, "%s_task_subtypes.txt", file_prefix);
  file = fopen(file_name, "w");
  if (file == NULL) error("Could not create file '%s'.", file_name);
  fprintf(file, "# subtype\tname\n");
  for (int subtype = 0; subtype < task_subtype_count; subtype++) {
    fprintf(file, "%i\t%s\n", subtype, subtaskID_names[subtype]);
  }
  fclose(file);
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
/**
 * @brief Create global communicators for each of the subtasks.
 */
void task_free_mpi_comms(void) {
  for (int i = 0; i < task_subtype_count; i++) {
    MPI_Comm_free(&subtaskMPI_comms[i]);
  }
}
#endif

/**
 * @brief dump all the tasks of all the known engines into a file for
 * postprocessing.
 *
 * Dumps the information to a file "thread_info-stepn.dat" where n is the
 * given step value, or "thread_info_MPI-stepn.dat", if we are running
 * under MPI. Note if running under MPI all the ranks are dumped into this
 * one file, which has an additional field to identify the rank.
 *
 * @param e the #engine
 * @param step the current step.
 */
void task_dump_all(struct engine *e, int step) {

#ifdef SWIFT_DEBUG_TASKS

  const ticks tic = getticks();

  /* Need this to convert ticks to seconds. */
  const unsigned long long cpufreq = clocks_get_cpufreq();

#ifdef WITH_MPI
  /* Make sure output file is empty, only on one rank. */
  char dumpfile[35];
  snprintf(dumpfile, sizeof(dumpfile), "thread_info_MPI-step%d.dat", step);
  FILE *file_thread;
  if (engine_rank == 0) {
    file_thread = fopen(dumpfile, "w");
    if (file_thread == NULL)
      error("Could not create/erase file '%s'.", dumpfile);
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
      if (file_thread == NULL)
        error("Could not open file '%s' for writing.", dumpfile);

      /* Add some information to help with the plots and conversion of ticks to
       * seconds. */
      fprintf(file_thread, " %03d 0 0 0 0 %lld %lld %lld %lld %lld 0 0 %lld\n",
              engine_rank, (long long int)e->tic_step,
              (long long int)e->toc_step, e->updates, e->g_updates,
              e->s_updates, cpufreq);
      int count = 0;
      for (int l = 0; l < e->sched.nr_tasks; l++) {
        if (!e->sched.tasks[l].implicit &&
            e->sched.tasks[l].tic > e->tic_step) {
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
  if (file_thread == NULL) error("Could not create file '%s'.", dumpfile);

  /* Add some information to help with the plots and conversion of ticks to
   * seconds. */
  fprintf(file_thread, " %d %d %d %d %lld %lld %lld %lld %lld %d %lld\n", -2,
          -1, -1, 1, (unsigned long long)e->tic_step,
          (unsigned long long)e->toc_step, e->updates, e->g_updates,
          e->s_updates, 0, cpufreq);
  for (int l = 0; l < e->sched.nr_tasks; l++) {
    if (!e->sched.tasks[l].implicit && e->sched.tasks[l].tic > e->tic_step) {
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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
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
 * and total time taken and the same numbers for the start of the task,
 * in millisec and then the fixed costs value.
 *
 * If header is set, only the fixed costs value is written into the output
 * file in a format that is suitable for inclusion in SWIFT (as
 * partition_fixed_costs.h).
 *
 * @param dumpfile name of the file for the output.
 * @param e the #engine
 * @param dump_tasks_threshold Fraction of the step time above whic any task
 * triggers a call to task_dump_all().
 * @param header whether to write a header include file.
 * @param allranks do the statistics over all ranks, if not just the current
 *                 one, only used if header is false.
 */
void task_dump_stats(const char *dumpfile, struct engine *e,
                     float dump_tasks_threshold, int header, int allranks) {

  const ticks function_tic = getticks();

  /* Need arrays for sum, min and max across all types and subtypes. */
  double sum[task_type_count][task_subtype_count];
  double tsum[task_type_count][task_subtype_count];
  double min[task_type_count][task_subtype_count];
  double tmin[task_type_count][task_subtype_count];
  double max[task_type_count][task_subtype_count];
  double tmax[task_type_count][task_subtype_count];
  int count[task_type_count][task_subtype_count];

  for (int j = 0; j < task_type_count; j++) {
    for (int k = 0; k < task_subtype_count; k++) {
      sum[j][k] = 0.0;
      tsum[j][k] = 0.0;
      count[j][k] = 0;
      min[j][k] = DBL_MAX;
      tmin[j][k] = DBL_MAX;
      max[j][k] = 0.0;
      tmax[j][k] = 0.0;
    }
  }

  double stepdt = (double)e->toc_step - (double)e->tic_step;
  double total[1] = {0.0};
  int dumped_plot_data = 0;
  for (int l = 0; l < e->sched.nr_tasks; l++) {
    int type = e->sched.tasks[l].type;

    /* Skip implicit tasks and tasks that have not ran. */
    if (!e->sched.tasks[l].implicit && e->sched.tasks[l].tic > 0) {
      int subtype = e->sched.tasks[l].subtype;

      double dt = e->sched.tasks[l].toc - e->sched.tasks[l].tic;
      sum[type][subtype] += dt;

      double tic = (double)e->sched.tasks[l].tic;
      tsum[type][subtype] += tic;
      count[type][subtype] += 1;
      if (dt < min[type][subtype]) {
        min[type][subtype] = dt;
      }
      if (tic < tmin[type][subtype]) {
        tmin[type][subtype] = tic;
      }
      if (dt > max[type][subtype]) {
        max[type][subtype] = dt;
      }
      if (tic > tmax[type][subtype]) {
        tmax[type][subtype] = tic;
      }
      total[0] += dt;

      /* Check if this is a problematic task and make a report. */
      if (dump_tasks_threshold > 0. && dt / stepdt > dump_tasks_threshold) {

        if (e->verbose)
          message(
              "Long running task detected: %s/%s using %.1f%% of step runtime",
              taskID_names[type], subtaskID_names[subtype],
              dt / stepdt * 100.0);

        if (!dumped_plot_data) {
#ifdef SWIFT_DEBUG_TASKS
          task_dump_all(e, e->step + 1);
#endif
          dumped_plot_data = 1;
        }
      }
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

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : tsum), tsum, size,
                     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task tsums");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : count), count, size,
                     MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task counts");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : min), min, size,
                     MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task minima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : tmin), tmin, size,
                     MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task minima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : max), max, size,
                     MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task maxima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : tmax), tmax, size,
                     MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task maxima");

    res = MPI_Reduce((engine_rank == 0 ? MPI_IN_PLACE : total), total, 1,
                     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to reduce task total time");
  }

  if (!allranks || (engine_rank == 0 && (allranks || header))) {
#endif

    FILE *dfile = fopen(dumpfile, "w");
    if (dfile == NULL) error("Could not create file '%s'.", dumpfile);
    if (header) {
      fprintf(dfile, "/* use as src/partition_fixed_costs.h */\n");
      fprintf(dfile, "#define HAVE_FIXED_COSTS 1\n");
    } else {
      fprintf(dfile,
              "# task ntasks min max sum mean percent mintic maxtic"
              " meantic fixed_cost\n");
    }

    for (int j = 0; j < task_type_count; j++) {
      const char *taskID = taskID_names[j];
      for (int k = 0; k < task_subtype_count; k++) {
        if (sum[j][k] > 0.0) {

          /* Fixed cost is in .1ns as we want to compare between runs in
           * some absolute units. */
          double mean = sum[j][k] / (double)count[j][k];
          int fixed_cost = (int)(clocks_from_ticks(mean) * 10000.f);
          if (header) {
            fprintf(dfile, "repartition_costs[%d][%d] = %10d; /* %s/%s */\n", j,
                    k, fixed_cost, taskID, subtaskID_names[k]);
          } else {
            double perc = 100.0 * sum[j][k] / total[0];
            double mintic = tmin[j][k] - e->tic_step;
            double maxtic = tmax[j][k] - e->tic_step;
            double meantic = tsum[j][k] / (double)count[j][k] - e->tic_step;
            fprintf(dfile,
                    "%15s/%-10s %10d %14.4f %14.4f %14.4f %14.4f %14.4f"
                    " %14.4f %14.4f %14.4f %10d\n",
                    taskID, subtaskID_names[k], count[j][k],
                    clocks_from_ticks(min[j][k]), clocks_from_ticks(max[j][k]),
                    clocks_from_ticks(sum[j][k]), clocks_from_ticks(mean), perc,
                    clocks_from_ticks(mintic), clocks_from_ticks(maxtic),
                    clocks_from_ticks(meantic), fixed_cost);
          }
        }
      }
    }
    fclose(dfile);
#ifdef WITH_MPI
  }
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - function_tic),
            clocks_getunit());
}

/**
 * @brief dump all the active tasks of all the known engines into files.
 *
 * Dumps the information into file "task_dump-stepn.dat" where n is the given
 * step value, or files "task_dump_MPI-stepn.dat_rank", if we are running
 * under MPI. Note if running under MPI all the ranks are dumped into separate
 * files to avoid interaction with other MPI calls that may be blocking at the
 * time. Very similar to task_dump_all() except for the additional fields used
 * in task debugging and we record tasks that have not ran (i.e !skip, but toc
 * == 0) and how many waits are still active.
 *
 * @param e the #engine
 */
void task_dump_active(struct engine *e) {

  const ticks tic = getticks();

  /* Need this to convert ticks to seconds. */
  unsigned long long cpufreq = clocks_get_cpufreq();
  char dumpfile[35];

#ifdef WITH_MPI
  snprintf(dumpfile, sizeof(dumpfile), "task_dump_MPI-step%d.dat_%d", e->step,
           e->nodeID);
#else
  snprintf(dumpfile, sizeof(dumpfile), "task_dump-step%d.dat", e->step);
#endif

  FILE *file_thread = fopen(dumpfile, "w");
  if (file_thread == NULL) error("Could not create file '%s'.", dumpfile);
  fprintf(file_thread,
          "# rank otherrank type subtype waits pair tic toc"
          " ci.hydro.count cj.hydro.count ci.grav.count cj.grav.count"
          " flags\n");

  /* Add some information to help with the plots and conversion of ticks to
   * seconds. */
  fprintf(file_thread, "%i 0 none none -1 0 %lld %lld %lld %lld %lld 0 %lld\n",
          engine_rank, (long long int)e->tic_step, (long long int)e->toc_step,
          e->updates, e->g_updates, e->s_updates, cpufreq);
  for (int l = 0; l < e->sched.nr_tasks; l++) {
    struct task *t = &e->sched.tasks[l];

    /* Not implicit and not skipped. */
    if (!t->implicit && !t->skip) {

      /* Get destination rank of MPI requests. */
      int paired = (t->cj != NULL);
      int otherrank = t->ci->nodeID;
      if (paired) otherrank = t->cj->nodeID;

      fprintf(file_thread, "%i %i %s %s %i %i %lli %lli %i %i %i %i %lli\n",
              engine_rank, otherrank, taskID_names[t->type],
              subtaskID_names[t->subtype], t->wait, paired,
              (long long int)t->tic, (long long int)t->toc,
              (t->ci != NULL) ? t->ci->hydro.count : 0,
              (t->cj != NULL) ? t->cj->hydro.count : 0,
              (t->ci != NULL) ? t->ci->grav.count : 0,
              (t->cj != NULL) ? t->cj->grav.count : 0, t->flags);
    }
  }
  fclose(file_thread);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Return the #task_categories of a given #task.
 *
 * @param t The #task.
 */
enum task_categories task_get_category(const struct task *t) {

  switch (t->type) {

    case task_type_cooling:
      return task_category_cooling;

    case task_type_csds:
      return task_category_csds;

    case task_type_star_formation:
    case task_type_star_formation_sink:
      return task_category_star_formation;

    case task_type_sink_formation:
      return task_category_sink;

    case task_type_drift_part:
    case task_type_drift_spart:
    case task_type_drift_sink:
    case task_type_drift_bpart:
    case task_type_drift_gpart:
      return task_category_drift;

    case task_type_sort:
    case task_type_stars_sort:
      return task_category_sort;

    case task_type_stars_resort:
      return task_category_resort;

    case task_type_send:
    case task_type_recv:
      return task_category_mpi;

    case task_type_pack:
    case task_type_unpack:
      return task_category_pack;

    case task_type_kick1:
    case task_type_kick2:
    case task_type_timestep:
    case task_type_collect:
      return task_category_time_integration;

    case task_type_timestep_limiter:
      return task_category_limiter;

    case task_type_timestep_sync:
      return task_category_sync;

    case task_type_ghost:
      return task_category_hydro_ghost;

    case task_type_extra_ghost:
    case task_type_end_hydro_force:
      return task_category_hydro;

    case task_type_stars_ghost:
    case task_type_stars_prep_ghost1:
    case task_type_hydro_prep_ghost1:
    case task_type_stars_prep_ghost2:
      return task_category_feedback;

    case task_type_bh_density_ghost:
    case task_type_bh_swallow_ghost2:
      return task_category_black_holes;

    case task_type_init_grav:
    case task_type_grav_long_range:
    case task_type_grav_mm:
    case task_type_grav_down:
    case task_type_end_grav_force:
      return task_category_gravity;

    case task_type_fof_self:
    case task_type_fof_pair:
    case task_type_fof_attach_self:
    case task_type_fof_attach_pair:
      return task_category_fof;

    case task_type_rt_in:
    case task_type_rt_ghost1:
    case task_type_rt_ghost2:
    case task_type_rt_transport_out:
    case task_type_rt_out:
    case task_type_rt_sort:
    case task_type_rt_advance_cell_time:
      return task_category_rt;

    case task_type_rt_tchem:
      return task_category_rt_tchem;

    case task_type_neutrino_weight:
      return task_category_neutrino;

    case task_type_self:
    case task_type_pair:
    case task_type_sub_self:
    case task_type_sub_pair: {
      switch (t->subtype) {

        case task_subtype_density:
          return task_category_hydro_density;

        case task_subtype_gradient:
        case task_subtype_force:
          return task_category_hydro;

        case task_subtype_limiter:
          return task_category_limiter;

        case task_subtype_grav:
        case task_subtype_external_grav:
          return task_category_gravity;

        case task_subtype_stars_density:
        case task_subtype_stars_prep1:
        case task_subtype_stars_prep2:
        case task_subtype_stars_feedback:
          return task_category_feedback;

        case task_subtype_bh_density:
        case task_subtype_bh_swallow:
        case task_subtype_do_gas_swallow:
        case task_subtype_do_bh_swallow:
        case task_subtype_bh_feedback:
          return task_category_black_holes;

        case task_subtype_sink_swallow:
        case task_subtype_sink_do_sink_swallow:
        case task_subtype_sink_do_gas_swallow:
          return task_category_sink;

        case task_subtype_rt_gradient:
        case task_subtype_rt_transport:
          return task_category_rt;

        default:
          return task_category_others;
      }
    }

    default:
      return task_category_others;
  }
}
