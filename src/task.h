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
#ifndef SWIFT_TASK_H
#define SWIFT_TASK_H

#include "../config.h"

/* Includes. */
#include "align.h"
#include "cell.h"
#include "cycle.h"

#define task_align 128

/**
 * @brief The different task types.
 *
 * Be sure to update the taskID_names array in tasks.c if you modify this list!
 * Also update the python task plotting scripts!
 */
enum task_types {
  task_type_none = 0,
  task_type_sort,
  task_type_self,
  task_type_pair,
  task_type_sub_self,
  task_type_sub_pair,
  task_type_init_grav,
  task_type_init_grav_out, /* Implicit */
  task_type_ghost_in,      /* Implicit */
  task_type_ghost,
  task_type_ghost_out, /* Implicit */
  task_type_extra_ghost,
  task_type_drift_part,
  task_type_drift_gpart,
  task_type_drift_gpart_out, /* Implicit */
  task_type_end_force,
  task_type_kick1,
  task_type_kick2,
  task_type_timestep,
  task_type_timestep_limiter,
  task_type_send,
  task_type_recv,
  task_type_grav_long_range,
  task_type_grav_mm,
  task_type_grav_down_in, /* Implicit */
  task_type_grav_down,
  task_type_grav_mesh,
  task_type_cooling,
  task_type_star_formation,
  task_type_logger,
  task_type_stars_ghost_in,
  task_type_stars_ghost,
  task_type_stars_ghost_out,
  task_type_stars_sort,
  task_type_fof_self,
  task_type_fof_pair,
  task_type_count
} __attribute__((packed));

/**
 * @brief The different task sub-types (for pairs, selfs and sub-tasks).
 */
enum task_subtypes {
  task_subtype_none = 0,
  task_subtype_density,
  task_subtype_gradient,
  task_subtype_force,
  task_subtype_limiter,
  task_subtype_grav,
  task_subtype_external_grav,
  task_subtype_tend,
  task_subtype_xv,
  task_subtype_rho,
  task_subtype_gpart,
  task_subtype_multipole,
  task_subtype_spart,
  task_subtype_stars_density,
  task_subtype_stars_feedback,
  task_subtype_count
} __attribute__((packed));

/**
 * @brief The type of particles/objects this task acts upon in a given cell.
 */
enum task_actions {
  task_action_none,
  task_action_part,
  task_action_gpart,
  task_action_spart,
  task_action_all,
  task_action_multipole,
  task_action_count
};

/**
 * @brief Names of the task types.
 */
extern const char *taskID_names[];

/**
 * @brief Names of the task sub-types.
 */
extern const char *subtaskID_names[];

/**
 *  @brief The MPI communicators for the different subtypes.
 */
#ifdef WITH_MPI
extern MPI_Comm subtaskMPI_comms[task_subtype_count];
#endif

/**
 * @brief A task to be run by the #scheduler.
 */
struct task {

  /*! Pointers to the cells this task acts upon */
  struct cell *ci, *cj;

  /*! List of tasks unlocked by this one */
  struct task **unlock_tasks;

  /*! Flags used to carry additional information (e.g. sort directions) */
  long long flags;

#ifdef WITH_MPI

  /*! Buffer for this task's communications */
  void *buff;

  /*! MPI request corresponding to this task */
  MPI_Request req;

#endif

  /*! Rank of a task in the order */
  int rank;

  /*! Weight of the task */
  float weight;

  /*! Number of tasks unlocked by this one */
  short int nr_unlock_tasks;

  /*! Number of unsatisfied dependencies */
  short int wait;

  /*! Type of the task */
  enum task_types type;

  /*! Sub-type of the task (for the tasks that have one */
  enum task_subtypes subtype;

  /*! Should the scheduler skip this task ? */
  char skip;

  /*! Is this task implicit (i.e. does not do anything) ? */
  char implicit;

#ifdef SWIFT_DEBUG_TASKS
  /*! ID of the queue or runner owning this task */
  short int rid;

  /*! Information about the direction of the pair task */
  short int sid;
#endif

  /*! Start and end time of this task */
  ticks tic, toc;

#ifdef SWIFT_DEBUG_CHECKS
  /* When was this task last run? */
  integertime_t ti_run;
#endif /* SWIFT_DEBUG_CHECKS */

} SWIFT_STRUCT_ALIGN;

/* Function prototypes. */
void task_unlock(struct task *t);
float task_overlap(const struct task *ta, const struct task *tb);
int task_lock(struct task *t);
void task_do_rewait(struct task *t);
void task_print(const struct task *t);
void task_dump_all(struct engine *e, int step);
void task_dump_stats(const char *dumpfile, struct engine *e, int header,
                     int allranks);
void task_get_full_name(int type, int subtype, char *name);
void task_get_group_name(int type, int subtype, char *cluster);

#ifdef WITH_MPI
void task_create_mpi_comms(void);
#endif
#endif /* SWIFT_TASK_H */
