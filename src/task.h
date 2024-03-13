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
#ifndef SWIFT_TASK_H
#define SWIFT_TASK_H

#include <config.h>

/* Includes. */
#include "align.h"
#include "cycle.h"
#include "timeline.h"

/* Forward declarations to avoid circular inclusion dependencies. */
struct cell;
struct engine;

#define task_align 128

/**
 * @brief The different task types.
 *
 * Be sure to update the taskID_names array in tasks.c if you modify this list!
 * Also update the python3 task plotting scripts!
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
  task_type_drift_spart,
  task_type_drift_sink,
  task_type_drift_bpart,
  task_type_drift_gpart,
  task_type_drift_gpart_out, /* Implicit */
  task_type_end_hydro_force,
  task_type_kick1,
  task_type_kick2,
  task_type_timestep,
  task_type_timestep_limiter,
  task_type_timestep_sync,
  task_type_collect,
  task_type_send,
  task_type_recv,
  task_type_pack,
  task_type_unpack,
  task_type_grav_long_range,
  task_type_grav_mm,
  task_type_grav_down_in, /* Implicit */
  task_type_grav_down,
  task_type_end_grav_force,
  task_type_cooling,
  task_type_cooling_in,  /* Implicit */
  task_type_cooling_out, /* Implicit */
  task_type_star_formation,
  task_type_star_formation_in,  /* Implicit */
  task_type_star_formation_out, /* Implicit */
  task_type_star_formation_sink,
  task_type_csds,
  task_type_stars_in,       /* Implicit */
  task_type_stars_out,      /* Implicit */
  task_type_stars_ghost_in, /* Implicit */
  task_type_stars_ghost,
  task_type_stars_ghost_out,   /* Implicit */
  task_type_stars_prep_ghost1, /* Implicit */
  task_type_hydro_prep_ghost1, /* Implicit */
  task_type_stars_prep_ghost2, /* Implicit */
  task_type_stars_sort,
  task_type_stars_resort,
  task_type_bh_in,  /* Implicit */
  task_type_bh_out, /* Implicit */
  task_type_bh_density_ghost,
  task_type_bh_swallow_ghost1, /* Implicit */
  task_type_bh_swallow_ghost2,
  task_type_bh_swallow_ghost3, /* Implicit */
  task_type_fof_self,
  task_type_fof_pair,
  task_type_fof_attach_self,
  task_type_fof_attach_pair,
  task_type_neutrino_weight,
  task_type_sink_in,     /* Implicit */
  task_type_sink_ghost1, /* Implicit */
  task_type_sink_ghost2, /* Implicit */
  task_type_sink_out,    /* Implicit */
  task_type_rt_in,       /* Implicit */
  task_type_rt_out,      /* Implicit */
  task_type_sink_formation,
  task_type_rt_ghost1,
  task_type_rt_ghost2,
  task_type_rt_transport_out, /* Implicit */
  task_type_rt_tchem,
  task_type_rt_advance_cell_time,
  task_type_rt_sort,
  task_type_rt_collect_times,
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
  task_subtype_part_swallow,
  task_subtype_bpart_merger,
  task_subtype_gpart,
  task_subtype_spart_density,
  task_subtype_part_prep1,
  task_subtype_spart_prep2,
  task_subtype_stars_density,
  task_subtype_stars_prep1,
  task_subtype_stars_prep2,
  task_subtype_stars_feedback,
  task_subtype_sf_counts,
  task_subtype_bpart_rho,
  task_subtype_bpart_feedback,
  task_subtype_bh_density,
  task_subtype_bh_swallow,
  task_subtype_do_gas_swallow,
  task_subtype_do_bh_swallow,
  task_subtype_bh_feedback,
  task_subtype_sink_do_sink_swallow,
  task_subtype_sink_swallow,
  task_subtype_sink_do_gas_swallow,
  task_subtype_rt_gradient,
  task_subtype_rt_transport,
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
  task_action_sink,
  task_action_bpart,
  task_action_all,
  task_action_multipole,
  task_action_count
};

/**
 * @brief The broad categories of the tasks.
 */
enum task_categories {
  task_category_drift,
  task_category_sort,
  task_category_resort,
  task_category_hydro,
  task_category_gravity,
  task_category_feedback,
  task_category_black_holes,
  task_category_cooling,
  task_category_star_formation,
  task_category_limiter,
  task_category_sync,
  task_category_time_integration,
  task_category_mpi,
  task_category_pack,
  task_category_fof,
  task_category_others,
  task_category_neutrino,
  task_category_sink,
  task_category_rt,
  task_category_csds,
  task_category_count
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
 * @brief Names of the task categories.
 */
extern const char *task_category_names[];

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
  int nr_unlock_tasks;

  /*! Number of unsatisfied dependencies */
  int wait;

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

  /* Total time spent running this task */
  ticks total_ticks;

#ifdef SWIFT_DEBUG_CHECKS
  /* When was this task last run? */
  integertime_t ti_run;
#endif /* SWIFT_DEBUG_CHECKS */

} SWIFT_STRUCT_ALIGN;

/* Function prototypes. */
void task_unlock(struct task *t);
float task_overlap(const struct task *ta, const struct task *tb);
int task_lock(struct task *t);
struct task *task_get_unique_dependent(const struct task *t);
void task_print(const struct task *t);
void task_dump_all(struct engine *e, int step);
void task_dump_stats(const char *dumpfile, struct engine *e,
                     float dump_tasks_threshold, int header, int allranks);
void task_dump_active(struct engine *e);
void task_get_full_name(int type, int subtype, char *name);
void task_create_name_files(const char *file_prefix);
void task_get_group_name(int type, int subtype, char *cluster);
enum task_categories task_get_category(const struct task *t);

#ifdef WITH_MPI
void task_create_mpi_comms(void);
void task_free_mpi_comms(void);
#endif
#endif /* SWIFT_TASK_H */
