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
#ifndef SWIFT_ENGINE_H
#define SWIFT_ENGINE_H

/* Config parameters. */
#include "../config.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Some standard headers. */
#include <pthread.h>
#include <stdio.h>

/* Includes. */
#include "clocks.h"
#include "collectgroup.h"
#include "cooling_struct.h"
#include "gravity_properties.h"
#include "parser.h"
#include "partition.h"
#include "potential.h"
#include "runner.h"
#include "scheduler.h"
#include "sourceterms_struct.h"
#include "space.h"
#include "task.h"
#include "units.h"

/* Some constants. */
enum engine_policy {
  engine_policy_none = 0,
  engine_policy_rand = (1 << 0),
  engine_policy_steal = (1 << 1),
  engine_policy_keep = (1 << 2),
  engine_policy_block = (1 << 3),
  engine_policy_cputight = (1 << 4),
  engine_policy_mpi = (1 << 5),
  engine_policy_setaffinity = (1 << 6),
  engine_policy_hydro = (1 << 7),
  engine_policy_self_gravity = (1 << 8),
  engine_policy_external_gravity = (1 << 9),
  engine_policy_cosmology = (1 << 10),
  engine_policy_drift_all = (1 << 11),
  engine_policy_reconstruct_mpoles = (1 << 12),
  engine_policy_cooling = (1 << 13),
  engine_policy_sourceterms = (1 << 14),
  engine_policy_stars = (1 << 15)
};
#define engine_maxpolicy 15
extern const char *engine_policy_names[];

#define engine_queue_scale 1.2
#define engine_maxproxies 64
#define engine_tasksreweight 1
#define engine_parts_size_grow 1.05
#define engine_redistribute_alloc_margin 1.2
#define engine_default_energy_file_name "energy"
#define engine_default_timesteps_file_name "timesteps"
#define engine_max_parts_per_ghost 1000

/* The rank of the engine as a global variable (for messages). */
extern int engine_rank;

/* Data structure for the engine. */
struct engine {

  /* Number of threads on which to run. */
  int nr_threads;

  /* The space with which the runner is associated. */
  struct space *s;

  /* The runner's threads. */
  struct runner *runners;

  /* The running policy. */
  int policy;

  /* The task scheduler. */
  struct scheduler sched;

  /* Common threadpool for all the engine's tasks. */
  struct threadpool threadpool;

  /* The minimum and maximum allowed dt */
  double dt_min, dt_max;

  /* Time of the simulation beginning */
  double timeBegin;

  /* Time of the simulation end */
  double timeEnd;

  /* The previous system time. */
  double timeOld;
  integertime_t ti_old;

  /* The current system time. */
  double time;
  integertime_t ti_current;

  /* The highest active bin at this time */
  timebin_t max_active_bin;

  /* Time step */
  double timeStep;

  /* Time base */
  double timeBase;
  double timeBase_inv;

  /* Minimal ti_end for the next time-step */
  integertime_t ti_end_min;

  /* Maximal ti_end for the next time-step */
  integertime_t ti_end_max;

  /* Maximal ti_beg for the next time-step */
  integertime_t ti_beg_max;

  /* Number of particles updated */
  size_t updates, g_updates, s_updates;

  /* Total numbers of particles in the system. */
  size_t total_nr_parts, total_nr_gparts;

  /* The internal system of units */
  const struct unit_system *internal_units;

  /* Snapshot information */
  double timeFirstSnapshot;
  double deltaTimeSnapshot;
  integertime_t ti_nextSnapshot;
  char snapshotBaseName[PARSER_MAX_LINE_SIZE];
  int snapshotCompression;
  struct unit_system *snapshotUnits;

  /* Statistics information */
  FILE *file_stats;
  double timeLastStatistics;
  double deltaTimeStatistics;

  /* Timesteps information */
  FILE *file_timesteps;

  /* The current step number. */
  int step;

  /* The number of particles updated in the previous step. */
  int count_step;

  /* Data for the threads' barrier. */
  pthread_barrier_t wait_barrier;
  pthread_barrier_t run_barrier;

  /* ID of the node this engine lives on. */
  int nr_nodes, nodeID;

  /* Proxies for the other nodes in this simulation. */
  struct proxy *proxies;
  int nr_proxies, *proxy_ind;

#ifdef SWIFT_DEBUG_TASKS
  /* Tic/toc at the start/end of a step. */
  ticks tic_step, toc_step;
#endif

#ifdef WITH_MPI
  /* CPU time of the last step. */
  double cputime_last_step;

  /* Step of last repartition. */
  int last_repartition;
#endif

  /* Wallclock time of the last time-step */
  float wallclock_time;

  /* Force the engine to rebuild? */
  int forcerebuild;

  /* Force the engine to repartition ? */
  int forcerepart;
  struct repartition *reparttype;

  /* Need to dump some statistics ? */
  int save_stats;

  /* Need to dump a snapshot ? */
  int dump_snapshot;

  /* How many steps have we done with the same set of tasks? */
  int tasks_age;

  /* Linked list for cell-task association. */
  struct link *links;
  int nr_links, size_links;

  /* Average number of tasks per cell. Used to estimate the sizes
   * of the various task arrays. */
  int tasks_per_cell;

  /* Are we talkative ? */
  int verbose;

  /* Physical constants definition */
  const struct phys_const *physical_constants;

  /* Properties of the hydro scheme */
  const struct hydro_props *hydro_properties;

  /* Properties of the self-gravity scheme */
  const struct gravity_props *gravity_properties;

  /* Properties of external gravitational potential */
  const struct external_potential *external_potential;

  /* Properties of the cooling scheme */
  const struct cooling_function_data *cooling_func;

  /* Properties of source terms */
  struct sourceterms *sourceterms;

  /* The (parsed) parameter file */
  const struct swift_params *parameter_file;

  /* Temporary struct to hold a group of deferable properties (in MPI mode
   * these are reduced together, but may not be required just yet). */
  struct collectgroup1 collect_group1;
};

/* Function prototypes. */
void engine_barrier(struct engine *e);
void engine_compute_next_snapshot_time(struct engine *e);
void engine_unskip(struct engine *e);
void engine_drift_all(struct engine *e);
void engine_drift_top_multipoles(struct engine *e);
void engine_reconstruct_multipoles(struct engine *e);
void engine_print_stats(struct engine *e);
void engine_dump_snapshot(struct engine *e);
void engine_init(struct engine *e, struct space *s,
                 const struct swift_params *params, int nr_nodes, int nodeID,
                 int nr_threads, int Ngas, int Ndm, int with_aff, int policy,
                 int verbose, struct repartition *reparttype,
                 const struct unit_system *internal_units,
                 const struct phys_const *physical_constants,
                 const struct hydro_props *hydro,
                 const struct gravity_props *gravity,
                 const struct external_potential *potential,
                 const struct cooling_function_data *cooling_func,
                 struct sourceterms *sourceterms);
void engine_launch(struct engine *e);
void engine_prepare(struct engine *e);
void engine_init_particles(struct engine *e, int flag_entropy_ICs,
                           int clean_h_values);
void engine_step(struct engine *e);
void engine_maketasks(struct engine *e);
void engine_split(struct engine *e, struct partition *initial_partition);
void engine_exchange_strays(struct engine *e, size_t offset_parts,
                            int *ind_part, size_t *Npart, size_t offset_gparts,
                            int *ind_gpart, size_t *Ngpart,
                            size_t offset_sparts, int *ind_spart,
                            size_t *Nspart);
void engine_rebuild(struct engine *e, int clean_h_values);
void engine_repartition(struct engine *e);
void engine_repartition_trigger(struct engine *e);
void engine_makeproxies(struct engine *e);
void engine_redistribute(struct engine *e);
void engine_print_policy(struct engine *e);
int engine_is_done(struct engine *e);
void engine_pin();
void engine_unpin();
void engine_clean(struct engine *e);
int engine_estimate_nr_tasks(struct engine *e);

#endif /* SWIFT_ENGINE_H */
