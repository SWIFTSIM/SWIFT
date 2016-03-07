/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "lock.h"
#include "proxy.h"
#include "runner.h"
#include "scheduler.h"
#include "space.h"
#include "task.h"
#include "partition.h"

/* Some constants. */
enum engine_policy {
  engine_policy_none = 0,
  engine_policy_rand = (1 << 0),
  engine_policy_steal = (1 << 1),
  engine_policy_keep = (1 << 2),
  engine_policy_block = (1 << 3),
  engine_policy_fixdt = (1 << 4),
  engine_policy_cputight = (1 << 5),
  engine_policy_mpi = (1 << 6),
  engine_policy_setaffinity = (1 << 7),
  engine_policy_hydro = (1 << 8),
  engine_policy_self_gravity = (1 << 9),
  engine_policy_external_gravity = (1 << 10)
};

extern const char *engine_policy_names[];

#define engine_queue_scale 1.2
#define engine_maxtaskspercell 96
#define engine_maxproxies 64
#define engine_tasksreweight 10

/* The rank of the engine as a global variable (for messages). */
extern int engine_rank;

/* The maximal number of timesteps in a simulation */
#define max_nr_timesteps (1 << 28)

/* Mini struct to link cells to density/force tasks. */
struct link {

  /* The task pointer. */
  struct task *t;

  /* The next pointer. */
  struct link *next;
};

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

  /* The minimum and maximum allowed dt */
  double dt_min, dt_max;

  /* Time of the simulation beginning */
  double timeBegin;

  /* Time of the simulation end */
  double timeEnd;

  /* The previous system time. */
  double timeOld;
  int ti_old;

  /* The current system time. */
  double time;
  int ti_current;

  /* Time step */
  double timeStep;

  /* Time base */
  double timeBase;

  /* File for statistics */
  FILE *file_stats;

  /* The current step number. */
  int step, nullstep;

  /* The number of particles updated in the previous step. */
  int count_step;

  /* Data for the threads' barrier. */
  pthread_mutex_t barrier_mutex;
  pthread_cond_t barrier_cond;
  volatile int barrier_running, barrier_launch, barrier_launchcount;

  /* ID of the node this engine lives on. */
  int nr_nodes, nodeID;

  /* Proxies for the other nodes in this simulation. */
  struct proxy *proxies;
  int nr_proxies, *proxy_ind;

  /* Tic at the start of a step. */
  ticks tic_step;

  /* Wallclock time of the last time-step */
  float wallclock_time;

  /* Force the engine to rebuild? */
  int forcerebuild;
  enum repartition_type forcerepart;

  /* How many steps have we done with the same set of tasks? */
  int tasks_age;

  /* Linked list for cell-task association. */
  struct link *links;
  int nr_links, size_links;

#ifdef WITH_MPI
  /* MPI data type for the particle transfers */
  MPI_Datatype part_mpi_type;
  MPI_Datatype xpart_mpi_type;
#endif
};

/* Function prototypes. */
void engine_barrier(struct engine *e, int tid);
void engine_init(struct engine *e, struct space *s, float dt, int nr_threads,
                 int nr_queues, int nr_nodes, int nodeID, int policy,
                 float timeBegin, float timeEnd, float dt_min, float dt_max);
void engine_launch(struct engine *e, int nr_runners, unsigned int mask,
                   unsigned int submask);
void engine_prepare(struct engine *e);
void engine_print(struct engine *e);
void engine_init_particles(struct engine *e);
void engine_step(struct engine *e);
void engine_maketasks(struct engine *e);
void engine_split(struct engine *e, struct partition *initial_partition);
int engine_exchange_strays(struct engine *e, int offset, int *ind, int N);
void engine_rebuild(struct engine *e);
void engine_repartition(struct engine *e);
void engine_makeproxies(struct engine *e);
void engine_redistribute(struct engine *e);
struct link *engine_addlink(struct engine *e, struct link *l, struct task *t);
void engine_print_policy(struct engine *e);
int engine_is_done(struct engine *e);

#endif /* SWIFT_ENGINE_H */
