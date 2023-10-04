/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SCHEDULER_H
#define SWIFT_SCHEDULER_H

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Some standard headers. */
#include <pthread.h>

/* Includes. */
#include "cell.h"
#include "inline.h"
#include "lock.h"
#include "queue.h"
#include "task.h"
#include "threadpool.h"

/* Some constants. */
#define scheduler_maxwait 3
#define scheduler_init_nr_unlocks 10000
#define scheduler_dosub 1
#define scheduler_maxsteal 10
#define scheduler_maxtries 2
#define scheduler_doforcesplit            \
  0 /* Beware: switching this on can/will \
       break engine_addlink as it assumes \
       a maximum number of tasks per cell. */

/* Flags . */
#define scheduler_flag_none 0
#define scheduler_flag_steal (1 << 1)

/* Data of a scheduler. */
struct scheduler {
  /* Scheduler flags. */
  unsigned int flags;

  /* Number of queues in this scheduler. */
  int nr_queues;

  /* Array of queues. */
  struct queue *queues;

  /* Total number of tasks. */
  int nr_tasks, size, tasks_next;

  /* Total number of waiting tasks. */
  int waiting;

  /* The task array. */
  struct task *tasks;

  /* The task indices. */
  int *tasks_ind;

  /* List of initial tasks. */
  int *tid_active;
  int active_count;

  /* The task unlocks. */
  struct task **volatile unlocks;
  int *volatile unlock_ind;
  volatile int nr_unlocks, size_unlocks, completed_unlock_writes;

  /* Lock for this scheduler. */
  swift_lock_type lock;

  /* Waiting queue. */
  pthread_mutex_t sleep_mutex;
  pthread_cond_t sleep_cond;

  /* The space associated with this scheduler. */
  struct space *space;

  /* Threadpool to use internally for mundane parallel work. */
  struct threadpool *threadpool;

  /* The node we are working on. */
  int nodeID;

  /* Maximum size of task messages, in bytes, to sent using non-buffered
   * MPI. */
  size_t mpi_message_limit;

  /* Total ticks spent running the tasks */
  ticks total_ticks;

  struct {
    /* Total ticks spent waiting for runners to come home. */
    ticks waiting_ticks;

    /* Total ticks spent by runners running tasks. */
    ticks active_ticks;
  } deadtime;

  /* Frequency of the dependency graph dumping. */
  int frequency_dependency;

  /* Specific cell to dump dependency graph for */
  long long dependency_graph_cellID;

  /* Frequency of the task levels dumping. */
  int frequency_task_levels;

#if defined(SWIFT_DEBUG_CHECKS)
  /* Stuff for the deadlock detector */

  /* How long to wait (in ms) before assuming we're in a deadlock */
  float deadlock_waiting_time_ms;

  /* Time at which last task was successfully retrieved from a queue */
  ticks last_successful_task_fetch;

  /* needed to dump queues on deadlock detection */
  struct engine *e;

#endif /* SWIFT_DEBUG_CHECKS */
};

/* Inlined functions (for speed). */
/**
 * @brief Add a regular task to the list of active tasks.
 *
 * @param s The #scheduler.
 * @param t The task to be added.
 */
__attribute__((always_inline)) INLINE static void scheduler_activate(
    struct scheduler *s, struct task *t) {
  if (atomic_cas(&t->skip, 1, 0)) {
    t->wait = 0;
    int ind = atomic_inc(&s->active_count);
    s->tid_active[ind] = t - s->tasks;
  }
}

/**
 * @brief Search a given linked list of task for a given subtype and activate
 * it.
 *
 * @param s The #scheduler.
 * @param link The first element in the linked list of links for the task of
 * interest.
 * @param subtype the task subtype to activate.
 */
__attribute__((always_inline)) INLINE static void
scheduler_activate_all_subtype(struct scheduler *s, struct link *link,
                               const enum task_subtypes subtype) {

  for (struct link *l = link; l != NULL; l = l->next) {
    if (l->t->subtype == subtype) scheduler_activate(s, l->t);
  }
}

/**
 * @brief Search and add an MPI send task to the list of active tasks.
 *
 * @param s The #scheduler.
 * @param link The first element in the linked list of links for the task of
 * interest.
 * @param subtype the task subtype to activate.
 * @param nodeID The nodeID of the foreign cell.
 *
 * @return The #link to the MPI send task.
 */
__attribute__((always_inline)) INLINE static struct link *
scheduler_activate_send(struct scheduler *s, struct link *link,
                        const enum task_subtypes subtype, const int nodeID) {
  struct link *l = NULL;
  for (l = link;
       l != NULL && !(l->t->cj->nodeID == nodeID && l->t->subtype == subtype);
       l = l->next)
    ;
  if (l == NULL) {
    error("Missing link to send task.");
  }
  scheduler_activate(s, l->t);
  return l;
}

/**
 * @brief Search and add an MPI recv task to the list of active tasks.
 *
 * @param s The #scheduler.
 * @param link The first element in the linked list of links for the task of
 * interest.
 * @param subtype the task subtype to activate.
 *
 * @return The #link to the MPI recv task.
 */
__attribute__((always_inline)) INLINE static struct link *
scheduler_activate_recv(struct scheduler *s, struct link *link,
                        const enum task_subtypes subtype) {
  struct link *l = NULL;
  for (l = link; l != NULL && l->t->subtype != subtype; l = l->next)
    ;
  if (l == NULL) {
    error("Missing link to recv task.");
  }
  scheduler_activate(s, l->t);
  return l;
}

/**
 * @brief Search and add an MPI pack task to the list of active tasks.
 *
 * @param s The #scheduler.
 * @param link The first element in the linked list of links for the task of
 * interest.
 * @param subtype the task subtype to activate.
 * @param nodeID The nodeID of the foreign cell.
 *
 * @return The #link to the MPI pack task.
 */
__attribute__((always_inline)) INLINE static struct link *
scheduler_activate_pack(struct scheduler *s, struct link *link,
                        enum task_subtypes subtype, int nodeID) {
  struct link *l = NULL;
  for (l = link;
       l != NULL && !(l->t->cj->nodeID == nodeID && l->t->subtype == subtype);
       l = l->next)
    ;
  if (l == NULL) {
    error("Missing link to pack task.");
  }
  scheduler_activate(s, l->t);
  return l;
}

/**
 * @brief Search and add an MPI unpack task to the list of active tasks.
 *
 * @param s The #scheduler.
 * @param link The first element in the linked list of links for the task of
 * interest.
 * @param subtype the task subtype to activate.
 *
 * @return The #link to the MPI unpack task.
 */
__attribute__((always_inline)) INLINE static struct link *
scheduler_activate_unpack(struct scheduler *s, struct link *link,
                          enum task_subtypes subtype) {
  struct link *l = NULL;
  for (l = link; l != NULL && l->t->subtype != subtype; l = l->next)
    ;
  if (l == NULL) {
    error("Missing link to unpack task.");
  }
  scheduler_activate(s, l->t);
  return l;
}

/* Function prototypes. */
void scheduler_clear_active(struct scheduler *s);
void scheduler_init(struct scheduler *s, struct space *space, int nr_tasks,
                    int nr_queues, unsigned int flags, int nodeID,
                    struct threadpool *tp);
struct task *scheduler_gettask(struct scheduler *s, int qid,
                               const struct task *prev);
void scheduler_enqueue(struct scheduler *s, struct task *t);
void scheduler_start(struct scheduler *s);
void scheduler_reset(struct scheduler *s, int nr_tasks);
void scheduler_ranktasks(struct scheduler *s);
void scheduler_reweight(struct scheduler *s, int verbose);
struct task *scheduler_addtask(struct scheduler *s, enum task_types type,
                               enum task_subtypes subtype, long long flags,
                               int implicit, struct cell *ci, struct cell *cj);
void scheduler_splittasks(struct scheduler *s, const int fof_tasks,
                          const int verbose);
struct task *scheduler_done(struct scheduler *s, struct task *t);
struct task *scheduler_unlock(struct scheduler *s, struct task *t);
void scheduler_addunlock(struct scheduler *s, struct task *ta, struct task *tb);
void scheduler_set_unlocks(struct scheduler *s);
void scheduler_dump_queue(struct scheduler *s);
void scheduler_print_tasks(const struct scheduler *s, const char *fileName);
void scheduler_clean(struct scheduler *s);
void scheduler_free_tasks(struct scheduler *s);
void scheduler_write_dependencies(struct scheduler *s, int verbose, int step);
void scheduler_write_cell_dependencies(struct scheduler *s, int verbose,
                                       int step);
void scheduler_write_task_level(const struct scheduler *s, int step);
void scheduler_dump_queues(struct engine *e);
void scheduler_report_task_times(const struct scheduler *s,
                                 const int nr_threads);

#endif /* SWIFT_SCHEDULER_H */
