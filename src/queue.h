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
#ifndef SWIFT_QUEUE_H
#define SWIFT_QUEUE_H

/* Includes. */
#include "cell.h"
#include "lock.h"
#include "task.h"

/* Some constants. */
#define queue_sizeinit 100
#define queue_sizegrow 2
#define queue_search_window 8
#define queue_incoming_size 10240
#define queue_struct_align 64

/* Constants dealing with task de-priorization. */
#define queue_lock_fail_reweight_factor 0.5

// All tasks.
#define queue_lock_fail_reweight_mask ((1ULL << task_type_count) - 1)

// All except drifts.
//#define queue_lock_fail_reweight_mask (((1ULL << task_type_count) - 1) & ~((1ULL << task_type_drift_part) | (1ULL << task_type_drift_gpart)))

// All except send and recv.
//#define queue_lock_fail_reweight_mask (((1ULL << task_type_count) - 1) & ~((1ULL << task_type_send) | (1ULL << task_type_recv)))

// All except send.
//#define queue_lock_fail_reweight_mask (((1ULL << task_type_count) - 1) & ~(1ULL << task_type_send))

// All except recv.
//#define queue_lock_fail_reweight_mask (((1ULL << task_type_count) - 1) & ~(1ULL << task_type_recv))


/* Counters. */
enum {
  queue_counter_swap = 0,
  queue_counter_count,
};
extern int queue_counter[queue_counter_count];

/** Struct containing a task offset and a weight, used to build the binary heap
 * of tasks in the queue. */
struct queue_entry {
  /* The offset of the task in the task list. */
  int tid;

  /* The weight of the task. This is stored in queue as well since it may well
   * be adjusted dynamically. */
  float weight;
};

/** The queue struct. */
struct queue {

  /* The lock to access this queue. */
  swift_lock_type lock;

  /* Size, count and next element. */
  int size, count;

  /* The actual tasks to which the indices refer. */
  struct task *tasks;

  /* The task indices and weights. */
  struct queue_entry *entries;

  /* DEQ for incoming tasks. */
  int *tid_incoming;
  volatile unsigned int first_incoming, last_incoming, count_incoming;

} __attribute__((aligned(queue_struct_align)));

/* Function prototypes. */
struct task *queue_gettask(struct scheduler *s, struct queue *q,
                           const struct task *prev, int blocking, int rid);
void queue_init(struct queue *q, struct task *tasks);
void queue_insert(struct queue *q, struct task *t);
void queue_clean(struct queue *q);

void queue_dump(int nodeID, int index, FILE *file, struct queue *q);

#endif /* SWIFT_QUEUE_H */
