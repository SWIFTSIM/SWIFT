/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_THREADPOOL_H
#define SWIFT_THREADPOOL_H

/* Config parameters. */
#include <config.h>

/* Standard headers */
#include <pthread.h>
#include <stddef.h>

/* Local includes. */
#include "barrier.h"
#include "cycle.h"

/* Local defines. */
#define threadpool_log_initial_size 1000
#define threadpool_default_chunk_ratio 7
#define threadpool_auto_chunk_size 0
#define threadpool_uniform_chunk_size -1

/* Key for thread-specific data (thread ID). */
extern pthread_key_t threadpool_tid;

/* Function type for mappings. */
typedef void (*threadpool_map_function)(void *map_data, int num_elements,
                                        void *extra_data);

/* Data for threadpool logging. */
struct mapper_log_entry {

  /* ID of the thread executing the chunk. */
  int tid;

  /* Size of the chunk processed. */
  int chunk_size;

  /* Pointer to the mapper function. */
  threadpool_map_function map_function;

  /*! Start and end time of this task */
  ticks tic, toc;
};

struct mapper_log {
  /* Log of threadpool mapper calls. */
  struct mapper_log_entry *log;

  /* Size of the allocated log. */
  int size;

  /* Number of entries in the log. */
  int count;
};

/* Task description for work-stealing queue. */
struct threadpool_queue_task {

  /* Pointer to the data to be processed. */
  void *data;

  /* Number of elements in this task. */
  int count;
};

/* Per-thread work-stealing queue. */
struct threadpool_queue {

  /* Mutex protecting this queue. */
  pthread_mutex_t lock;

  /* Buffer of tasks (power-of-2 sized for fast indexing). */
  struct threadpool_queue_task *tasks;

  /* Capacity of the task buffer (always power of 2). */
  int capacity;

  /* Head index (for stealing from the front). */
  int head;

  /* Tail index (for owner pop from the back). */
  int tail;
};

/* Per-threadpool queue state. */
struct threadpool_queue_state {

  /* The threadpool this state belongs to. */
  struct threadpool *tp;

  /* Array of per-thread queues. */
  struct threadpool_queue *queues;

  /* Mutex protecting the sleep condition variable. */
  pthread_mutex_t sleep_lock;

  /* Condition variable for waking sleeping threads. */
  pthread_cond_t sleep_cond;

  /* Number of threads currently sleeping. */
  volatile int sleeping_count;

  /* Is the queue system currently active for this threadpool? */
  int active;

  /* Number of tasks currently in flight (being processed or queued). */
  volatile int tasks_in_flight;

  /* Cached map function for queue mode. */
  threadpool_map_function map_function;

  /* Cached extra data for queue mode. */
  void *map_extra_data;
};

/* Data of a threadpool. */
struct threadpool {

  /* The threads themselves. */
  pthread_t *threads;

  /* This is where threads go to rest. */
  swift_barrier_t wait_barrier;
  swift_barrier_t run_barrier;

  /* Current map data and count. */
  void *map_data, *map_extra_data;
  volatile size_t map_data_count, map_data_size, map_data_stride;
  volatile ptrdiff_t map_data_chunk;
  volatile threadpool_map_function map_function;

  /* Number of threads in this pool. */
  int num_threads;

  /* Counter for the number of threads that are done. */
  volatile int num_threads_running;

  /* Flag indicating if queue mode is active. */
  volatile int use_queue;

  /* Queue state for work-stealing mode. */
  struct threadpool_queue_state *queue_state;

#ifdef SWIFT_DEBUG_THREADPOOL
  struct mapper_log *logs;
#endif
};

/* Function prototypes. */
void threadpool_init(struct threadpool *tp, int num_threads);
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, int chunk,
                    void *extra_data);
int threadpool_gettid(void);
void threadpool_clean(struct threadpool *tp);

/* Queue-based threadpool prototypes. */
void threadpool_queue_init(struct threadpool *tp);
void threadpool_queue_clean(struct threadpool *tp);
void threadpool_queue_chomp(struct threadpool *tp, int thread_id);
void threadpool_map_with_queue(struct threadpool *tp,
                               threadpool_map_function map_function,
                               void *map_data, size_t N, int stride, int chunk,
                               void *extra_data);
void threadpool_queue_add(struct threadpool *tp, void **ptrs, int n_ptrs);
int threadpool_queue_get_waiting(struct threadpool *tp);

#ifdef HAVE_SETAFFINITY
void threadpool_set_affinity_mask(cpu_set_t *entry_affinity);
#endif
#ifdef SWIFT_DEBUG_THREADPOOL
void threadpool_reset_log(struct threadpool *tp);
void threadpool_dump_log(struct threadpool *tp, const char *filename,
                         int reset);
void threadpool_log(struct threadpool *tp, int tid, size_t chunk_size,
                    ticks tic, ticks toc);
#endif

#endif /* SWIFT_THREADPOOL_H */
