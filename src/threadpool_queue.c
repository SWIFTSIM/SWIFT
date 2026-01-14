/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 William Roper (w.roper@sussex.ac.uk)
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
#include <limits.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "threadpool.h"

/* Local headers. */
#include "atomic.h"
#include "clocks.h"
#include "error.h"
#include "minmax.h"

/**
 * @brief Initialize queue support for a threadpool (called from
 * threadpool_init).
 *
 * @param tp The #threadpool.
 */
void threadpool_queue_init(struct threadpool *tp) {

  /* Allocate the state structure. */
  struct threadpool_queue_state *state =
      (struct threadpool_queue_state *)malloc(sizeof(*state));
  if (state == NULL) error("Failed to allocate threadpool queue state.");

  /* Initialize the state. */
  state->tp = tp;

  /* Initialize the shared queue. */
  /* We use the first queue slot in the structure as the shared queue. */
  state->queues = (struct threadpool_queue *)calloc(
      1, sizeof(struct threadpool_queue));
  if (state->queues == NULL) error("Failed to allocate shared queue.");

  if (pthread_mutex_init(&state->queues[0].lock, NULL) != 0)
    error("Failed to initialize queue mutex.");
  state->queues[0].tasks = NULL;
  state->queues[0].capacity = 0;
  state->queues[0].head = 0;
  state->queues[0].tail = 0;

  /* Initialize synchronization primitives. */
  if (pthread_mutex_init(&state->sleep_lock, NULL) != 0)
    error("Failed to initialize sleep lock.");
  if (pthread_cond_init(&state->sleep_cond, NULL) != 0)
    error("Failed to initialize sleep condition variable.");

  /* Initialize counters. */
  state->sleeping_count = 0;
  state->is_sleeping = NULL; /* Not used in single queue mode */
  state->active = 0;
  state->tasks_in_flight = 0;
  state->map_function = NULL;
  state->map_extra_data = NULL;

  /* Store in threadpool. */
  tp->queue_state = state;
}

/**
 * @brief Destroy the queue state for a threadpool.
 *
 * @param tp The #threadpool.
 */
void threadpool_queue_clean(struct threadpool *tp) {

  /* Get the state. */
  struct threadpool_queue_state *state = tp->queue_state;
  if (state == NULL) return;

  /* Clean up the shared queue. */
  pthread_mutex_destroy(&state->queues[0].lock);
  free(state->queues[0].tasks);
  free(state->queues);

  /* Clean up synchronization primitives. */
  pthread_mutex_destroy(&state->sleep_lock);
  pthread_cond_destroy(&state->sleep_cond);

  /* Free the state. */
  free(state);
  tp->queue_state = NULL;
}

/**
 * @brief Ensure a queue has capacity for at least the requested number of
 * tasks.
 *
 * @param queue The #threadpool_queue to resize.
 * @param need The minimum number of tasks the queue must hold.
 */
static void threadpool_queue_reserve(struct threadpool_queue *queue, int need) {

  /* Already have enough capacity. */
  if (queue->capacity >= need) return;

  /* Find the next power-of-2 size. */
  int new_capacity = queue->capacity ? queue->capacity : 64;
  while (new_capacity < need) new_capacity <<= 1;

  /* Allocate new buffer. */
  struct threadpool_queue_task *new_tasks =
      (struct threadpool_queue_task *)malloc(
          sizeof(struct threadpool_queue_task) * new_capacity);
  if (new_tasks == NULL) error("Failed to allocate queue task buffer.");

  /* Copy existing tasks to new buffer. */
  const int count = queue->tail - queue->head;
  for (int i = 0; i < count; i++) {
    const int old_index = (queue->head + i) & (queue->capacity - 1);
    new_tasks[i] = queue->tasks[old_index];
  }

  /* Replace the old buffer. */
  free(queue->tasks);
  queue->tasks = new_tasks;
  queue->capacity = new_capacity;
  queue->head = 0;
  queue->tail = count;
}

/**
 * @brief Push a task onto a queue (must hold the queue lock).
 *
 * @param queue The #threadpool_queue.
 * @param data Pointer to the task data.
 * @param count Number of elements in this task.
 */
static void threadpool_queue_push_locked(struct threadpool_queue *queue,
                                         void *data, int count) {

  /* Ensure we have space. */
  const int current_count = queue->tail - queue->head;
  if (queue->capacity - current_count < 1)
    threadpool_queue_reserve(queue, queue->capacity ? 2 * queue->capacity : 64);

  /* Add the task. */
  const int index = queue->tail & (queue->capacity - 1);
  queue->tasks[index].data = data;
  queue->tasks[index].count = count;
  queue->tail++;
}

/**
 * @brief Get tasks from the shared queue.
 *
 * @param state The #threadpool_queue_state.
 * @param tasks Array to receive the tasks.
 * @param max_tasks Maximum number of tasks to get.
 *
 * @return Number of tasks obtained.
 */
static int threadpool_queue_get_task(struct threadpool_queue_state *state,
                                     struct threadpool_queue_task *tasks,
                                     int max_tasks) {

  struct threadpool_queue *queue = &state->queues[0];
  int popped = 0;

  pthread_mutex_lock(&queue->lock);

  /* Check if queue has tasks. */
  while (popped < max_tasks && queue->tail - queue->head > 0) {
    const int index = queue->head & (queue->capacity - 1);
    tasks[popped] = queue->tasks[index];
    queue->head++;
    popped++;
  }

  pthread_mutex_unlock(&queue->lock);
  return popped;
}

/**
 * @brief Worker thread main loop for queue-based processing.
 *
 * @param tp The #threadpool.
 * @param thread_id The ID of this worker thread.
 */
void threadpool_queue_chomp(struct threadpool *tp, int thread_id) {

  /* Store thread ID as thread-specific data. */
  int local_tid = thread_id;
  pthread_setspecific(threadpool_tid, &local_tid);

  /* Get the queue state for this threadpool. */
  struct threadpool_queue_state *state = tp->queue_state;

  /* Buffers for batch processing. */
  struct threadpool_queue_task tasks[128];
  void *batch_buffer[128];

  while (1) {

    /* Determine the desired chunk size. */
    ptrdiff_t chunk_size;
    if (tp->map_data_chunk == threadpool_uniform_chunk_size) {
      chunk_size = ((thread_id + 1) * tp->map_data_size / tp->num_threads) -
                   (thread_id * tp->map_data_size / tp->num_threads);
    } else {
      chunk_size =
          (tp->map_data_size - tp->map_data_count) / (2 * tp->num_threads);
      if (chunk_size > tp->map_data_chunk) chunk_size = tp->map_data_chunk;
    }
    if (chunk_size < 1) chunk_size = 1;
    if (chunk_size > INT_MAX) chunk_size = INT_MAX;

    /* Clamp the number of tasks we fetch from the queue. */
    int desired_queue_tasks = (int)chunk_size;
    if (desired_queue_tasks > 128) desired_queue_tasks = 128;

    /* Try to get tasks from the shared queue first (donated work). */
    int n_tasks =
        threadpool_queue_get_task(state, tasks, desired_queue_tasks);

    /* If we got queue tasks, process them. */
    if (n_tasks > 0) {
      int batch_count = 0;

      /* Execute the tasks. */
      for (int i = 0; i < n_tasks; i++) {

        if (tasks[i].count == 1) {
          batch_buffer[batch_count++] = tasks[i].data;
        } else {
          /* Flush batch if necessary. */
          if (batch_count > 0) {
#ifdef SWIFT_DEBUG_THREADPOOL
            ticks tic = getticks();
#endif
            state->map_function(batch_buffer, batch_count,
                                state->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
            threadpool_log(tp, thread_id, batch_count, tic, getticks());
#endif
            batch_count = 0;
          }

          /* Process large chunk. */
#ifdef SWIFT_DEBUG_THREADPOOL
          ticks tic = getticks();
#endif
          state->map_function(tasks[i].data, tasks[i].count,
                              state->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
          threadpool_log(tp, thread_id, tasks[i].count, tic, getticks());
#endif
        }
      }

      /* Flush remaining batched items. */
      if (batch_count > 0) {
#ifdef SWIFT_DEBUG_THREADPOOL
        ticks tic = getticks();
#endif
        state->map_function(batch_buffer, batch_count, state->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
        threadpool_log(tp, thread_id, batch_count, tic, getticks());
#endif
      }

      /* Decrement in-flight counter. */
      atomic_sub(&state->tasks_in_flight, n_tasks);
      continue;
    }

    /* Queue is empty, try to get a chunk from the main range (threadpool_map style). */
    size_t task_ind = atomic_add(&tp->map_data_count, chunk_size);
    if (task_ind < tp->map_data_size) {
      /* We got a valid range! */
      if (task_ind + chunk_size > tp->map_data_size)
        chunk_size = tp->map_data_size - task_ind;

#ifdef SWIFT_DEBUG_THREADPOOL
      ticks tic = getticks();
#endif
      state->map_function((char *)tp->map_data + (tp->map_data_stride * task_ind),
                          (int)chunk_size, state->map_extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
      threadpool_log(tp, thread_id, chunk_size, tic, getticks());
#endif

      /* Update in-flight counter to reflect the work we just picked up. 
       * Note: In threadpool_map, we don't track in-flight tasks this way, 
       * but for the queue sleep logic to work, we need to know threads are busy. 
       * However, the loop continues and eventually we hit the sleep check. 
       * The trick is that `tasks_in_flight` tracks *queue* tasks. 
       * The `map_data_count` tracks the main array. 
       * A thread is "done" only when `map_data_count >= size` AND `tasks_in_flight == 0`.
       */
      continue;
    }

    /* No immediate work available. Check if we are done or need to sleep. */
    pthread_mutex_lock(&state->sleep_lock);

    /* Re-check both conditions under lock. */
    if (state->tasks_in_flight == 0 && tp->map_data_count >= tp->map_data_size) {
      if (state->sleeping_count > 0) pthread_cond_broadcast(&state->sleep_cond);
      pthread_mutex_unlock(&state->sleep_lock);
      break;
    }

    /* Sleep until work or termination. */
    state->sleeping_count++;
    pthread_cond_wait(&state->sleep_cond, &state->sleep_lock);
    state->sleeping_count--;
    pthread_mutex_unlock(&state->sleep_lock);
  }
}

/**
 * @brief Add work to the shared queue.
 *
 * @param tp The #threadpool.
 * @param data_ptrs Array of pointers to data items to process.
 * @param count Number of pointers in the array.
 */
void threadpool_queue_add(struct threadpool *tp, void **data_ptrs, int count) {

#ifdef SWIFT_DEBUG_CHECKS
  if (count <= 0) {
    error("Trying to add no work to threadpool queue.");
  }
#endif

  struct threadpool_queue_state *state = tp->queue_state;

  /* Add items to the shared queue. */
  pthread_mutex_lock(&state->queues[0].lock);
  for (int i = 0; i < count; i++)
    threadpool_queue_push_locked(&state->queues[0], data_ptrs[i], 1);
  pthread_mutex_unlock(&state->queues[0].lock);

  /* Update the in-flight counter. */
  atomic_add(&state->tasks_in_flight, count);

  /* Wake any sleeping threads. */
  pthread_mutex_lock(&state->sleep_lock);
  if (state->sleeping_count > 0) pthread_cond_broadcast(&state->sleep_cond);
  pthread_mutex_unlock(&state->sleep_lock);
}

/**
 * @brief Get the number of currently sleeping threads.
 *
 * @param tp The #threadpool.
 * @return Number of sleeping threads.
 */
int threadpool_queue_get_sleeping_count(struct threadpool *tp) {
  struct threadpool_queue_state *state = tp->queue_state;
  return (state == NULL) ? 0 : state->sleeping_count;
}

/**
 * @brief Map a function over data using a shared queue.
 *
 * @param tp The #threadpool.
 * @param map_function The function to apply to each data chunk.
 * @param map_data Array of data to process.
 * @param count Total number of elements in map_data.
 * @param stride Size in bytes of each element.
 * @param chunk_size Chunk size hint.
 * @param extra_data Additional data passed to the map function.
 */
void threadpool_map_with_queue(struct threadpool *tp,
                               threadpool_map_function map_function,
                               void *map_data, size_t count, int stride,
                               int chunk_size, void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic_total = getticks();
#endif

  struct threadpool_queue_state *state = tp->queue_state;

  /* Set the map data and signal the threads (copied from threadpool_map). */
  tp->map_data_stride = stride;
  tp->map_data_size = count;
  tp->map_data_count = 0;
  if (chunk_size == threadpool_auto_chunk_size) {
    tp->map_data_chunk =
        max((count / (tp->num_threads * threadpool_default_chunk_ratio)), 1U);
  } else if (chunk_size == threadpool_uniform_chunk_size) {
    tp->map_data_chunk = threadpool_uniform_chunk_size;
  } else {
    tp->map_data_chunk = chunk_size;
  }
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;

  /* Attach the map function to the queue state as well. */
  state->map_function = map_function;
  state->map_extra_data = extra_data;
  
  /* Reset in-flight tasks (queue is initially empty). */
  state->tasks_in_flight = 0;

  /* Activate and run. */
  state->active = 1;
  tp->use_queue = 1;

  swift_barrier_wait(&tp->run_barrier);
  threadpool_queue_chomp(tp, tp->num_threads - 1);
  swift_barrier_wait(&tp->wait_barrier);

  state->active = 0;
  tp->use_queue = 0;

#ifdef SWIFT_DEBUG_THREADPOOL
  threadpool_log(tp, -1, count, tic_total, getticks());
#endif
}