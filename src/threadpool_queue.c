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
  state->queues = (struct threadpool_queue *)calloc(
      tp->num_threads, sizeof(struct threadpool_queue));
  if (state->queues == NULL) error("Failed to allocate per-thread queues.");

  /* Initialize each per-thread queue. */
  for (int i = 0; i < tp->num_threads; i++) {
    if (pthread_mutex_init(&state->queues[i].lock, NULL) != 0)
      error("Failed to initialize queue mutex.");
    state->queues[i].tasks = NULL;
    state->queues[i].capacity = 0;
    state->queues[i].head = 0;
    state->queues[i].tail = 0;
  }

  /* Initialize synchronization primitives. */
  if (pthread_mutex_init(&state->sleep_lock, NULL) != 0)
    error("Failed to initialize sleep lock.");
  if (pthread_cond_init(&state->sleep_cond, NULL) != 0)
    error("Failed to initialize sleep condition variable.");

  /* Initialize counters. */
  state->sleeping_count = 0;
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
static void threadpool_queue_destroy_state(struct threadpool *tp) {

  /* Get the state. */
  struct threadpool_queue_state *state = tp->queue_state;
  if (state == NULL) return;

  /* Clean up per-thread queues. */
  for (int i = 0; i < tp->num_threads; i++) {
    pthread_mutex_destroy(&state->queues[i].lock);
    free(state->queues[i].tasks);
  }
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
 * Resizes the queue buffer to a power-of-2 size that can hold at least
 * @c need tasks. The buffer size always grows by powers of 2 for efficient
 * indexing with bitwise AND.
 *
 * Note that the queue remains allocated once grown instead of being freed
 * at the end of a call to threadpool_map_with_queue to avoid repeated
 * allocations. This is fine since of the course of the run the maximum queue
 * size will be relatively stable (and likely grow rather than shrink). It will
 * be freed when the threadpool is destroyed and threadpool_clean is called.
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
 * @brief Try to pop a task from the back of the owner's queue.
 *
 * @param queue The #threadpool_queue to pop from.
 * @param data (output) Pointer to receive the task data.
 * @param count (output) Pointer to receive the task count.
 *
 * @return 1 if a task was popped, 0 if the queue was empty.
 */
static int threadpool_queue_owner_pop(struct threadpool_queue *queue,
                                      void **data, int *count) {

  int success = 0;
  pthread_mutex_lock(&queue->lock);

  if (queue->tail - queue->head > 0) {
    queue->tail--;
    const int index = queue->tail & (queue->capacity - 1);
    *data = queue->tasks[index].data;
    *count = queue->tasks[index].count;
    success = 1;
  }

  pthread_mutex_unlock(&queue->lock);
  return success;
}

/**
 * @brief Try to steal a task from the front of another thread's queue.
 *
 * @param queue The #threadpool_queue to steal from.
 * @param data (output) Pointer to receive the task data.
 * @param count (output) Pointer to receive the task count.
 *
 * @return 1 if a task was stolen, 0 if the queue was empty.
 */
static int threadpool_queue_steal(struct threadpool_queue *queue, void **data,
                                  int *count) {

  int success = 0;
  pthread_mutex_lock(&queue->lock);

  if (queue->tail - queue->head > 0) {
    const int index = queue->head & (queue->capacity - 1);
    *data = queue->tasks[index].data;
    *count = queue->tasks[index].count;
    queue->head++;
    success = 1;
  }

  pthread_mutex_unlock(&queue->lock);
  return success;
}

/**
 * @brief Try to get a task for the given thread.
 *
 * First tries to pop from the thread's own queue, then attempts to steal
 * from other threads' queues in round-robin order.
 *
 * @param state The #threadpool_queue_state.
 * @param thread_id The ID of the requesting thread.
 * @param data (output) Pointer to receive the task data.
 * @param count (output) Pointer to receive the task count.
 *
 * @return 1 if a task was obtained, 0 if all queues are empty.
 */
static int threadpool_queue_get_task(struct threadpool_queue_state *state,
                                     int thread_id, void **data, int *count) {

  /* Try own queue first. */
  if (threadpool_queue_owner_pop(&state->queues[thread_id], data, count))
    return 1;

  /* Try stealing from other threads. */
  for (int offset = 1; offset < state->tp->num_threads; offset++) {
    const int victim_id = (thread_id + offset) % state->tp->num_threads;
    if (threadpool_queue_steal(&state->queues[victim_id], data, count))
      return 1;
  }

  return 0;
}

/**
 * @brief Worker thread main loop for queue-based processing.
 *
 * Continuously tries to get tasks from the work-stealing queues, processes
 * them, and sleeps when no work is available. Terminates when all queues
 * are empty and no tasks are in flight.
 *
 * @param tp The #threadpool.
 * @param thread_id The ID of this worker thread.
 */
static void threadpool_queue_chomp(struct threadpool *tp, int thread_id) {

  /* Store thread ID as thread-specific data. */
  int local_tid = thread_id;
  pthread_setspecific(threadpool_tid, &local_tid);

  /* Get the queue state for this threadpool. */
  struct threadpool_queue_state *state = tp->queue_state;

  while (1) {

    /* Try to get work. */
    void *data;
    int count;
    if (threadpool_queue_get_task(state, thread_id, &data, &count)) {

#ifdef SWIFT_DEBUG_THREADPOOL
      ticks tic = getticks();
#endif

      /* Execute the task. */
      state->map_function(data, count, state->map_extra_data);

#ifdef SWIFT_DEBUG_THREADPOOL
      threadpool_log(tp, thread_id, count, tic, getticks());
#endif

      /* Decrement the in-flight counter. */
      atomic_dec(&state->tasks_in_flight);
      continue;
    }

    /* No immediate work available. Check if we're done.
     * We rely on the atomic tasks_in_flight counter to avoid TOCTOU races
     * when checking individual queues. If this counter is zero, all work
     * is complete. */
    if (state->tasks_in_flight == 0) {

      /* Wake all sleeping threads so they can exit and reach the barrier. */
      pthread_mutex_lock(&state->sleep_lock);
      if (state->sleeping_count > 0) pthread_cond_broadcast(&state->sleep_cond);
      pthread_mutex_unlock(&state->sleep_lock);
      break;
    }

    /* Sleep until new work arrives. */
    pthread_mutex_lock(&state->sleep_lock);
    state->sleeping_count++;

    /* Spurious wake-ups are harmless - we'll just loop again. */
    pthread_cond_wait(&state->sleep_cond, &state->sleep_lock);

    state->sleeping_count--;
    pthread_mutex_unlock(&state->sleep_lock);
  }
}

/**
 * @brief Worker loop for queue mode (called from threadpool_runner).
 *
 * This is the public entry point called from threadpool_runner() when
 * queue mode is active.
 *
 * @param tp The #threadpool.
 * @param thread_id The ID of the calling thread.
 */
void threadpool_queue_run(struct threadpool *tp, int thread_id) {

  /* Run the internal queue worker loop. */
  threadpool_queue_chomp(tp, thread_id);
}

/**
 * @brief Clean up queue support for a threadpool (called from
 * threadpool_clean).
 *
 * @param tp The #threadpool.
 */
void threadpool_queue_clean(struct threadpool *tp) {
  threadpool_queue_destroy_state(tp);
}

/**
 * @brief Add work to the current thread's queue.
 *
 * This function allows threads to dynamically add new work items during
 * execution. Work is added to the calling thread's queue and can be stolen
 * by other idle threads.
 *
 * @param tp The #threadpool.
 * @param data_ptrs Array of pointers to data items to process.
 * @param count Number of pointers in the array.
 */
void threadpool_queue_add(struct threadpool *tp, void **data_ptrs, int count) {

  if (count <= 0) return;

  /* Get the queue state. */
  struct threadpool_queue_state *state = tp->queue_state;
  if (state == NULL)
    error("Cannot add tasks to queue: queue mode not initialized!");

  /* Get the calling thread's ID. */
  const int thread_id = threadpool_gettid();

  /* Add each item as a single-element task. */
  pthread_mutex_lock(&state->queues[thread_id].lock);
  for (int i = 0; i < count; i++)
    threadpool_queue_push_locked(&state->queues[thread_id], data_ptrs[i], 1);
  pthread_mutex_unlock(&state->queues[thread_id].lock);

  /* Update the in-flight counter. */
  atomic_add(&state->tasks_in_flight, count);

  /* Wake a sleeping thread if any are waiting. */
  pthread_mutex_lock(&state->sleep_lock);
  if (state->sleeping_count > 0) pthread_cond_signal(&state->sleep_cond);
  pthread_mutex_unlock(&state->sleep_lock);
}

/**
 * @brief Map a function over data using work-stealing queues.
 *
 * Similar to threadpool_map(), but uses work-stealing queues to enable
 * dynamic load balancing. Threads can steal work from each other, and
 * new work can be added during execution via threadpool_queue_add().
 *
 * @param tp The #threadpool.
 * @param map_function The function to apply to each data chunk.
 * @param map_data Array of data to process.
 * @param count Total number of elements in map_data.
 * @param stride Size in bytes of each element.
 * @param chunk_size Chunk size hint (see threadpool_map() for details).
 * @param extra_data Additional data passed to the map function.
 */
void threadpool_map_with_queue(struct threadpool *tp,
                               threadpool_map_function map_function,
                               void *map_data, size_t count, int stride,
                               int chunk_size, void *extra_data) {

#ifdef SWIFT_DEBUG_THREADPOOL
  ticks tic_total = getticks();
#endif

  /* Single-threaded path - just call the function directly. */
  if (tp->num_threads == 1) {
    if (count <= INT_MAX) {
      map_function(map_data, count, extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
      tp->map_function = map_function;
      threadpool_log(tp, 0, count, tic_total, getticks());
#endif
    } else {
      /* Handle counts larger than INT_MAX. */
      size_t chunk = INT_MAX;
      size_t offset = 0;
      while (offset < count) {
#ifdef SWIFT_DEBUG_THREADPOOL
        ticks tic = getticks();
#endif
        const int current_chunk = (int)min(chunk, count - offset);
        map_function((char *)map_data + (stride * offset), current_chunk,
                     extra_data);
#ifdef SWIFT_DEBUG_THREADPOOL
        threadpool_log(tp, 0, current_chunk, tic, getticks());
#endif
        offset += (size_t)current_chunk;
      }
    }
    return;
  }

  /* Get the queue state. */
  struct threadpool_queue_state *state = tp->queue_state;
  if (state == NULL)
    error("Queue state not initialized! This should not happen.");
  state->map_function = map_function;
  state->map_extra_data = extra_data;

  /* Determine the chunk size for initial task distribution. */
  int initial_chunk_size;
  if (chunk_size == threadpool_auto_chunk_size) {
    initial_chunk_size =
        max((count / (tp->num_threads * threadpool_default_chunk_ratio)), 1U);
  } else if (chunk_size == threadpool_uniform_chunk_size) {
    initial_chunk_size = (count + tp->num_threads - 1) / tp->num_threads;
  } else {
    initial_chunk_size = (chunk_size > 0) ? chunk_size : 1;
  }

  /* Distribute initial tasks round-robin across threads. */
  size_t offset = 0;
  int target_thread = 0;
  int task_count = 0;

  while (offset < count) {
    const int current_chunk =
        (int)min((size_t)initial_chunk_size, count - offset);
    void *chunk_data = (char *)map_data + (stride * offset);

    pthread_mutex_lock(&state->queues[target_thread].lock);
    threadpool_queue_push_locked(&state->queues[target_thread], chunk_data,
                                 current_chunk);
    pthread_mutex_unlock(&state->queues[target_thread].lock);

    task_count++;
    offset += (size_t)current_chunk;
    target_thread = (target_thread + 1) % tp->num_threads;
  }

  /* Update the in-flight counter. */
  atomic_add(&state->tasks_in_flight, task_count);

  /* Activate queue mode and set up the threadpool for execution. */
  state->active = 1;
  tp->use_queue = 1;
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;

  /* Start the worker threads and participate in the work. */
  swift_barrier_wait(&tp->run_barrier);
  threadpool_queue_chomp(tp, tp->num_threads - 1);
  swift_barrier_wait(&tp->wait_barrier);

  /* Deactivate queue mode. */
  state->active = 0;
  tp->use_queue = 0;

#ifdef SWIFT_DEBUG_THREADPOOL
  threadpool_log(tp, -1, count, tic_total, getticks());
#endif
}
