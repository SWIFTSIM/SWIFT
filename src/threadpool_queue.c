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
  state->is_sleeping = (int *)calloc(tp->num_threads, sizeof(int));
  if (state->is_sleeping == NULL) error("Failed to allocate is_sleeping array.");
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
  free((void *)state->is_sleeping);
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
 * @brief Try to pop tasks from the back of the owner's queue (LIFO).
 *
 * The owner thread pops from the tail (back) of its own queue in LIFO (Last In
 * First Out) order. This provides good cache locality by working on
 * recently-added tasks first. The queue is logically a circular buffer indexed
 * by head and tail, where head <= tail. When a task is popped, tail is
 * decremented.
 *
 * @param queue The #threadpool_queue to pop from.
 * @param tasks Array to receive the tasks.
 * @param max_tasks Maximum number of tasks to pop.
 *
 * @return Number of tasks popped.
 */
static int threadpool_queue_owner_pop(struct threadpool_queue *queue,
                                      struct threadpool_queue_task *tasks,
                                      int max_tasks) {

  /* Lock before we try to access the queue. */
  int popped = 0;
  pthread_mutex_lock(&queue->lock);

  /* Check if queue has tasks (tail > head means non-empty). */
  while (popped < max_tasks && queue->tail - queue->head > 0) {
    queue->tail--;
    const int index = queue->tail & (queue->capacity - 1);
    tasks[popped] = queue->tasks[index];
    popped++;
  }

  /* Release the lock, we're open for theft. */
  pthread_mutex_unlock(&queue->lock);
  return popped;
}

/**
 * @brief Try to steal tasks from the front of another thread's queue (FIFO).
 *
 * Stealing threads take from the head (front) of the victim's queue in FIFO
 * (First In First Out) order, opposite from the owner which pops from the tail.
 * This reduces contention - the owner and thieves work from opposite ends of
 * the queue. Stealing the oldest tasks also promotes better load balancing, as
 * these are often larger chunks that were added first. When a task is stolen,
 * head is incremented.
 *
 * @param queue The #threadpool_queue to steal from.
 * @param tasks Array to receive the tasks.
 * @param max_tasks Maximum number of tasks to steal.
 *
 * @return Number of tasks stolen.
 */
static int threadpool_queue_steal(struct threadpool_queue *queue,
                                  struct threadpool_queue_task *tasks,
                                  int max_tasks) {

  /* Lock before we try to access the queue. */
  int stolen = 0;
  pthread_mutex_lock(&queue->lock);

  /* Check if queue has tasks (tail > head means non-empty). */
  while (stolen < max_tasks && queue->tail - queue->head > 0) {
    const int index = queue->head & (queue->capacity - 1);
    tasks[stolen] = queue->tasks[index];
    queue->head++;
    stolen++;
  }

  /* Mischief managed, release the lock. */
  pthread_mutex_unlock(&queue->lock);
  return stolen;
}

/**
 * @brief Try to get tasks for the given thread.
 *
 * First tries to pop from the thread's own queue, then attempts to steal
 * from other threads' queues in round-robin order.
 *
 * @param state The #threadpool_queue_state.
 * @param thread_id The ID of the requesting thread.
 * @param tasks Array to receive the tasks.
 * @param max_tasks Maximum number of tasks to get.
 *
 * @return Number of tasks obtained.
 */
static int threadpool_queue_get_task(struct threadpool_queue_state *state,
                                     int thread_id,
                                     struct threadpool_queue_task *tasks,
                                     int max_tasks) {

  /* Try own queue first. */
  int got = threadpool_queue_owner_pop(&state->queues[thread_id], tasks,
                                       max_tasks);
  if (got > 0) return got;

  /* Try stealing from other threads. */
  for (int offset = 1; offset < state->tp->num_threads; offset++) {
    const int victim_id = (thread_id + offset) % state->tp->num_threads;
    got = threadpool_queue_steal(&state->queues[victim_id], tasks, max_tasks);
    if (got > 0)
      /* Heist successful! */
      return got;
  }

  /* No work available. */
  return 0;
}

/**
 * @brief Worker thread main loop for queue-based processing.
 *
 * This is called from threadpool_runner if the threadpool is configured
 * to use queues.
 *
 * Continuously tries to get tasks from the work-stealing queues, processes
 * them, and sleeps when no work is available. Terminates when all queues
 * are empty and no tasks are in flight.
 *
 * @param tp The #threadpool.
 * @param thread_id The ID of this worker thread.
 */
/**
 * @brief Worker thread main loop for queue-based processing.
 *
 * This is called from threadpool_runner if the threadpool is configured
 * to use queues.
 *
 * Continuously tries to get tasks from the work-stealing queues, processes
 * them, and sleeps when no work is available. Terminates when all queues
 * are empty and no tasks are in flight.
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
    int desired;
    if (tp->map_data_chunk == threadpool_uniform_chunk_size) {
      desired = 128; /* High throughput for uniform mode */
    } else {
      /* Dynamic sizing: scale with remaining work (Guided scheduling) */
      const int remaining = state->tasks_in_flight;
      desired = remaining / (2 * tp->num_threads);

      /* Clamp to the limit set in map_with_queue */
      if (desired > tp->map_data_chunk) desired = tp->map_data_chunk;
      if (desired < 1) desired = 1;
    }

    /* Clamp to local buffer size */
    if (desired > 128) desired = 128;

    /* Try to get tasks in bulk. */
    int n_tasks = threadpool_queue_get_task(state, thread_id, tasks, desired);

    if (n_tasks > 0) {
      int batch_count = 0;

      /* Execute the tasks. */
      for (int i = 0; i < n_tasks; i++) {

        /* If this is a single item, we can batch it. */
        if (tasks[i].count == 1) {
          batch_buffer[batch_count++] = tasks[i].data;
        } else {
          /* If we have a chunk, flush the batch first. */
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

          /* Now process the large chunk directly. */
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

      /* Decrement in-flight counter by the number of tasks retrieved. */
      atomic_sub(&state->tasks_in_flight, n_tasks);
      continue;
    }

    /* No immediate work available. Check if we are done or need to sleep. */
    pthread_mutex_lock(&state->sleep_lock);

    /* Re-check under lock to avoid lost wakeup race. */
    if (state->tasks_in_flight == 0) {
      if (state->sleeping_count > 0) pthread_cond_broadcast(&state->sleep_cond);
      pthread_mutex_unlock(&state->sleep_lock);
      break;
    }

    /* Sleep until work or termination. */
    state->sleeping_count++;
    state->is_sleeping[thread_id] = 1;
    pthread_cond_wait(&state->sleep_cond, &state->sleep_lock);
    state->is_sleeping[thread_id] = 0;
    state->sleeping_count--;
    pthread_mutex_unlock(&state->sleep_lock);
  }
}

/**
 * @brief Add work to a specific thread's queue.
 *
 * This function allows threads to dynamically add new work items during
 * execution. Work is added to the target thread's queue and can be stolen
 * by other idle threads if that thread doesn't pick it up first.
 *
 * @param tp The #threadpool.
 * @param data_ptrs Array of pointers to data items to process.
 * @param count Number of pointers in the array.
 * @param tid The target thread ID to add the work to.
 */
void threadpool_queue_add(struct threadpool *tp, void **data_ptrs, int count,
                          int tid) {

#ifdef SWIFT_DEBUG_CHECKS
  /* No work to add, what are we doing? */
  if (count <= 0) {
    error("Trying to add no work to threadpool queue.");
  }
#endif

  /* Get the queue state. */
  struct threadpool_queue_state *state = tp->queue_state;

  /* Add each item to the queue. */
  pthread_mutex_lock(&state->queues[tid].lock);
  for (int i = 0; i < count; i++)
    threadpool_queue_push_locked(&state->queues[tid], data_ptrs[i], 1);
  pthread_mutex_unlock(&state->queues[tid].lock);

  /* Update the in-flight counter. */
  atomic_add(&state->tasks_in_flight, count);

  /* Wake any sleeping threads. */
  pthread_mutex_lock(&state->sleep_lock);
  if (state->sleeping_count > 0) pthread_cond_broadcast(&state->sleep_cond);
  pthread_mutex_unlock(&state->sleep_lock);
}

/**
 * @brief Find a thread that is currently waiting for work.
 *
 * @param tp The #threadpool.
 * @return The ID of a waiting thread, or -1 if none are waiting.
 */
int threadpool_queue_get_waiting_tid(struct threadpool *tp) {
  struct threadpool_queue_state *state = tp->queue_state;
  if (state == NULL || state->sleeping_count == 0) return -1;

  for (int i = 0; i < tp->num_threads; i++) {
    if (state->is_sleeping[i]) return i;
  }
  return -1;
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

  /* Get the queue state. */
  struct threadpool_queue_state *state = tp->queue_state;

  /* Attach the map function and extra data to the state. */
  state->map_function = map_function;
  state->map_extra_data = extra_data;

  /* Determine the chunk size for initial task distribution. */
  int initial_chunk_size;
  if (chunk_size == threadpool_auto_chunk_size) {
    tp->map_data_chunk =
        max((count / (tp->num_threads * threadpool_default_chunk_ratio)), 1U);
    initial_chunk_size = tp->map_data_chunk;
  } else if (chunk_size == threadpool_uniform_chunk_size) {
    tp->map_data_chunk = threadpool_uniform_chunk_size;
    initial_chunk_size = (count + tp->num_threads - 1) / tp->num_threads;
  } else {
    tp->map_data_chunk = chunk_size;
    initial_chunk_size = (chunk_size > 0) ? chunk_size : 1;
  }

  /* Distribute initial tasks using block distribution. */
  int task_count = 0;
  for (int i = 0; i < tp->num_threads; i++) {
    const size_t thread_offset = (i * count) / tp->num_threads;
    const size_t thread_count =
        ((i + 1) * count) / tp->num_threads - thread_offset;

    /* No work for this thread. */
    if (thread_count == 0) continue;

    /* Lock the queue while we add tasks (almost certainly unnecessary here as
     * no other thread should be touching this queue yet, but for safety). */
    pthread_mutex_lock(&state->queues[i].lock);

    /* Break the thread's work into chunks and add to the queue. */
    size_t local_offset = 0;
    while (local_offset < thread_count) {
      const int current_chunk =
          (int)min((size_t)initial_chunk_size, thread_count - local_offset);
      void *chunk_data =
          (char *)map_data + (stride * (thread_offset + local_offset));

      threadpool_queue_push_locked(&state->queues[i], chunk_data,
                                   current_chunk);

      task_count++;
      local_offset += (size_t)current_chunk;
    }

    /* Done with this queue. */
    pthread_mutex_unlock(&state->queues[i].lock);
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
