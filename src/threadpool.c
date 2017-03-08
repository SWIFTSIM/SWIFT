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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "threadpool.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"

void *threadpool_runner(void *data) {

  /* Our threadpool. */
  struct threadpool *tp = (struct threadpool *)data;

  /* Main loop. */
  while (1) {

    /* Let the controller know that this thread is waiting. */
    pthread_mutex_lock(&tp->thread_mutex);
    tp->num_threads_waiting += 1;
    if (tp->num_threads_waiting == tp->num_threads) {
      pthread_cond_signal(&tp->control_cond);
    }

    /* Wait for the controller. */
    pthread_cond_wait(&tp->thread_cond, &tp->thread_mutex);
    tp->num_threads_waiting -= 1;
    tp->num_threads_running += 1;
    if (tp->num_threads_running == tp->num_threads) {
      pthread_cond_signal(&tp->control_cond);
    }
    pthread_mutex_unlock(&tp->thread_mutex);

    /* The index of the mapping task we will work on next. */
    while (1) {
      /* Desired chunk size. */
      size_t chunk_size =
          (tp->map_data_size - tp->map_data_count) / (2 * tp->num_threads);
      if (chunk_size > tp->map_data_chunk) chunk_size = tp->map_data_chunk;
      if (chunk_size < 1) chunk_size = 1;

      /* Get a chunk and check its size. */
      size_t task_ind = atomic_add(&tp->map_data_count, chunk_size);
      if (task_ind >= tp->map_data_size) break;
      if (task_ind + chunk_size > tp->map_data_size)
        chunk_size = tp->map_data_size - task_ind;

      /* Call the mapper function. */
      tp->map_function((char *)tp->map_data + (tp->map_data_stride * task_ind),
                       chunk_size, tp->map_extra_data);
    }
  }
}

/**
 * @brief Initialises the #threadpool with a given number of threads.
 *
 * @param tp The #threadpool.
 * @param num_threads The number of threads.
 */
void threadpool_init(struct threadpool *tp, int num_threads) {

  /* Initialize the thread counters. */
  tp->num_threads = num_threads;
  tp->num_threads_waiting = 0;

  /* If there is only a single thread, do nothing more as of here as
     we will just do work in the (blocked) calling thread. */
  if (num_threads == 1) return;

  /* Init the threadpool mutexes. */
  if (pthread_mutex_init(&tp->thread_mutex, NULL) != 0)
    error("Failed to initialize mutexex.");
  if (pthread_cond_init(&tp->control_cond, NULL) != 0 ||
      pthread_cond_init(&tp->thread_cond, NULL) != 0)
    error("Failed to initialize condition variables.");

  /* Set the task counter to zero. */
  tp->map_data_size = 0;
  tp->map_data_count = 0;
  tp->map_data_stride = 0;
  tp->map_data_chunk = 0;
  tp->map_function = NULL;

  /* Allocate the threads. */
  if ((tp->threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads)) ==
      NULL) {
    error("Failed to allocate thread array.");
  }

  /* Create and start the threads. */
  pthread_mutex_lock(&tp->thread_mutex);
  for (int k = 0; k < num_threads; k++) {
    if (pthread_create(&tp->threads[k], NULL, &threadpool_runner, tp) != 0)
      error("Failed to create threadpool runner thread.");
  }

  /* Wait for all the threads to be up and running. */
  while (tp->num_threads_waiting < tp->num_threads) {
    pthread_cond_wait(&tp->control_cond, &tp->thread_mutex);
  }
  pthread_mutex_unlock(&tp->thread_mutex);
}

/**
 * @brief Map a function to an array of data in parallel using a #threadpool.
 *
 * The function @c map_function is called on each element of @c map_data
 * in parallel.
 *
 * @param tp The #threadpool on which to run.
 * @param map_function The function that will be applied to the map data.
 * @param map_data The data on which the mapping function will be called.
 * @param N Number of elements in @c map_data.
 * @param stride Size, in bytes, of each element of @c map_data.
 * @param chunk Number of map data elements to pass to the function at a time.
 * @param extra_data Addtitional pointer that will be passed to the mapping
 *        function, may contain additional data.
 */
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, int chunk,
                    void *extra_data) {

  /* If we just have a single thread, call the map function directly. */
  if (tp->num_threads == 1) {
    map_function(map_data, N, extra_data);
    return;
  }

  /* Set the map data and signal the threads. */
  pthread_mutex_lock(&tp->thread_mutex);
  tp->map_data_stride = stride;
  tp->map_data_size = N;
  tp->map_data_count = 0;
  tp->map_data_chunk = chunk;
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_running = 0;
  pthread_cond_broadcast(&tp->thread_cond);

  /* Wait for all the threads to be up and running. */
  while (tp->num_threads_running < tp->num_threads) {
    pthread_cond_wait(&tp->control_cond, &tp->thread_mutex);
  }

  /* Wait for all threads to be done. */
  while (tp->num_threads_waiting < tp->num_threads) {
    pthread_cond_wait(&tp->control_cond, &tp->thread_mutex);
  }
  pthread_mutex_unlock(&tp->thread_mutex);
}

/**
 * @brief Frees up the memory allocated for this #threadpool.
 */
void threadpool_clean(struct threadpool *tp) { free(tp->threads); }
