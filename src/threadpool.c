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
#include <string.h>
#include <stdlib.h>

/* This object's header. */
#include "threadpool.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"

void *threadpool_runner(void *data) {

  /* Our threadpool. */
  struct threadpool *tp = (struct threadpool *)data;

  /* Signal to the threadpool that this thread is up and running. */
  pthread_mutex_lock(&tp->control_mutex);
  tp->num_threads += 1;
  pthread_cond_signal(&tp->control_cond);
  pthread_mutex_unlock(&tp->control_mutex);

  /* Main loop. */
  while (1) {

    /* Wait for a signal. */
    pthread_mutex_lock(&tp->thread_mutex);
    do {
      pthread_cond_wait(&tp->thread_cond, &tp->thread_mutex);
    } while (tp->map_function == NULL);
    pthread_mutex_unlock(&tp->thread_mutex);

    /* The index of the mapping task we will work on next. */
    size_t task_ind;
    while ((task_ind = atomic_inc(&tp->map_data_count)) < tp->map_data_size) {
      tp->map_function(tp->map_data + tp->map_data_stride * task_ind,
                       tp->map_extra_data);
    }

    /* Signal to the threadpool that we are done for now. */
    pthread_mutex_lock(&tp->control_mutex);
    tp->num_threads_done += 1;
    pthread_cond_signal(&tp->control_cond);
    pthread_mutex_unlock(&tp->control_mutex);
  }
}

void threadpool_init(struct threadpool *tp, int num_threads) {

  /* We will use this to count how many threads are up and running. */
  tp->num_threads = 0;

  /* Init the threadpool mutexes. */
  if (pthread_mutex_init(&tp->control_mutex, NULL) != 0 ||
      pthread_mutex_init(&tp->thread_mutex, NULL) != 0)
    error("Failed to initialize mutexex.");
  if (pthread_cond_init(&tp->control_cond, NULL) != 0 ||
      pthread_cond_init(&tp->thread_cond, NULL) != 0)
    error("Failed to initialize condition variables.");

  /* Set the task counter to zero. */
  tp->map_data_size = 0;
  tp->map_data_count = 0;
  tp->map_data_stride = 0;
  tp->map_function = NULL;

  /* Allocate the threads. */
  if ((tp->threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads)) ==
      NULL) {
    error("Failed to allocate thread array.");
  }

  /* Create and start the threads. */
  for (int k = 0; k < num_threads; k++) {
    if (pthread_create(&tp->threads[k], NULL, &threadpool_runner, tp) != 0)
      error("Failed to create threadpool runner thread.");
  }

  /* Wait for all the threads to be up and running. */
  pthread_mutex_lock(&tp->control_mutex);
  while (tp->num_threads < num_threads) {
    pthread_cond_wait(&tp->control_cond, &tp->control_mutex);
  }
  pthread_mutex_unlock(&tp->control_mutex);
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
 * @param extra_data Addtitional pointer that will be passed to the mapping
 *        function, may contain additional data.
 */

void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, void *extra_data) {

  /* Set the map data and signal the threads. */
  pthread_mutex_lock(&tp->thread_mutex);
  tp->map_data_stride = stride;
  tp->map_data_size = N;
  tp->map_data_count = 0;
  tp->map_function = map_function;
  tp->map_data = map_data;
  tp->map_extra_data = extra_data;
  tp->num_threads_done = 0;
  pthread_cond_broadcast(&tp->thread_cond);
  pthread_mutex_unlock(&tp->thread_mutex);

  /* Wait for the threads to come home. */
  pthread_mutex_lock(&tp->control_mutex);
  while (tp->num_threads_done < tp->num_threads) {
    pthread_cond_wait(&tp->control_cond, &tp->control_mutex);
  }
  pthread_mutex_unlock(&tp->control_mutex);
}
