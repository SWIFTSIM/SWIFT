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
#include "../config.h"

/* Some standard headers. */
#include <pthread.h>

/* Function type for mappings. */
typedef void (*threadpool_map_function)(void *map_data, void *extra_data);

/* Data of a threadpool. */
struct threadpool {

  /* Number of threads in this pool. */
  int num_threads;

  /* The threads themselves. */
  pthread_t *threads;

  /* This is where threads go to rest. */
  pthread_mutex_t control_mutex, thread_mutex;
  pthread_cond_t control_cond, thread_cond;

  /* Current map data and count. */
  void *map_data, *map_extra_data;
  size_t map_data_count, map_data_size, map_data_stride;
  threadpool_map_function map_function;
  
  /* Counter for the number of threads that are done. */
  int num_threads_done;
};

/* Function prototypes. */
void threadpool_init(struct threadpool *tp, int num_threads);
void threadpool_map(struct threadpool *tp, threadpool_map_function map_function,
                    void *map_data, size_t N, int stride, void *extra_data);

#endif /* SWIFT_THREADPOOL_H */
