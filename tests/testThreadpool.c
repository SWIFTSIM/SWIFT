/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

#include "../config.h"

// Standard includes.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Local includes.
#include "../src/atomic.h"
#include "../src/threadpool.h"

void map_function_first(void *map_data, int num_elements, void *extra_data) {
  const int *inputs = (int *)map_data;
  for (int ind = 0; ind < num_elements; ind++) {
    int input = inputs[ind];
    usleep(rand() % 1000000);
    printf("   map_function_first: got input %i.\n", input);
    fflush(stdout);
  }
}

void map_function_second(void *map_data, int num_elements, void *extra_data) {
  const int *inputs = (int *)map_data;
  for (int ind = 0; ind < num_elements; ind++) {
    int input = inputs[ind];
    usleep(rand() % 1000000);
    printf("   map_function_second: got input %i.\n", input);
    fflush(stdout);
  }
}

int main(int argc, char *argv[]) {

  // Some constants for this test.
  const int N = 20;
  const int num_runs = 2;

  // Create threadpools with different numbers of threads.
  for (int num_thread = 1; num_thread <= 16; num_thread *= 4) {
    printf("# Creating threadpool with %d threads\n", num_thread);
    struct threadpool tp;
    threadpool_init(&tp, num_thread);

    // Main loop.
    for (int run = 0; run < num_runs; run++) {

      // Run over a set of integers and print them.
      int data[N];
      for (int k = 0; k < N; k++) data[k] = k;
      printf("1..processing integers from 0..%i.\n", N);
      fflush(stdout);
      threadpool_map(&tp, map_function_first, data, N, sizeof(int), 1, NULL);

      // Do the same thing again, with less jobs than threads.
      printf("2..processing integers from 0..%i.\n", N / 2);
      fflush(stdout);
      threadpool_map(&tp, map_function_second, data, N / 2, sizeof(int), 1,
                     NULL);

      // Do the same thing again, with a chunk size of two.
      printf("3..processing integers from 0..%i.\n", N);
      fflush(stdout);
      threadpool_map(&tp, map_function_first, data, N, sizeof(int), 2, NULL);
    }

/* If logging was enabled, dump the log. */
#ifdef SWIFT_DEBUG_THREADPOOL
    char filename[80];
    sprintf(filename, "threadpool_log-%d.txt", num_thread);
    printf("# Dumping log\n");
    threadpool_dump_log(&tp, filename, 1);
#endif

    /* Be clean */
    threadpool_clean(&tp);
    printf("\n");
  }

  return 0;
}
