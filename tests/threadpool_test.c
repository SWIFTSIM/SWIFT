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

// Standard includes.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// Local includes.
#include "../src/threadpool.h"
#include "../src/atomic.h"

void map_function_first(void *map_data, void *extra_data) {
  const int input = *(int *)map_data;
  usleep(rand() % 1000000);
  printf("map_function_first: got input %i.\n", input);
  fflush(stdout);
}

void map_function_second(void *map_data, void *extra_data) {
  const int input = *(int *)map_data;
  usleep(rand() % 1000000);
  printf("map_function_second: got input %i.\n", input);
  fflush(stdout);
}

int main(int argc, char *argv[]) {

  // Some constants for this test.
  const int num_threads = 16;
  const int N = 20;
  const int num_runs = 2;

  // Create a threadpool with 8 threads.
  struct threadpool tp;
  threadpool_init(&tp, num_threads);

  // Main loop.
  for (int run = 0; run < num_runs; run++) {

    // Run over a set of integers and print them.
    int data[N];
    for (int k = 0; k < N; k++) data[k] = k;
    printf("processing integers from 0..%i.\n", N);
    fflush(stdout);
    threadpool_map(&tp, map_function_first, data, N, sizeof(int), NULL);

    // Do the same thing again, with less jobs than threads.
    printf("processing integers from 0..%i.\n", num_threads / 2);
    fflush(stdout);
    threadpool_map(&tp, map_function_second, data, num_threads / 2, sizeof(int),
                   NULL);
  }
}
