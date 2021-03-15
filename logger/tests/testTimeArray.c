/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

#include "logger_time.h"

#include <stdlib.h>
#include <time.h>

#define NUMBER_OF_ELEMENT 10000
#define TIME_BASE 0.04
#define OFFSET_BASE 1000

int main(int argc, char *argv[]) {

  /* Check that we are really testing the reallocation */
  if (NUMBER_OF_ELEMENT < LOGGER_TIME_INIT_SIZE) {
    error("Not testing the reallocation.");
  }

  /* Fix the random seed in order to reproduce the results */
  srand(100);

  /* Initialize the time array */
  struct time_array times;
  time_array_init(&times, /* initial_size */ 1024);

  /* Add elements */
  for (size_t i = 0; i < NUMBER_OF_ELEMENT; i++) {
    integertime_t int_time = i;
    double time = i * TIME_BASE;
    size_t offset = i * OFFSET_BASE;

    time_array_append(&times, int_time, time, offset);
  }

  /* Check the elements */
  for (size_t i = 0; i < NUMBER_OF_ELEMENT; i++) {
    integertime_t int_time = i;
    double time = i * TIME_BASE;
    size_t offset = i * OFFSET_BASE;

    /* Ensure that we can get the correct offset when looking
       in between the records. */
    int r = rand() % OFFSET_BASE;
    size_t read_offset = offset + r;

    /* The offset cannot be larger than the largest one */
    if (i == NUMBER_OF_ELEMENT - 1) {
      read_offset = offset;
    }

    /* Get the index from the offset */
    size_t ind = time_array_get_index(&times, read_offset);

    /* Check the values obtained */
    assert(i == ind);
    assert(int_time == times.records[ind].int_time);
    assert(time == times.records[ind].time);
    assert(offset == times.records[ind].offset);
  }

  time_array_free(&times);
  return 0;
}
