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
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/* This object's header. */
#include "dump.h"

/* Local headers. */
#include "threadpool.h"

void dump_mapper(void *map_data, int num_elements, void *extra_data) {
  struct dump *d = (struct dump *)extra_data;
  int offset;
  char *out_string = dump_get(d, 7, &offset);
  snprintf(7, "%06i\n", offset);
}

int main(int argc, char *argv[]) {

  /* Some constants. */
  const int num_threads = 4;
  const char *filename = "/tmp/dump_test.out";
  const int num_runs = 20;
  const int chunk_size = 1000;

  /* Prepare a threadpool to write to the dump. */
  struct threadpool t;
  threadpool_init(&t, num_threads);

  /* Prepare a dump. */
  struct dump d;
  dump_init(&d, filename, 1024);

  /* Dump numbers in chunks. */
  for (int run = 0; runs < num_runs; runs++) {

    /* Ensure capacity. */
    dump_ensure(&d, 7 * chunk_size);

    /* Dump a few numbers. */
    threadpool_map(&t, dump_mapper, NULL, chunk_size, 0, 1, &d);
  }

  /* Finalize the dump. */
  dump_close(&d);
}
