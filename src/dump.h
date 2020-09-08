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
#ifndef SWIFT_DUMP_H
#define SWIFT_DUMP_H

/* Config parameters. */
#include "../config.h"

#ifdef HAVE_POSIX_FALLOCATE /* Are we on a sensible platform? */

/* Standard headers */
#include <stdlib.h>

/* Local headers. */
#include "atomic.h"

/** The dump struct. */
struct dump {

  /* The memory-mapped data of this dump. */
  void *data;

  /* The size of the memory-mapped data, in bytes. */
  size_t size;

  /* The number of bytes that have been dumped. */
  atomic_size_t count;

  /* The offset of the data within the current file. */
  size_t file_offset;

  /* The file with which this memory is associated. */
  int fd;

  /* Mask containing the significant bits for page addresses. */
  size_t page_mask;
};

/* Function prototypes. */
void dump_init(struct dump *d, const char *filename, size_t size);
void dump_restart(struct dump *d, const char *filename);
void dump_ensure(struct dump *d, size_t required_size, size_t increase_size);
void dump_sync(struct dump *d);
void dump_close(struct dump *d);
void *dump_get(struct dump *d, size_t count, size_t *offset);

#endif /* HAVE_POSIX_FALLOCATE */

#endif /* SWIFT_DUMP_H */
