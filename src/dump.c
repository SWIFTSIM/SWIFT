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

#ifdef HAVE_POSIX_FALLOCATE

/* This object's header. */
#include "dump.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"

/* Some standard headers. */
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/**
 * @brief Obtain a chunk of memory from a dump.
 *
 * @param d The #dump.
 * @param count The number of bytes requested.
 * @param offset The offset of the returned memory address within the dump file.
 * @return A pointer to the memory-mapped chunk of data.
 */
void *dump_get(struct dump *d, size_t count, size_t *offset) {
  size_t local_offset = atomic_add(&d->count, count);
#ifdef SWIFT_DEBUG_CHECKS
  if (d->count > d->size) error("Dump file is too small.");
#endif
  *offset = local_offset + d->file_offset;
  return (char *)d->data + local_offset;
}

/**
 * @brief Ensure that at least size bytes are available in the #dump.
 *
 * @param d The #dump.
 * @param required_size The required size for the #dump
 * @param increase_size If not enough size, increase by this amount
 */
void dump_ensure(struct dump *d, size_t required_size, size_t increase_size) {

  /* If we have enough space already, just bail. */
  if (d->size - d->count > required_size) return;

  /* Unmap the current data. */
  if (munmap(d->data, d->size) != 0) {
    error("Failed to unmap %zi bytes of dump data (%s).", d->size,
          strerror(errno));
  }

  /* Update the size and count. */
  const size_t trunc_count = d->count & d->page_mask;
  d->file_offset += trunc_count;
  d->count -= trunc_count;
  d->size = (d->count + increase_size + ~d->page_mask) & d->page_mask;

  /* Re-allocate the file size. */
  if (posix_fallocate(d->fd, d->file_offset, d->size) != 0) {
    error("Failed to pre-allocate the dump file.");
  }

  /* Re-map starting at the end of the file. */
  if ((d->data = mmap(NULL, d->size, PROT_WRITE, MAP_SHARED, d->fd,
                      d->file_offset)) == MAP_FAILED) {
    error("Failed to allocate map of size %zi bytes (%s).", d->size,
          strerror(errno));
  }
}

/**
 * @brief Flush the #dump to disk.
 */
void dump_sync(struct dump *d) {
  if (msync(d->data, d->count, MS_SYNC) != 0)
    error("Failed to sync memory-mapped data.");
}

/**
 * @brief Finalize the #dump.
 */
void dump_close(struct dump *d) {
  /* Unmap the data in memory. */
  if (munmap(d->data, d->count) != 0) {
    error("Failed to unmap dump data (%s).", strerror(errno));
  }

  /* Truncate the file to the correct length. */
  if (ftruncate(d->fd, d->file_offset + d->count) != 0) {
    error("Failed to truncate dump file (%s).", strerror(errno));
  }

  /* Close the memory-mapped file. */
  if (close(d->fd) != 0) error("Failed to close memory-mapped file.");
}

/**
 * @brief Initialize a file dump.
 *
 * @param d The #dump to initialize.
 * @param filename The fully qualified name of the file in which to dump,
 *                 note that it will be overwritten.
 * @param size The initial buffer size for this #dump.
 */
void dump_init(struct dump *d, const char *filename, size_t size) {

  /* Create the output file.
     The option O_RDWR seems to be required by mmap.
  */
  if ((d->fd = open(filename, O_CREAT | O_RDWR, 0660)) == -1) {
    error("Failed to create dump file '%s' (%s).", filename, strerror(errno));
  }

  /* Adjust the size to be at least the page size. */
  const size_t page_mask = ~(sysconf(_SC_PAGE_SIZE) - 1);
  size = (size + ~page_mask) & page_mask;

  /* Pre-allocate the file size. */
  if (posix_fallocate(d->fd, 0, size) != 0) {
    error("Failed to pre-allocate the dump file.");
  }

  /* Map memory to the created file. */
  if ((d->data = mmap(NULL, size, PROT_WRITE, MAP_SHARED, d->fd, 0)) ==
      MAP_FAILED) {
    error("Failed to allocate map of size %zi bytes (%s).", size,
          strerror(errno));
  }

  /* Init some counters. */
  d->size = size;
  d->count = 0;
  d->file_offset = 0;
  d->page_mask = page_mask;
}

/**
 * @brief Restart a file dump.
 *
 * @param d The #dump to restart.
 * @param filename The fully qualified name of the file in which to dump,
 *                 note that it will be overwritten.
 */
void dump_restart(struct dump *d, const char *filename) {
  /* Create the output file.
     The option O_RDWR seems to be required by mmap.
  */
  if ((d->fd = open(filename, O_RDWR, 0660)) == -1) {
    error("Failed to open dump file '%s' (%s).", filename, strerror(errno));
  }

  /* Adjust the size to be at least the page size. */
  const size_t page_mask = ~(sysconf(_SC_PAGE_SIZE) - 1);
  size_t size = (d->size + ~page_mask) & page_mask;

  /* Pre-allocate the file size. */
  if (posix_fallocate(d->fd, 0, size) != 0) {
    error("Failed to pre-allocate the dump file.");
  }

  /* Map memory to the created file. */
  if ((d->data = mmap(NULL, size, PROT_WRITE, MAP_SHARED, d->fd,
                      d->file_offset)) == MAP_FAILED) {
    error("Failed to allocate map of size %zi bytes (%s).", d->size,
          strerror(errno));
  }
}

#endif
