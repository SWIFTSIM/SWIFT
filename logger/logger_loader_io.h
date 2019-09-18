/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
/**
 * @file logger_loader_io.h
 * @brief This file contains basic IO function.
 */
#ifndef LOGGER_LOGGER_LOADER_IO_H
#define LOGGER_LOGGER_LOADER_IO_H

#include "logger_header.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

size_t logger_loader_io_get_file_size(int fd);
void *logger_loader_io_mmap_file(char *filename, size_t *file_size,
                                 int read_only);
void logger_loader_io_munmap_file(void *map, size_t file_size);

/**
 * @brief read a mask with its offset.
 *
 * @param h #header file structure.
 * @param data Pointer to the data to read.
 * @param mask (output) mask read from the data.
 * @param diff_offset (output) offset difference to previous/next corresponding
 * record.
 *
 * @return memory after the record header.
 */
__attribute__((always_inline)) INLINE static void *logger_loader_io_read_mask(
    const struct header *h, void *data, size_t *mask, size_t *diff_offset) {
  /* read mask */
  if (mask) {
    *mask = 0;
    memcpy(mask, data, LOGGER_MASK_SIZE);
  }
  data += LOGGER_MASK_SIZE;

  /* read offset */
  if (diff_offset) {
    *diff_offset = 0;
    memcpy(diff_offset, data, LOGGER_OFFSET_SIZE);
  }
  data += LOGGER_OFFSET_SIZE;

  return data;
}

/**
 * @brief read a single value from a file.
 *
 * @param data Pointer to the data to read.
 * @param size size of the data to read.
 * @param p pointer where to store the data.

 * @return memory after the data written.
 */
__attribute__((always_inline)) INLINE static void *logger_loader_io_read_data(
    void *data, const size_t size, void *p) {
  memcpy(p, data, size);
  return data + size;
};

/**
 * @brief write a single value in a file.
 *
 * @param data Pointer to the data to read.
 * @param size size of the data to write.
 * @param p pointer to the data.
 *
 * @return memory after the data written.
 */
__attribute__((always_inline)) INLINE static void *logger_loader_io_write_data(
    void *data, const size_t size, const void *p) {
  memcpy(data, p, size);

  return data + size;
};

#endif  // LOGGER_LOGGER_LOADER_IO_H
