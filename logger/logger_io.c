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
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "logger_header.h"
#include "logger_io.h"
#include "logger_tools.h"

/**
 * @brief get the size of a file
 *
 * @param fd file id
 * @param size out: file size
 *
 */
void io_get_file_size(int fd, size_t *size) {
  struct stat s;
  int status = fstat(fd, &s);
  if (status != 0) error("Unable to get file size (%s)", strerror(errno));
  *size = s.st_size;
}

/**
 * @brief Open a file and map it
 *
 * @param filename file to read
 * @param fd out: file id
 * @param map out: file mapping
 *
 */
void io_open_file(char *filename, int *fd, void **map) {
  /* open file */
  *fd = open(filename, O_RDWR);
  if (*fd == -1)
    error("Unable to open file %s (%s)", filename, strerror(errno));

  /* get file size */
  size_t size = 0;
  io_get_file_size(*fd, &size);

  /* map memory */
  *map = mmap(NULL, size, PROT_WRITE | PROT_READ, MAP_SHARED, *fd, 0);
  if (map == MAP_FAILED)
    error("Failed to allocate map of size %zi bytes. (%s)", size,
          strerror(errno));
}

/**
 * @brief Close a file and unmap it
 *
 * @param fd file id
 * @param map file mapping
 *
 */
void io_close_file(int *fd, void **map) {
  /* get file size */
  size_t size = 0;
  io_get_file_size(*fd, &size);

  /* unmap */
  if (munmap(*map, size) != 0) {
    error("Unable to unmap the file (%s)", strerror(errno));
  }

  close(*fd);
}
