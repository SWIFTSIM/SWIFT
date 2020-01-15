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
#include "logger_loader_io.h"
#include "logger_tools.h"

/**
 * @brief get the size of a file.
 *
 * @param fd file id.
 *
 * @return file size.
 */
size_t logger_loader_io_get_file_size(int fd) {
  struct stat s;
  int status = fstat(fd, &s);
  if (status != 0) error("Unable to get file size (%s).", strerror(errno));
  return s.st_size;
}

/**
 * @brief Map a file.
 *
 * #logger_loader_io_munmap_file should be called to unmap the file.
 *
 * @param map The #mapped_file.
 * @param filename file to read.
 * @param read_only Open the file in read only mode?
 *
 */
void logger_loader_io_mmap_file(struct mapped_file *map, const char *filename,
                                int read_only) {
  /* open the file. */
  int fd;

  if (read_only)
    fd = open(filename, O_RDONLY);
  else
    fd = open(filename, O_RDWR);

  if (fd == -1)
    error("Unable to open file %s (%s).", filename, strerror(errno));

  /* get the file size. */
  map->mmap_size = logger_loader_io_get_file_size(fd);

  /* map the memory. */
  int mode = PROT_READ;
  if (!read_only) mode |= PROT_WRITE;

  map->map = mmap(NULL, map->mmap_size, mode, MAP_SHARED, fd, 0);
  if (map->map == MAP_FAILED)
    error("Failed to allocate map of size %zi bytes (%s).", map->mmap_size,
          strerror(errno));

  /* Close the file. */
  close(fd);
}

/**
 * @brief Unmap a file.
 *
 * @param map The #mapped_file.
 *
 */
void logger_loader_io_munmap_file(struct mapped_file *map) {
  /* unmap the file. */
  if (munmap(map->map, map->mmap_size) != 0) {
    error("Unable to unmap the file (%s).", strerror(errno));
  }

  /* Reset values */
  map->map = NULL;
  map->mmap_size = 0;
}
