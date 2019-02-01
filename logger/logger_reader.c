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
#include "logger_header.h"
#include "logger_io.h"

/**
 * @brief Reverse offset in dump file
 *
 * @param filename string filename of the dump file
 * @param verbose Verbose level
 */
void reverse_offset(char *filename, int verbose) {
  struct header h;

  /* open file */
  int fd;
  void *map;
  io_open_file(filename, &fd, &map);

  /* read header */
  header_read(&h, map);

  if (verbose > 0) {
    header_print(&h);
  }

  /* check offset direction */
  if (h.forward_offset) {
    error("Offset are already reversed");
  }

  /* compute file size */
  size_t sz;
  io_get_file_size(fd, &sz);

  size_t offset;

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (verbose > 0) {
    message("Check offsets...\n");
  }
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    message("Check done\n");
  }
#endif

  /* reverse header offset */
  header_change_offset_direction(&h, map);

  offset = h.offset_first;

  /* reverse chunks */
  if (verbose > 0) {
    message("Reversing offsets...\n");
  }
  while (offset < sz) {
    tools_reverse_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    message("Reversing done\n");
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (verbose > 0) {
    message("Check offsets...\n");
  }
  offset = h.offset_first;
  while (offset < sz) {
    tools_check_offset(&h, map, &offset);
  }
  if (verbose > 0) {
    message("Check done\n");
  }
#endif

  /* free internal variables */
  header_free(&h);

  io_close_file(&fd, &map);
}
