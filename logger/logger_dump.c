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
#include "logger_dump.h"
#include "logger_reader.h"
#include "logger_io.h"

/**
 * @brief Initialize the #logger_dump.
 *
 * If required this function will also reverse the offsets.
 * @param dump The #logger_dump.
 * @param filename the dump's filename.
 * @param reader The #logger_reader.
 */
void logger_dump_init(
    struct logger_dump *dump, char *filename,
    struct logger_reader *reader) {

  /* Set the pointer to the reader. */
  dump->reader = reader;

  /* Open file, map it and get its size. */
  if (reader->verbose > 1)
    message("Mapping the dump file.");
  dump->dump.map = io_mmap_file(filename, &dump->dump.file_size);

  /* Read header. */
  if (reader->verbose > 1)
    message("Reading the header.");
  header_read(&dump->header, dump);

  /* Print the header. */
  if (reader->verbose > 0) {
    header_print(&dump->header);
  }

  /* Check if the offset are corrupted */
  if (header_are_offset_corrupted(&dump->header)) {
    error("The offsets have been corrupted");
  }

  /* Reverse offset direction */
  if (reader->verbose > 1)
    message("Checking if offsets need to be reversed.");
  if (header_are_offset_backward(&dump->header)) {
    logger_dump_reverse_offset(dump);
  }

  /* Initialize the time array */
  if (reader->verbose > 1)
    message("Reading the time stamps.");
  time_array_init(&dump->times, dump);

  if (reader->verbose > 0) {
    time_array_print(&dump->times);
  }

}

/**
 * @brief Free the allocated memory and unmap the file.
 *
 * @param dump The #logger_dump.
 */
void logger_dump_free(struct logger_dump *dump) {
  io_munmap_file(dump->dump.map, dump->dump.file_size);

}


/**
 * @brief Reverse offset in dump file
 *
 * @param dump The #logger_dump
 */
void logger_dump_reverse_offset(struct logger_dump *dump) {

  struct header *header = &dump->header;
  const struct logger_reader *reader = dump->reader;

  if (!header_are_offset_backward(header)) {
    error("The offset are already reversed.");
  }


#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (reader->verbose > 0) {
    message("Check offsets...\n");
  }
  size_t offset_debug = header->offset_first;
  while (offset_debug < dump->dump.file_size) {
    tools_check_offset(header, dump->dump.map, &offset_debug);
  }
  if (reader->verbose > 0) {
    message("Check done\n");
  }
#endif

  /* reverse header offset */
  header_change_offset_direction(header, logger_offset_corrupted);
  
  size_t offset = header->offset_first;

  /* reverse chunks */
  if (reader->verbose > 0) {
    message("Reversing offsets...\n");
  }
  while (offset < dump->dump.file_size) {
    tools_reverse_offset(header, dump->dump.map, &offset);
  }
  if (reader->verbose > 0) {
    message("Reversing done\n");
  }

  /* reverse header offset */
  header_change_offset_direction(header, logger_offset_forward);

#ifdef SWIFT_DEBUG_CHECKS
  /* check offset */
  if (reader->verbose > 0) {
    message("Check offsets...\n");
  }
  offset_debug = header->offset_first;
  while (offset_debug < dump->dump.file_size) {
    tools_check_offset(header, dump->dump.map, &offset_debug);
  }
  if (reader->verbose > 0) {
    message("Check done\n");
  }
#endif

}
