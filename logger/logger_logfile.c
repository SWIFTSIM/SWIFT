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
#include "logger_logfile.h"
#include "logger_reader.h"
#include "logger_io.h"

/**
 * @brief Initialize the #logger_logfile.
 *
 * If required this function will also reverse the offsets.
 * @param log The #logger_logfile.
 * @param filename the log's filename.
 * @param reader The #logger_reader.
 */
void logger_logfile_init(
    struct logger_logfile *log, char *filename,
    struct logger_reader *reader) {

  /* Set the pointer to the reader. */
  log->reader = reader;

  /* Set pointers to zero */
  time_array_init_to_zero(&log->times);

  /* Open file, map it and get its size. */
  if (reader->verbose > 1)
    message("Mapping the log file.");
  log->log.map = io_mmap_file(filename, &log->log.file_size);

  /* Read the header. */
  if (reader->verbose > 1)
    message("Reading the header.");
  header_read(&log->header, log);

  /* Print the header. */
  if (reader->verbose > 0) {
    header_print(&log->header);
  }

  /* Check if the offset are corrupted */
  if (header_is_corrupted(&log->header)) {
    error("The offsets have been corrupted");
  }

  /* Reverse the offsets direction */
  if (header_is_backward(&log->header)) {
    logger_logfile_reverse_offset(log);
  }

  /* Initialize the time array */
  if (reader->verbose > 1)
    message("Reading the time stamps.");
  time_array_init(&log->times, log);

  /* Print the time array */
  if (reader->verbose > 0) {
    time_array_print(&log->times);
  }

}

/**
 * @brief Free the allocated memory and unmap the file.
 *
 * @param log The #logger_logfile.
 */
void logger_logfile_free(struct logger_logfile *log) {
  io_munmap_file(log->log.map, log->log.file_size);

}


/**
 * @brief Reverse offset in log file
 *
 * @param log The #logger_logfile
 */
void logger_logfile_reverse_offset(struct logger_logfile *log) {

  /* Get pointers */
  struct header *header = &log->header;
  const struct logger_reader *reader = log->reader;

  /* Check if the offsets need to be reversed */
  if (!header_is_backward(header)) {
    error("The offsets are already reversed.");
  }


#ifdef SWIFT_DEBUG_CHECKS
  if (reader->verbose > 0) {
    message("Check record's headers...");
  }

  /* check that the record offset points to another record. */
  for(size_t offset_debug = header->offset_first_record;
      offset_debug < log->log.file_size;
      offset_debug = tools_check_record_consistency(reader, offset_debug)) {}

  if (reader->verbose > 0) {
    message("Record's headers are correct.");
  }
#endif

  message("WARNING: Modifying the logfile, do not kill the job!");

  /* Set the offset direction to a corrupted status. */
  header_change_offset_direction(header, logger_offset_corrupted);
  
  if (reader->verbose > 0) {
    message("Reversing offsets...");
  }

  /* reverse the record's offset. */
  for(size_t offset = header->offset_first_record;
      offset < log->log.file_size;
      offset = tools_reverse_offset(header, log->log.map, offset)) {}

  if (reader->verbose > 0) {
    message("Reversing done");
  }

  /* Now that the offset are effectively reversed, can set the direction to
     forward. */
  header_change_offset_direction(header, logger_offset_forward);

  message("WARNING: Modification done, you can now safely kill the job.");

#ifdef SWIFT_DEBUG_CHECKS
  if (reader->verbose > 0) {
    message("Check record's headers...");
  }

  /* check that the record offset points to another record. */
  for(size_t offset_debug = header->offset_first_record;
      offset_debug < log->log.file_size;
      offset_debug = tools_check_record_consistency(reader, offset_debug)) {}

  if (reader->verbose > 0) {
    message("Record's headers are correct.");
  }
#endif

}
