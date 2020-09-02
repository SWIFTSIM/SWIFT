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

#include "logger_header.h"
#include "logger_logfile.h"
#include "logger_reader.h"
#include "swift.h"

int main(int argc, char *argv[]) {

  /*
    First generate the file.
  */

  message("Generating the dump.");
  /* Create required structures. */
  struct logger_writer log;
  struct swift_params params;
  char filename[200] = "testLogfileHeader.yml";

  /* Read parameters. */
  parser_read_file(filename, &params);

  /* Initialize the engine */
  struct engine e;
  e.policy =
      engine_policy_hydro | engine_policy_stars | engine_policy_self_gravity;

  /* Initialize the logger. */
  logger_init(&log, &e, &params);

  /* get dump filename. */
  char dump_filename[PARSER_MAX_LINE_SIZE];
  strcpy(dump_filename, log.base_name);
  strcat(dump_filename, "_0000.dump");

  /* Write file header. */
  logger_write_file_header(&log);

  /* Copy the masks */
  const int logger_count_mask_tmp = log.logger_count_mask;
  struct mask_data *mask_data = (struct mask_data *)malloc(
      log.logger_count_mask * sizeof(struct mask_data));
  memcpy(mask_data, log.logger_mask_data,
         log.logger_count_mask * sizeof(struct mask_data));

  /* clean memory. */
  logger_free(&log);
  /*
    Then read the file.
  */

  message("Reading the header.");
  /* Generate required structure for reading. */
  struct logger_reader reader;
  struct logger_logfile *logfile = &reader.log;
  logfile->reader = &reader;

  /* Set verbose level. */
  reader.verbose = 1;

  /* Read the header */
  logger_logfile_init_from_file(logfile, dump_filename, &reader,
                                /* only_header */ 1);
  /*
    Finally check everything.
  */

  struct header *h = &logfile->header;
  message("Checking versions.");
  assert(h->major_version == logger_major_version);
  assert(h->minor_version == logger_minor_version);

  message("Checking offset of first record");
  assert(h->offset_first_record == logfile->log.mmap_size);

  message("Checking number of masks");
  /* Get info about the masks */
  int count_mask = 0;
  size_t mask = 0;
  /* Skip the mask that appears multiple times */
  for (int i = 0; i < logger_count_mask_tmp; i++) {
    if (mask_data[i].mask & ~mask) {
      mask |= mask_data[i].mask;
      count_mask += 1;
    }
  }
  assert(h->masks_count == count_mask);

  message("Checking masks");
  mask = 0;
  int j = 0;
  for (int i = 0; i < count_mask; i++) {
    /* skip the duplicate masks */
    if (mask_data[i].mask & ~mask) {
      mask |= mask_data[i].mask;
      assert(mask_data[j].size == h->masks[i].size);
      assert(mask_data[j].mask == h->masks[i].mask);
      assert(strcmp(mask_data[j].name, h->masks[i].name) == 0);
      j += 1;
    }
  }

  message("Checking offset direction");
  assert(h->offset_direction == logger_offset_backward);

  free(mask_data);
  logger_logfile_free(logfile);
  return 0;
}
