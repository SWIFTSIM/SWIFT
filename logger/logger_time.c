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
#include "logger_time.h"
#include "logger_io.h"
#include "logger_logfile.h"
#include "logger_reader.h"

/**
 * @brief read a time record.
 *
 * @param int_time integer time read.
 * @param time time read.
 * @param reader The #logger_reader.
 * @param offset position in the file.
 *
 * @return position after the time record.
 */
size_t time_read(integertime_t *int_time, double *time, const struct logger_reader *reader,
               size_t offset) {

  /* Initialize variables. */
  const struct header *h = &reader->log.header;
  void *map = h->log->log.map;

  size_t mask = 0;
  size_t prev_offset = 0;
  *int_time = 0;
  *time = 0;

  /* read record header */
  offset = io_read_mask(h, map, offset, &mask, &prev_offset);

#ifdef SWIFT_DEBUG_CHECKS

  /* check if time mask is present in log file header. */
  int ind = header_get_field_index(h, "timestamp");
  if (ind == -1) error("File header does not contain a mask for time");

  /* check if reading a time record. */
  if (h->masks[ind].mask != mask) error("Not a time record");
#endif

  /* read the record */
  offset = io_read_data(map, sizeof(unsigned long long int), int_time, offset);
  offset = io_read_data(map, sizeof(double), time, offset);

  return offset;
}

/**
 * @brief get offset of first time record
 *
 * @param h file #header
 * @return offset of first time record
 *
 */
size_t time_offset_first_record(const struct header *h) {

  /* Initialize a few variables */
  size_t offset = h->offset_first_record;
  void *map = h->log->log.map;

  /* Check that the first record is really a time record. */
  int i = header_get_field_index(h, "timestamp");

  if (i == -1) error("Time mask not present in the log file header");

  size_t mask = 0;
  io_read_mask(h, map, offset, &mask, NULL);

  if (mask != h->masks[i].mask) error("Log file should begin by timestep");

  return h->offset_first_record;
}

/**
 * @brief Initialize an empty time array.
 *
 * @param t #time_array to initialize.
 */
void time_array_init_to_zero(struct time_array *t) {
  t->next = NULL;
  t->prev = NULL;
}

/**
 * @brief Initialize a time array.
 *
 * @param t #time_array to initialize.
 * @param log The #logger_logfile.
 */
void time_array_init(struct time_array *t, struct logger_logfile *log) {

  /* Initialize a few variables. */
  t->next = NULL;
  t->prev = NULL;

  integertime_t int_time = 0;
  double time = 0;

  /* get file size */
  size_t file_size = log->log.file_size;

  /* get first time stamp */
  size_t offset = time_offset_first_record(&log->header);
  while (offset < file_size) {
    /* read current time record and store it. */
    t->offset = offset;
    size_t tmp_offset = offset;
    tmp_offset = time_read(&int_time, &time, log->reader, tmp_offset);
    t->int_time = int_time;
    t->time = time;

    /* get next record */
    int test = tools_get_next_record(&log->header, log->log.map, &offset,
				    log->log.file_size);
    if (test == -1) break;

    /* allocate next time_array. */
    struct time_array *tmp = malloc(sizeof(struct time_array));
    tmp->prev = t;
    tmp->next = NULL;
    t->next = tmp;
    t = tmp;
  }

  /* free unused time_array */
  struct time_array *tmp = t->prev;
  tmp->next = NULL;
  free(t);
}

/**
 * @brief access the time of a given record (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the record
 *
 * @return integer time of the record
 */
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->int_time;
}

/**
 * @brief access the time of a given record (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the record
 *
 * @return time of the record
 */
double time_array_get_time(const struct time_array *t, const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->time;
}

/**
 * @brief access the #time_array of a given record (by its offset)
 *
 * @param t #time_array to access.
 * @param offset offset of the record.
 *
 * @return pointer to the requested #time_array
 */
struct time_array *time_array_get_time_array(const struct time_array *t,
                                             const size_t offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!t) error("NULL pointer");
#endif
  const struct time_array *tmp;
  /* Find the time_array with the correct offset */
  for(tmp = t; tmp->next && tmp->offset <= offset; tmp = tmp->next)
    {}

  /* If providing the offset of a time_array, need to give it back. */
  if (tmp->offset == offset) return (struct time_array *) tmp;

  return (struct time_array *) tmp->prev;
}

/**
 * @brief free memory of a #time_array
 *
 * @param t #time_array to free
 */
void time_array_free(struct time_array *t) {
  struct time_array *tmp;
  for(tmp = t; t->next; t = tmp) {
    tmp = t->next;
    free(t);
  }
}

/**
 * @brief print a #time_array
 *
 * @param t #time_array to print
 */
void time_array_print(const struct time_array *t) {
  size_t threshold = 4;

  size_t n = time_array_count(t);
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lli (%g)", n, t->int_time, t->time);

  /* Loop over all elements. */
  for(size_t i = 1; i < n; i++) {
    if (!t->next)
      error("Next pointer not initialized %zi", i);

    t = t->next;
    /* Skip the times at the center of the array. */
    if (i < threshold || i > up_threshold)
      printf(", %lli (%g)", t->int_time, t->time);

    if (i == threshold) printf(", ...");
  }

  printf("]\n");
}

/**
 * @brief print a #time_array (offset)
 *
 * @param t #time_array to print
 */
void time_array_print_offset(const struct time_array *t) {
  size_t threshold = 4;

  size_t n = time_array_count(t);
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lu", n, t->offset);

  /* Loop over all elements. */
  for(size_t i = 1; i < n; i++) {
    t = t->next;
    /* Skip the offset in the middle of the array. */
    if (i < threshold || i > up_threshold) printf(", %lu", t->offset);

    if (i == threshold) printf(", ...");
  }

  printf("]\n");
}

/**
 * @brief count number of element in #time_array
 *
 * @param t #time_array to count
 *
 * @return number of element
 */
size_t time_array_count(const struct time_array *t) {
  size_t count = 1;
  for(const struct time_array *tmp = t; tmp->next; tmp = tmp->next) {
    count += 1;
  }

  return count;
}
