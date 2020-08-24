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

#include "logger_loader_io.h"
#include "logger_logfile.h"
#include "logger_reader.h"

/**
 * @brief Check if enough space is available and increase it if required.
 *
 * @param t The #time_array.
 */
void time_array_ensure_size(struct time_array *t) {
  /* Check if we still have some place. */
  if (t->size < t->capacity) return;

  /* Increase the size */
  t->capacity *= 2;

  /* Allocate the new array */
  struct time_record *tmp = malloc(sizeof(struct time_record) * t->capacity);
  if (tmp == NULL) error("Failed to allocate the time records.");

  /* Copy the memory */
  memcpy(tmp, t->records, sizeof(struct time_record) * t->size);

  /* Cleanup the memory */
  free(t->records);

  /* Set the pointer to the new array */
  t->records = tmp;
}

/**
 * @brief Add an element to the #time_array.
 *
 * @param t The #time_array.
 * @param int_time The time in integer.
 * @param time The time in double.
 * @param offset The offset of the record.
 */
void time_array_append(struct time_array *t, const integertime_t int_time,
                       const double time, const size_t offset) {

  /* Increase the available space if required */
  time_array_ensure_size(t);

  /* Copy the values */
  t->records[t->size].time = time;
  t->records[t->size].int_time = int_time;
  t->records[t->size].offset = offset;

  /* Increase the size used. */
  t->size += 1;
}

/**
 * @brief read a time record.
 *
 * @param int_time integer time read.
 * @param time time read.
 * @param reader The #logger_reader.
 * @param offset position in the file.
 *
 */
size_t time_read(integertime_t *int_time, double *time,
                 const struct logger_reader *reader, size_t offset) {

  /* Initialize variables. */
  const struct header *h = &reader->log.header;
  void *map = h->log->log.map;

  size_t mask = 0;
  size_t prev_offset = 0;
  *int_time = 0;
  *time = 0;

  /* read record header. */
  map =
      logger_loader_io_read_mask(h, (char *)map + offset, &mask, &prev_offset);

#ifdef SWIFT_DEBUG_CHECKS

  /* check if time mask is present in log file header. */
  int ind = header_get_field_index(h, "timestamp");
  if (ind == -1) error("File header does not contain a mask for time.");

  /* check if reading a time record. */
  if (h->masks[ind].mask != mask) error("Not a time record.");
#endif

  /* read the record. */
  map =
      logger_loader_io_read_data(map, sizeof(unsigned long long int), int_time);
  map = logger_loader_io_read_data(map, sizeof(double), time);

  return (char *)map - (char *)h->log->log.map;
}

/**
 * @brief get offset of first time record
 *
 * @param h file #header
 * @return offset of first time record
 *
 */
size_t time_offset_first_record(const struct header *h) {

  /* Initialize a few variables. */
  size_t offset = h->offset_first_record;
  void *map = h->log->log.map;

  /* Check that the first record is really a time record. */
  int i = header_get_field_index(h, "timestamp");

  if (i == -1) error("Time mask not present in the log file header.");

  size_t mask = 0;
  logger_loader_io_read_mask(h, (char *)map + offset, &mask, NULL);

  if (mask != h->masks[i].mask) error("Log file should begin by timestep.");

  return h->offset_first_record;
}

/**
 * @brief Initialize an empty time array.
 *
 * @param t #time_array to initialize.
 */
void time_array_init(struct time_array *t) {
  /* Allocate the arrays */
  t->records = malloc(sizeof(struct time_record) * LOGGER_TIME_INIT_SIZE);
  if (t->records == NULL) error("Failed to initialize the time records.");

  /* Initialize the sizes */
  t->size = 0;
  t->capacity = LOGGER_TIME_INIT_SIZE;
}

/**
 * @brief Initialize a time array from a file.
 *
 * @param t #time_array to initialize.
 * @param log The #logger_logfile.
 */
void time_array_populate(struct time_array *t, struct logger_logfile *log) {

  /* Initialize a few variables. */
  integertime_t int_time = 0;
  double time = 0;

  /* get file size. */
  size_t file_size = log->log.mmap_size;

  /* get first timestamp. */
  size_t offset = time_offset_first_record(&log->header);
  while (offset < file_size) {
    /* read current time record and store it. */
    size_t tmp_offset = offset;
    time_read(&int_time, &time, log->reader, tmp_offset);
    time_array_append(t, int_time, time, offset);

    /* get next record. */
    int test = tools_get_next_record(&log->header, log->log.map, &offset,
                                     log->log.mmap_size);
    if (test == -1) break;
  }
}

/**
 * @brief access the time of a given record (by its offset).
 *
 * @param t #time_array to access.
 * @param offset offset of the record.
 *
 * @return integer time of the record.
 */
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset) {
  size_t ind = time_array_get_index(t, offset);
  return t->records[ind].int_time;
}

/**
 * @brief access the time of a given record (by its offset).
 *
 * @param t #time_array to access.
 * @param offset offset of the record.
 *
 * @return time of the record.
 */
double time_array_get_time(const struct time_array *t, const size_t offset) {
  size_t ind = time_array_get_index(t, offset);
  return t->records[ind].time;
}

/**
 * @brief Find the index of the last time record written before a given offset.
 *
 * @param t #time_array to access.
 * @param offset offset of the record.
 *
 * @return The index of the last time record.
 */
size_t time_array_get_index(const struct time_array *t, const size_t offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!t) error("NULL pointer.");

  if (offset < t->records[0].offset || offset > t->records[t->size - 1].offset)
    error("Offset outside of range. %zi > %zi > %zi",
          t->records[t->size - 1].offset, offset, t->records[0].offset);
#endif

  /* right will contain the index at the end of the loop */
  size_t left = 0;
  size_t right = t->size - 1;

  /* Find the time_array with the correct offset through a bisection method. */
  // TODO use interpolation search (same for the other binary searches)
  while (left <= right) {
    size_t center = (left + right) / 2;
    const size_t offset_center = t->records[center].offset;

    if (offset > offset_center) {
      left = center + 1;
    } else if (offset < offset_center) {
      right = center - 1;
    } else {
      return center;
    }
  }

  /* Avoid the sentinel */
  if (right == t->size - 1) {
    right = right - 1;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (t->records[right].offset > offset ||
      t->records[right + 1].offset <= offset) {
    error("Found the wrong element");
  }

#endif

  return right;
}

/**
 * @brief Find the index of the last time record written before a given time.
 *
 * @param t #time_array to access.
 * @param time The time requested.
 *
 * @return The index of the last time record.
 */
size_t time_array_get_index_from_time(const struct time_array *t,
                                      const double time) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!t) error("NULL pointer.");

  if (time < t->records[0].time || time > t->records[t->size - 1].time)
    error("Time outside of range (%g > %g).", time,
          t->records[t->size - 1].time);
#endif

  /* right will contain the index at the end of the loop */
  size_t left = 0;
  size_t right = t->size - 1;

  /* Find the time_array with the correct time through a bisection method. */
  while (left <= right) {
    size_t center = (left + right) / 2;
    const double time_center = t->records[center].time;

    if (time > time_center) {
      left = center + 1;
    } else if (time < time_center) {
      right = center - 1;
    } else {
      return center;
    }
  }

  /* Avoid the sentinel */
  if (right == t->size - 1) {
    right = right - 1;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (t->records[right].time > time || t->records[right + 1].time <= time) {
    error("Found the wrong element");
  }

#endif

  return right;
}

/**
 * @brief free memory of a #time_array
 *
 * @param t #time_array to free
 */
void time_array_free(struct time_array *t) {
  /* Free the arrays */
  free(t->records);
  t->records = NULL;

  /* Reset the counters */
  t->size = 0;
  t->capacity = 0;
}

/**
 * @brief print a #time_array
 *
 * @param t #time_array to print
 */
void time_array_print(const struct time_array *t) {
#ifdef SWIFT_DEBUG_CHECKS
  const size_t threshold = 1000;
#else
  const size_t threshold = 5;
#endif

  size_t n = t->size;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lli (%g)", n, t->records[0].int_time,
         t->records[0].time);

  /* Loop over all elements. */
  for (size_t i = 1; i < n; i++) {
    /* Skip the times at the center of the array. */
    if (i < threshold || i > up_threshold)
      printf(", %zi (%g)", t->records[i].offset, t->records[i].time);

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
  const size_t threshold = 4;

  size_t n = t->size;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lu", n, t->records[0].offset);

  /* Loop over all elements. */
  for (size_t i = 1; i < n; i++) {
    /* Skip the offset in the middle of the array. */
    if (i < threshold || i > up_threshold)
      printf(", %lu", t->records[i].offset);

    if (i == threshold) printf(", ...");
  }

  printf("]\n");
}
