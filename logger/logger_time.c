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

/**
 * @brief read a time stamp
 *
 * @param timestamp timestamp read
 * @param time time read
 * @param h #header file structure
 * @param map file mapping
 * @param offset In: position in the file, Out: shifted by the timestamp
 */
void time_read(integertime_t *timestamp, double *time, const struct header *h,
               void *map, size_t *offset) {
  size_t mask = 0;
  size_t prev_offset = 0;
  *timestamp = 0;
  *time = 0;

  /* read chunck header */
  io_read_mask(h, map, offset, &mask, &prev_offset);

#ifdef SWIFT_DEBUG_CHECKS

  /* check if timestamp is present */
  int ind = header_get_field_index(h, "timestamp");
  if (ind == -1) error("Header does not contain a timestamp");

  /* check if timestamp */
  if (h->masks[ind].mask != mask) error("Not a timestamp");
#endif

  /* read data */
  // TODO
  io_read_data(map, sizeof(unsigned long long int), timestamp, offset);
  io_read_data(map, sizeof(double), time, offset);
}

/**
 * @brief get offset of first timestamp
 *
 * @param h file #header
 * @param map file mapping
 * @param offset out: offset of first timestamp
 *
 */
void time_first_timestamp(const struct header *h, void *map, size_t *offset) {
  *offset = h->offset_first;

  int i = header_get_field_index(h, "timestamp");

  if (i == -1) error("Time stamp not present in header");

  size_t tmp = *offset;
  size_t mask = 0;
  io_read_mask(h, map, offset, &mask, NULL);

  if (mask != h->masks[i].mask) error("Dump should begin by timestep");

  *offset = tmp;
}

/**
 * @brief Initialize a time array
 *
 * @param t #time_array to initialize
 * @param h #header file structure
 * @param map file mapping
 * @param fd file id
 */
void time_array_init(struct time_array *t, const struct header *h, void *map,
                     int fd) {

  t->next = NULL;
  t->prev = NULL;

  /* get first time stamp */
  size_t offset = 0;
  time_first_timestamp(h, map, &offset);

  integertime_t timestamp = 0;
  double time = 0;

  /* get file size */
  size_t file_size;
  io_get_file_size(fd, &file_size);

  while (offset < file_size) {

    /* read time */
    t->offset = offset;
    size_t tmp_offset = offset;
    time_read(&timestamp, &time, h, map, &tmp_offset);
    t->timestamp = timestamp;
    t->time = time;

    /* get next chunk */
    int test = tools_get_next_chunk(h, map, &offset, fd);
    if (test == -1) break;

    /* allocate next time_array */
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
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return integer time of the chunk
 */
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->timestamp;
}

/**
 * @brief access the time of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return time of the chunk
 */
double time_array_get_time(struct time_array *t, const size_t offset) {
  const struct time_array *tmp = time_array_get_time_array(t, offset);
  return tmp->time;
}

/**
 * @brief access the #time_array of a given chunk (by its offset)
 *
 * @param t #time_array to access
 * @param offset offset of the chunk
 *
 * @return pointer to the requested #time_array
 */
struct time_array *time_array_get_time_array(struct time_array *t,
                                             const size_t offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!t) return 0;
#endif
  while (t->next) {
    if (t->offset > offset) break;

    t = t->next;
  }

  if (t->prev == NULL) return t;

  struct time_array *tmp = t->prev;
  return tmp;
}

/**
 * @brief free memory of a #time_array
 *
 * @param t #time_array to free
 */
void time_array_free(struct time_array *t) {
  struct time_array *tmp;
  while (t->next) {
    tmp = t->next;
    free(t);
    t = tmp;
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
  size_t i = 0;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lli (%g)", n, t->timestamp, t->time);

  while (t->next) {
    i += 1;
    t = t->next;
    if (i < threshold || i > up_threshold)
      printf(", %lli (%g)", t->timestamp, t->time);

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
  size_t i = 0;
  size_t up_threshold = n - threshold;

  printf("Times (size %lu): [%lu", n, t->offset);

  while (t->next) {
    i += 1;
    t = t->next;
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
  while (t->next) {
    t = t->next;
    count += 1;
  }

  return count;
}
