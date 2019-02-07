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
#ifndef __LOGGER_LOGGER_TIMELINE_H__
#define __LOGGER_LOGGER_TIMELINE_H__

#include "logger_header.h"
#include "logger_tools.h"

typedef char timebin_t;
typedef long long integertime_t;

/**
 * @brief This structure contains all the timestamp.
 *
 * In order to obtain easily the time step of a chunk,
 * this structure is required. It contains all the time step
 * with their timestamp and position in the file.
 *
 * This structure is initialized with @time_array_init and
 * freed with @time_array_free.
 *
 * The time step of an offset can be obtained with
 * @time_array_get_integertime, @time_array_get_time and
 * @time_array_get_time_array.
 *
 * The size of the time array can be accessed with
 * @time_array_count.
 */
struct time_array {
  /* Pointer to next element */
  void *next;

  /* Pointer to prev element */
  void *prev;

  /* Integertime of this timestamp */
  integertime_t timestamp;

  /* Double time of this timestamp */
  double time;

  /* Offset in the file of this timestamp */
  size_t offset;
};

void time_read(integertime_t *timestamp, double *time, const struct header *h,
               void *map, size_t *offset);
void time_array_init(struct time_array *t, const struct header *h, void *map,
                     int fd);
integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset);
double time_array_get_time(struct time_array *t, const size_t offset);
struct time_array *time_array_get_time_array(struct time_array *t,
                                             const size_t offset);
void time_array_free(struct time_array *t);
void time_array_print(const struct time_array *t);
void time_array_print_offset(const struct time_array *t);
size_t time_array_count(const struct time_array *t);
void time_first_timestamp(const struct header *h, void *map, size_t *offset);

#endif  // __LOGGER_LOGGER_TIMELINE_H__
