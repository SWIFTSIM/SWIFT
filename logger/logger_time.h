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
#ifndef LOGGER_LOGGER_TIMELINE_H
#define LOGGER_LOGGER_TIMELINE_H

#include "logger_header.h"
#include "logger_tools.h"

typedef int8_t timebin_t;
typedef long long integertime_t;

struct logger_reader;

#define LOGGER_TIME_INIT_SIZE 1024

/**
 * @brief This structure contains all the information present in a time record.
 */
struct time_record {
  /* Integertime of the records. */
  integertime_t int_time;

  /* Double time of the records. */
  double time;

  /* Offset in the file of the time records. */
  size_t offset;
};

/**
 * @brief This structure contains all the time record.
 *
 * In order to obtain easily the time step of a record,
 * this structure is required. It contains all the time step
 * with their integer time, double time and position in the file.
 *
 * This structure is initialized with #time_array_init and #time_array_populate,
 * and freed with #time_array_free.
 *
 * The time step of an offset can be obtained with
 * #time_array_get_integertime, #time_array_get_time and
 * #time_array_get_index.
 */
struct time_array {

  /* The complete list of time record */
  struct time_record *records;

  /* Number of element in the arrays. */
  size_t size;

  /* Maximum number of element available */
  size_t capacity;
};

void time_array_append(struct time_array *t, const integertime_t int_time,
                       const double time, const size_t offset);
size_t time_read(integertime_t *int_time, double *time,
                 const struct logger_reader *reader, size_t offset);

void time_array_init(struct time_array *t);
void time_array_populate(struct time_array *t, struct logger_logfile *log);

integertime_t time_array_get_integertime(struct time_array *t,
                                         const size_t offset);

double time_array_get_time(const struct time_array *t, const size_t offset);

size_t time_array_get_index(const struct time_array *t, const size_t offset);

void time_array_free(struct time_array *t);

void time_array_print(const struct time_array *t);

void time_array_print_offset(const struct time_array *t);

size_t time_offset_first_record(const struct header *h);

#endif  // LOGGER_LOGGER_TIMELINE_H
