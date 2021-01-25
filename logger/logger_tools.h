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
/**
 * @brief This file contains functions that help to navigate in the logs.
 */
#ifndef LOGGER_LOGGER_TOOLS_H
#define LOGGER_LOGGER_TOOLS_H

#include "../config.h"

/* Swift include */
#include "../src/dimension.h"
#include "../src/error.h"
#include "../src/inline.h"
#include "../src/logger.h"
#include "../src/logger_io.h"
#include "../src/part_type.h"

#ifdef HAVE_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#endif

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRING_SIZE 200

/* Define the size of all the fields. */
#define member_size(type, member) sizeof(((type *)0)->member)

struct header;
struct logger_reader;

/**
 * @brief Structure dealing with reading / interpolating a field.
 */
struct logger_field {
  /* Value of the field. */
  void *field;

  /* Value of its first derivative (NULL if not existing). */
  void *first_deriv;

  /* Value of its second derivative (NULL if not existing). */
  void *second_deriv;
};

/**
 * @brief Defines a type of modules (e.g. gravity, stars,
 * chemistry, star_formation, ...)
 */
enum module_field {
  field_module_default = 0,
  field_module_chemistry,
  field_module_star_formation,
};

/**
 * @brief Structure containing all the information
 * of a field.
 */
struct field_information {

  /* Module containing the field */
  enum module_field module;

  /* Name of the field */
  const char *name;

  /* Index of the field in the local array. */
  int local_index;

  /* Index of the field in the global array. */
  int global_index;

  /* Index of the first derivative in the local array. */
  int local_index_first;

  /* Index of the first derivative in the global array. */
  int global_index_first;

  /* Index of the second derivative in the local array. */
  int local_index_second;

  /* Index of the second derivative in the global array. */
  int global_index_second;
};

int tools_get_number_fields(enum part_type type);
void tools_get_list_fields(struct field_information *fields,
                           enum part_type type, const struct header *h);

int tools_get_next_record(const struct header *h, void *map, size_t *offset,
                          size_t file_size);
int _tools_get_next_record_backward(const struct header *h, void *map,
                                    size_t *offset, size_t file_size);
int _tools_get_next_record_forward(const struct header *h, void *map,
                                   size_t *offset);
size_t tools_reverse_offset(const struct header *h, void *map, size_t offset);
size_t tools_check_record_consistency(const struct logger_reader *reader,
                                      size_t offset);

void tools_print_progress(float percentage, int remaining_time,
                          const char *message);

#ifndef HAVE_PYTHON
#define error_python(s, ...) error(s, ##__VA_ARGS__);
#else
/**
 * @brief Print the python trace back
 */
__attribute__((always_inline)) INLINE static void logger_loader_print_traceback(
    void) {

  /* Import the traceback module */
  PyObject *pyth_module = PyImport_ImportModule("traceback");
  PyObject_CallMethod(pyth_module, "print_stack", "");

  Py_DECREF(pyth_module);
}

#define error_python(s, ...)         \
  ({                                 \
    logger_loader_print_traceback(); \
    error(s, ##__VA_ARGS__);         \
  })
#endif

#endif  // LOGGER_LOGGER_TOOLS_H
