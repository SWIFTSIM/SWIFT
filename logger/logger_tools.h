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
#ifndef __LOGGER_LOGGER_TOOLS_H__
#define __LOGGER_LOGGER_TOOLS_H__

#include "../config.h"

/* Swift include */
#include "../src/dimension.h"
#include "../src/inline.h"
#include "../src/logger.h"
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

struct header;
struct logger_reader;

#define error(s, ...)                                                        \
  ({                                                                         \
    fflush(stdout);                                                          \
    fprintf(stderr, "%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
            ##__VA_ARGS__);                                                  \
    abort();                                                                 \
  })

#ifdef SWIFT_DEBUG_CHECKS
#define message(s, ...)                                             \
  ({                                                                \
    printf("%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__, \
           ##__VA_ARGS__);                                          \
  })
#else
({								    \
    printf("%s(): " s "\n", __FUNCTION__,			    \
           ##__VA_ARGS__);                                          \
  })
#endif

int tools_get_next_record(const struct header *h, void *map, size_t *offset,
                         size_t file_size);
int _tools_get_next_record_backward(const struct header *h, void *map,
                                   size_t *offset, size_t file_size);
int _tools_get_next_record_forward(const struct header *h, void *map,
                                  size_t *offset);
size_t tools_reverse_offset(const struct header *h, void *map, size_t offset);
size_t tools_check_record_consistency(const struct logger_reader *reader, size_t offset);

#endif  //__LOGGER_LOGGER_TOOLS_H__
