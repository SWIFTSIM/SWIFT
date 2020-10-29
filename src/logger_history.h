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
#ifndef SWIFT_LOGGER_HISTORY_H
#define SWIFT_LOGGER_HISTORY_H

#include "../config.h"

/* Standard includes */
#include <stdint.h>

/* Local include */
#include "error.h"
#include "part_type.h"

#if defined(WITH_LOGGER)

/* Forward declaration */
struct xpart;
struct part;
struct gpart;
struct spart;
struct bpart;
struct engine;
struct swift_params;

/**
 * @brief Contains the information concerning
 * a particle for the index files.
 */
struct logger_index_data {
  /* Id of the particle. */
  int64_t id;

  /* Offset of the particle in the file. */
  uint64_t offset;
};

/**
 * @brief Structure dealing with the changes in the number
 * of particles (e.g. creation, deletion, transformation).
 */
struct logger_history {

  /* Number of elements currently stored */
  uint64_t size;

  /* Size of the current buffer */
  size_t capacity;

  /* Buffer containing the particles */
  struct logger_index_data *data;
};

void logger_history_init(struct logger_history *hist);
void logger_history_reset(struct logger_history *hist);
void logger_history_free(struct logger_history *hist);
void logger_history_log(struct logger_history *hist, const long long id,
                        const uint64_t last_offset);
void logger_history_write(struct logger_history *hist, struct engine *e,
                          FILE *f);

void logger_history_dump(const struct logger_history *hist, FILE *stream);
void logger_history_restore(struct logger_history *hist, FILE *stream);

#endif  // WITH_LOGGER
#endif  // SWIFT_LOGGER_HISTORY_H
