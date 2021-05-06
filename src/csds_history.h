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
#ifndef SWIFT_CSDS_HISTORY_H
#define SWIFT_CSDS_HISTORY_H

#include "../config.h"

/* Standard includes */
#include <stdint.h>

/* Local include */
#include "error.h"
#include "lock.h"
#include "part_type.h"

#if defined(WITH_CSDS)

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
struct csds_index_data {
  /* Id of the particle. */
  int64_t id;

  /* Offset of the particle in the file. */
  uint64_t offset;
};

/**
 * @brief Structure dealing with the changes in the number
 * of particles (e.g. creation, deletion, transformation).
 */
struct csds_history {

  /* Number of elements currently stored */
  uint64_t size;

  /* Size of the current buffer */
  size_t capacity;

  /* Buffer containing the particles */
  struct csds_index_data *data;

  /*! Spin lock for logging events. */
  swift_lock_type lock;
};

void csds_history_init(struct csds_history *hist);
void csds_history_reset(struct csds_history *hist);
void csds_history_free(struct csds_history *hist);
void csds_history_log(struct csds_history *hist, const long long id,
                      const uint64_t last_offset);
void csds_history_write(struct csds_history *hist, struct engine *e, FILE *f);

void csds_history_dump(const struct csds_history *hist, FILE *stream);
void csds_history_restore(struct csds_history *hist, FILE *stream);

#endif  // WITH_CSDS
#endif  // SWIFT_CSDS_HISTORY_H
