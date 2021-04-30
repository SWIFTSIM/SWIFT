/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "csds_history.h"

/* Standard includes */
#include <string.h>

/* Local include */
#include "csds_io.h"
#include "part.h"

#if defined(WITH_CSDS)

#define CSDS_HISTORY_INIT_SIZE 1024

/**
 * @brief Initialize the structure for the first time.
 *
 * @param hist The #csds_history.
 */
void csds_history_init(struct csds_history *hist) {

  /* Set the counters to their initial value */
  hist->size = 0;
  hist->capacity = CSDS_HISTORY_INIT_SIZE;
  lock_init(&hist->lock);

  hist->data = (struct csds_index_data *)swift_malloc(
      "csds_history",
      sizeof(struct csds_index_data) * CSDS_HISTORY_INIT_SIZE);
  if (hist->data == NULL) {
    error("Failed to allocate memory for the csds_history.");
  }
}

/**
 * @brief Reset the structure (for example just after a dump).
 *
 * @param hist The #csds_history.
 * @param params The #swift_params.
 * @param already_allocated Are the data already allocated? (Need to free it?)
 */
void csds_history_reset(struct csds_history *hist) {

  swift_free("csds_history", hist->data);

  csds_history_init(hist);
}

/**
 * @brief Free the structure (e.g. just before exiting).
 *
 * @param hist The #csds_history.
 */
void csds_history_free(struct csds_history *hist) {
  /* Set the counters to 0 */
  hist->size = 0;
  hist->capacity = 0;
  if (lock_destroy(&hist->lock) != 0) error("Error destroying lock");

  /* Free the memory */
  if (hist->data != NULL) {
    swift_free("csds_history", hist->data);
    hist->data = NULL;
  }
}

/**
 * @brief Log a the particle information into the #csds_history.
 *
 * @param hist The #csds_history.
 * @param data The data from the particle.
 */
void csds_history_log(struct csds_history *hist, const long long id,
                        const uint64_t last_offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (id < 0) {
    error(
        "Negative ID for a particle. "
        "Are you trying to log a gpart linked to another type of particles?");
  }
#endif
  const struct csds_index_data data = {id, last_offset};

  /* Lock the history */
  lock_lock(&hist->lock);

  /* Check if enough space is left */
  if (hist->size == hist->capacity) {
    /* Compute the previous amount of memory */
    const size_t memsize = sizeof(struct csds_index_data) * hist->capacity;

    /* Increase the capacity of the array */
    hist->capacity *= 2;

    /* Allocate the new array and copy the content of the previous one */
    struct csds_index_data *tmp =
        (struct csds_index_data *)swift_malloc("csds_history", 2 * memsize);

    memcpy(tmp, hist->data, memsize);

    /* Free the previous array and switch the pointers */
    swift_free("csds_history", hist->data);
    hist->data = tmp;
  }

  /* Save the new particle */
  hist->data[hist->size] = data;

  /* Increase the element counter */
  hist->size += 1;

  /* Unlock the history. */
  if (lock_unlock(&hist->lock) != 0)
    error("Impossible to unlock CSDS history.");
}

/**
 * @brief Write the history into an index file.
 *
 * @param hist The #csds_history.
 * @param e The #engine.
 * @param f The file where to write the history.
 */
void csds_history_write(struct csds_history *hist, struct engine *e,
                          FILE *f) {
  /* Generate the structures for writing the index file */
  const int num_fields = 2;
  struct io_props list[2];
  list[0] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           hist->data, id, "Field not used");
  list[1] = io_make_output_field("Offset", UINT64, 1, UNIT_CONV_NO_UNITS, 0.f,
                                 hist->data, offset, "Field not used");

  write_index_array(e, f, list, num_fields, hist->size);

  /* Reset the csds history */
  csds_history_reset(hist);
}

void csds_history_dump(const struct csds_history *hist, FILE *stream) {
  restart_write_blocks((void *)hist, sizeof(struct csds_history), 1, stream,
                       "csds_history", "csds_history");

  if (hist->size != 0)
    restart_write_blocks((void *)hist->data, sizeof(struct csds_index_data),
                         hist->size, stream, "csds_history_data",
                         "csds_history_data");
}

void csds_history_restore(struct csds_history *hist, FILE *stream) {
  restart_read_blocks((void *)hist, sizeof(struct csds_history), 1, stream,
                      NULL, "csds_history");

  hist->data = malloc(hist->capacity * sizeof(struct csds_index_data));
  if (hist->data == NULL) {
    error("Failed to allocate array for CSDS history");
  }

  if (hist->size != 0)
    restart_read_blocks((void *)hist->data, sizeof(struct csds_index_data),
                        hist->size, stream, NULL, "csds_history_data");
}

#endif  // WITH_CSDS
