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
#include "logger_history.h"

/* Standard includes */
#include <string.h>

/* Local include */
#include "logger_io.h"
#include "part.h"

#if defined(WITH_LOGGER)

#define LOGGER_HISTORY_INIT_SIZE 1024

/**
 * @brief Initialize the structure for the first time.
 *
 * @param hist The #logger_history.
 */
void logger_history_init(struct logger_history *hist) {

  /* Set the counters to their initial value */
  hist->size = 0;
  hist->capacity = LOGGER_HISTORY_INIT_SIZE;
  lock_init(hist->lock);

  hist->data = (struct logger_index_data *)swift_malloc(
      "logger_history",
      sizeof(struct logger_index_data) * LOGGER_HISTORY_INIT_SIZE);
  if (hist->data == NULL) {
    error("Failed to allocate memory for the logger_history.");
  }
}

/**
 * @brief Reset the structure (for example just after a dump).
 *
 * @param hist The #logger_history.
 * @param params The #swift_params.
 * @param already_allocated Are the data already allocated? (Need to free it?)
 */
void logger_history_reset(struct logger_history *hist) {

  swift_free("logger_history", hist->data);

  logger_history_init(hist);
}

/**
 * @brief Free the structure (e.g. just before exiting).
 *
 * @param hist The #logger_history.
 */
void logger_history_free(struct logger_history *hist) {
  /* Set the counters to 0 */
  hist->size = 0;
  hist->capacity = 0;
  lock_destroy(hist->lock);

  /* Free the memory */
  if (hist->data != NULL) {
    swift_free("logger_history", hist->data);
    hist->data = NULL;
  }
}

/**
 * @brief Log a the particle information into the #logger_history.
 *
 * @param hist The #logger_history.
 * @param data The data from the particle.
 */
void logger_history_log(struct logger_history *hist, const long long id,
                        const uint64_t last_offset) {

#ifdef SWIFT_DEBUG_CHECKS
  if (id < 0) {
    error(
        "Negative ID for a particle. "
        "Are you trying to log a gpart linked to another type of particles?");
  }
#endif
  const struct logger_index_data data = {id, last_offset};

  /* Lock the history */
  lock_lock(hist->lock);

  /* Check if enough space is left */
  if (hist->size == hist->capacity) {
    /* Compute the previous amount of memory */
    const size_t memsize = sizeof(struct logger_index_data) * hist->capacity;

    /* Increase the capacity of the array */
    hist->capacity *= 2;

    /* Allocate the new array and copy the content of the previous one */
    struct logger_index_data *tmp =
        (struct logger_index_data *)swift_malloc("logger_history", 2 * memsize);

    memcpy(tmp, hist->data, memsize);

    /* Free the previous array and switch the pointers */
    swift_free("logger_history", hist->data);
    hist->data = tmp;
  }

  /* Save the new particle */
  hist->data[hist->size] = data;

  /* Increase the element counter */
  hist->size += 1;

  /* Unlock the history. */
  lock_unlock(hist->lock);
}

/**
 * @brief Write the history into an index file.
 *
 * @param hist The #logger_history.
 * @param e The #engine.
 * @param f The file where to write the history.
 */
void logger_history_write(struct logger_history *hist, struct engine *e,
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

  /* Reset the logger history */
  logger_history_reset(hist);
}

void logger_history_dump(const struct logger_history *hist, FILE *stream) {
  restart_write_blocks((void *)hist, sizeof(struct logger_history), 1, stream,
                       "logger_history", "logger_history");

  if (hist->size != 0)
    restart_write_blocks((void *)hist->data, sizeof(struct logger_index_data),
                         hist->size, stream, "logger_history_data",
                         "logger_history_data");
}

void logger_history_restore(struct logger_history *hist, FILE *stream) {
  restart_read_blocks((void *)hist, sizeof(struct logger_history), 1, stream,
                      NULL, "logger_history");

  hist->data = malloc(hist->capacity * sizeof(struct logger_index_data));
  if (hist->data == NULL) {
    error("Failed to allocate array for logger history");
  }

  if (hist->size != 0)
    restart_read_blocks((void *)hist->data, sizeof(struct logger_index_data),
                        hist->size, stream, NULL, "logger_history_data");
}

#endif  // WITH_LOGGER
