/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_NONE_STARS_LOGGER_H
#define SWIFT_NONE_STARS_LOGGER_H

#ifdef WITH_LOGGER

#include "logger_io.h"

/**
 * List of all possible mask.
 * Outside the module, only stars_logger_field_count is used.
 */
enum stars_logger_fields {
  stars_logger_field_count,
};

/* Name of each possible mask. */
extern const char *stars_logger_field_names[stars_logger_field_count];

/**
 * @brief Initialize the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #stars_logger_fields!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int stars_logger_writer_populate_mask_data(
    struct mask_data *mask_data) {
  return 0;
}

/**
 * @brief Generates the mask and compute the size of the record.
 *
 * WARNING: The order should be the same in all the functions and
 * #stars_logger_fields!
 *
 * @param masks The list of masks (same order than in #stars_logger_init).
 * @param part The #spart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void stars_logger_compute_size_and_mask(
    const struct mask_data *masks, const struct spart *part,
    const int write_all, size_t *buffer_size, unsigned int *mask) {}

/**
 * @brief Write a particle to the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #stars_logger_fields!
 *
 * @param masks The list of masks (same order than in #stars_logger_init).
 * @param p The #spart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *stars_logger_write_particle(
    const struct mask_data *mask_data, const struct spart *p,
    unsigned int *mask, char *buff) {
  return NULL;
}

#endif /* WITH_LOGGER */
#endif /* SWIFT_NONE_STARS_LOGGER_H */
