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
#ifndef SWIFT_CHEMISTRY_NONE_CHEMISTRY_LOGGER_H
#define SWIFT_CHEMISTRY_NONE_CHEMISTRY_LOGGER_H

#include "logger_io.h"

#ifdef WITH_LOGGER

/*
 * List of all possible mask.
 * Outside the module, only chemistry_logger_field_count is used.
 */
enum chemistry_logger_fields_part {
  chemistry_logger_field_part_count = 0,
};
enum chemistry_logger_fields_spart {
  chemistry_logger_field_spart_count = 0,
};

/* Name of each possible mask. */
extern const char
    *chemistry_logger_field_names_part[chemistry_logger_field_part_count];
extern const char
    *chemistry_logger_field_names_spart[chemistry_logger_field_spart_count];

/**
 * @brief Initialize the logger for the #part.
 *
 * WARNING: The order should be the same in all the functions and
 * #chemistry_logger_fields_part!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int chemistry_logger_writer_populate_mask_data_part(
    struct mask_data *mask_data) {
  return chemistry_logger_field_part_count;
}

/**
 * @brief Initialize the logger for the #spart.
 *
 * WARNING: The order should be the same in all the functions and
 * #chemistry_logger_fields_spart!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int chemistry_logger_writer_populate_mask_data_spart(
    struct mask_data *mask_data) {
  return chemistry_logger_field_spart_count;
}

/**
 * @brief Generates the mask and compute the size of the record for the #part.
 *
 * WARNING: The order should be the same in all the functions and
 * #chemistry_logger_fields_part!
 *
 * @param masks The list of masks (same order than in
 * #chemistry_logger_writer_populate_mask_data_part).
 * @param part The #part that will be written.
 * @param xpart The #xpart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void chemistry_logger_compute_size_and_mask_part(
    const struct mask_data *masks, const struct part *part,
    const struct xpart *xpart, const int write_all, size_t *buffer_size,
    unsigned int *mask) {}

/**
 * @brief Generates the mask and compute the size of the record for the #spart.
 *
 * WARNING: The order should be the same in all the functions and
 * #chemistry_logger_fields_spart!
 *
 * @param masks The list of masks (same order than in
 * #chemistry_logger_writer_populate_mask_data_spart).
 * @param spart The #spart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void chemistry_logger_compute_size_and_mask_spart(
    const struct mask_data *masks, const struct spart *spart,
    const int write_all, size_t *buffer_size, unsigned int *mask) {}

/**
 * @brief Write a #part to the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_logger_fields_part!
 *
 * @param masks The list of masks (same order than in
 * #chemistry_logger_writer_populate_mask_data_part).
 * @param p The #part to write.
 * @param xp The #xpart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *chemistry_logger_write_particle(
    const struct mask_data *mask_data, const struct part *p,
    const struct xpart *xp, unsigned int *mask, char *buff) {
  return buff;
}

/**
 * @brief Write a #spart to the logger.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_logger_fields_spart!
 *
 * @param masks The list of masks (same order than in
 * #chemistry_logger_writer_populate_mask_data_spart).
 * @param sp The #spart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *chemistry_logger_write_sparticle(
    const struct mask_data *mask_data, const struct spart *sp,
    unsigned int *mask, char *buff) {
  return buff;
}

#endif  // WITH_LOGGER
#endif  // SWIFT_CHEMISTRY_NONE_CHEMISTRY_LOGGER_H
