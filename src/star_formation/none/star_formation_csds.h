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
#ifndef SWIFT_STAR_FORMATION_NONE_STAR_FORMATION_CSDS_H
#define SWIFT_STAR_FORMATION_NONE_STAR_FORMATION_CSDS_H

#include "csds_io.h"

#ifdef WITH_CSDS

/*
 * List of all possible mask.
 * Outside the module, only star_formation_csds_field_count is used.
 */
enum star_formation_csds_fields_spart {
  star_formation_csds_field_count = 0,
};

/* Name of each possible mask. */
extern const char
    *star_formation_csds_field_names[star_formation_csds_field_count];

/**
 * @brief Initialize the csds for the #spart.
 *
 * WARNING: The order should be the same in all the functions and
 * #star_formation_csds_fields_spart!
 *
 * @param mask_data Data for each type of mask.
 *
 * @return Number of masks used.
 */
INLINE static int star_formation_csds_writer_populate_mask_data(
    struct mask_data *mask_data) {
  return star_formation_csds_field_count;
}

/**
 * @brief Generates the mask and compute the size of the record for the #spart.
 *
 * WARNING: The order should be the same in all the functions and
 * #star_formation_csds_fields_spart!
 *
 * @param masks The list of masks (same order than in
 * #star_formation_csds_writer_populate_mask_data_spart).
 * @param spart The #spart that will be written.
 * @param write_all Are we forcing to write all the fields?
 *
 * @param buffer_size (out) The requested size for the buffer.
 * @param mask (out) The mask that will be written.
 */
INLINE static void star_formation_csds_compute_size_and_mask(
    const struct mask_data *masks, const struct spart *spart,
    const int write_all, size_t *buffer_size, unsigned int *mask) {}

/**
 * @brief Write a #spart to the csds.
 *
 * WARNING: The order should be the same in all the functions and
 * #hydro_csds_fields_spart!
 *
 * @param masks The list of masks (same order than in
 * #star_formation_csds_writer_populate_mask_data_spart).
 * @param sp The #spart to write.
 * @param mask The mask to use for this record.
 * @param buff The buffer where to write the particle.
 *
 * @return The buffer after the data.
 */
INLINE static char *star_formation_csds_write_sparticle(
    const struct mask_data *mask_data, const struct spart *sp,
    unsigned int *mask, char *buff) {
  return buff;
}

#endif  // WITH_CSDS
#endif  // SWIFT_STAR_FORMATION_NONE_STAR_FORMATION_CSDS_H
