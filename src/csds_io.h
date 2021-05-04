/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_CSDS_IO_H
#define SWIFT_CSDS_IO_H

/* Config parameters. */
#include "../config.h"

#ifdef WITH_CSDS

/* Includes. */
#include "engine.h"
#include "io_properties.h"
#include "part.h"
#include "units.h"

/* This enum defines the type of particle to use
   with a given mask.
   The values should be the same than in part_type.h. */
enum mask_type {
  mask_type_gas = 0,
  mask_type_dark_matter = 1,
  /* Only need a single type of dm. */
  mask_type_stars = 4,
  mask_type_black_hole = 5,
  mask_type_timestep = -1,
} __attribute__((packed));

struct mask_data {
  /* Number of bytes for a mask. */
  int size;

  /* Mask value. */
  unsigned int mask;

  /* Type of particle (follow part_type.h and -1 for timestamp). */
  enum mask_type type;

  /* Name of the mask. */
  char name[100];

  /* Variables used only for the reader. */
  struct {
    /* Index of the fields containing the first derivative (< 0 for none) */
    int first_deriv;

    /* Index of the fields containing the second derivative (< 0 for none) */
    int second_deriv;
  } reader;
};

void write_index_array(const struct engine* e, FILE* f, struct io_props* props,
                       size_t n_props, size_t N);
/**
 * @brief Initialize the mask_data with a given field.
 *
 * @param name The name of the field.
 * @param size The size of the field.
 *
 * @return The new mask_data.
 */
INLINE static struct mask_data csds_create_mask_entry(const char* name,
                                                      int size) {
  struct mask_data mask;
  /* Copy the fields */
  strcpy(mask.name, name);
  mask.size = size;
  mask.mask = 0;

  return mask;
}

/**
 * @brief Add a given field to the current mask.
 *
 * @param mask_data The mask_data corresponding to the field that we wish to
 * write.
 * @param buffer_size (in) The current size of the future buffer. (out) The
 * updated size.
 *
 * @return The mask of the current field.
 */
INLINE static size_t csds_add_field_to_mask(struct mask_data mask_data,
                                            size_t* buffer_size) {

  *buffer_size += mask_data.size;
  return mask_data.mask;
}

/**
 * @brief Check if a field should be written according to the mask set in
 * #csds_add_field_to_mask.
 *
 * @param mask_data The mask_data corresponding to the current field.
 * @param mask The mask used for the current record.
 */
INLINE static int csds_should_write_field(struct mask_data mask_data,
                                          unsigned int* mask) {

  const int test = mask_data.mask & *mask;
  if (test) {
    *mask &= ~mask_data.mask;
  }

  return test;
}

void csds_write_index_file(struct csds_writer* log, struct engine* e);
void csds_write_description(struct csds_writer* log, struct engine* e);

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param xparts The extra particle array.
 * @param list (out) The parameters to write.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static int hydro_write_index(
    const struct part* parts, const struct xpart* xparts,
    struct io_props* list) {

  /* List what we want to write */
  list[0] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           parts, id, "Field not used");
  list[1] =
      io_make_output_field("Offset", UINT64, 1, UNIT_CONV_NO_UNITS, 0.f, xparts,
                           csds_data.last_offset, "Field not used");

  return 2;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param gparts The gparticle array.
 * @param list (out) The parameters to write.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static int darkmatter_write_index(
    const struct gpart* gparts, struct io_props* list) {

  /* List what we want to write */
  list[0] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           gparts, id_or_neg_offset, "Field not used");
  list[1] =
      io_make_output_field("Offset", UINT64, 1, UNIT_CONV_NO_UNITS, 0.f, gparts,
                           csds_data.last_offset, "Field not used");

  return 2;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param sparts The sparticle array.
 * @param list (out) The parameters to write.
 *
 * In this version, we only want the ids and the offset.
 */
__attribute__((always_inline)) INLINE static int stars_write_index(
    const struct spart* sparts, struct io_props* list) {

  /* List what we want to write */
  list[0] =
      io_make_output_field("ParticleIDs", ULONGLONG, 1, UNIT_CONV_NO_UNITS, 0.f,
                           sparts, id, "Field not used");
  list[1] =
      io_make_output_field("Offset", UINT64, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
                           csds_data.last_offset, "Field not used");

  return 2;
}

#endif

#endif /* SWIFT_CSDS_IO_H */
