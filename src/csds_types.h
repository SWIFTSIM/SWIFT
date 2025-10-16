/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumniepfl.ch)
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
#ifndef SWIFT_CSDS_TYPES_H
#define SWIFT_CSDS_TYPES_H

/* Includes required for the structs below */
#include "error.h"
#include "io_properties.h"
#include "units.h"

/* Include the CSDS */
#include "csds/src/logfile_writer.h"

/* Forward declarations needed for the function pointers in csds_field */
struct part;
struct xpart;
struct gpart;
struct spart;
struct engine;

/**
 * @brief Enum defining the type of particle associated with a given CSDS mask.
 * The values should match those in part_type.h.
 */
enum mask_for_type {
  mask_for_gas = 0,
  mask_for_dark_matter = 1,
  /* Only need a single type of dm. */
  mask_for_sinks = 3,
  mask_for_stars = 4,
  mask_for_black_hole = 5,
  mask_for_timestep = -1,
} __attribute__((packed));

/**
 * @brief Structure describing a single CSDS data field.
 * This structure definition is required by struct csds_writer to determine
 * the size of the fixed_fields array elements.
 */
struct csds_field {
  /* Name of the field */
  char name[CSDS_STRING_SIZE];

  /* Mask value. */
  unsigned int mask;

  /* Type of particle (follow part_type.h and -1 for timestamp). */
  enum mask_for_type type;

  /* The offset of the field within the particle */
  size_t offset;

  /* The size of the field */
  size_t size;

  /* Do we use the xpart or the normal one? */
  int use_xpart;

  /* Conversion functions (NULL if none) */
  void *(*conversion_hydro)(const struct part *, const struct xpart *xp,
                            const struct engine *e, void *buffer);
  void *(*conversion_grav)(const struct gpart *, const struct engine *e,
                            void *buffer);
  void *(*conversion_stars)(const struct spart *, const struct engine *e,
                             void *buffer);
};

#endif /* SWIFT_CSDS_TYPES_H */
