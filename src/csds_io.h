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

struct csds_field {
  /* Name of the field */
  char name[csds_string_length];

  /* Mask value. */
  unsigned int mask;

  /* Type of particle (follow part_type.h and -1 for timestamp). */
  enum mask_type type;

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

/**
 * @brief Define a field common to all the simulations
 * (e.g. timestamp or special flag).
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param size_field The size of the field in bytes (sizeof)
 */
#define csds_define_common_field(csds_field, field_name, size_field) \
  {                                                                  \
    csds_field.offset = 0;                                           \
    csds_field.size = size_field;                                    \
    if (strlen(field_name) >= csds_string_length)                    \
      error("Name %s too long", field_name);                         \
    strcpy(csds_field.name, field_name);                             \
    csds_field.mask = 0;                                             \
    csds_field.use_xpart = -1;                                       \
    csds_field.conversion_hydro = NULL;                              \
    csds_field.conversion_grav = NULL;                               \
    csds_field.conversion_stars = NULL;                              \
  }

/**
 * @brief Define a field that simply requires a memcpy (e.g. the IDs).
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param part The type of particles (e.g. struct part)
 * @param field The field to write (e.g. id)
 */
#define csds_define_standard_field(csds_field, field_name, part, field) \
  {                                                                     \
    csds_field.offset = offsetof(part, field);                          \
    part *tmp;                                                          \
    csds_field.size = sizeof(tmp->field);                               \
    if (strlen(field_name) >= csds_string_length)                       \
      error("Name %s too long", field_name);                            \
    strcpy(csds_field.name, field_name);                                \
    csds_field.mask = 0;                                                \
    csds_field.use_xpart = -1;                                          \
    csds_field.conversion_hydro = NULL;                                 \
    csds_field.conversion_grav = NULL;                                  \
    csds_field.conversion_stars = NULL;                                 \
  }

/**
 * @brief Same than the previous function but for the hydro.
 * This function specifies the particle (part vs xpart)
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param part The type of particles (e.g. struct part)
 * @param field The field to write (e.g. id)
 * @param saving_xpart Does the field belongs to the xpart?
 */
#define csds_define_hydro_standard_field(csds_field, field_name, part, field, \
                                         saving_xpart)                        \
  {                                                                           \
    csds_define_standard_field(csds_field, field_name, part, field);          \
    csds_field.use_xpart = saving_xpart;                                      \
  }

/**
 * @brief Define a field from a function.
 * This function should never be called by the user.
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param conversion_func The conversion function.
 * @param field_size The size of the field to write.
 * @param part_type The type of particles.
 */
#define csds_define_field_from_function_general(                    \
    csds_field, field_name, conversion_func, field_size, part_type) \
  {                                                                 \
    if (strlen(field_name) >= csds_string_length)                   \
      error("Name %s too long", field_name);                        \
    strcpy(csds_field.name, field_name);                            \
    csds_field.size = field_size;                                   \
    csds_field.mask = 0;                                            \
    csds_field.use_xpart = -1;                                      \
    csds_field.conversion_hydro = NULL;                             \
    csds_field.conversion_grav = NULL;                              \
    csds_field.conversion_stars = NULL;                             \
    csds_field.conversion_##part_type = conversion_func;            \
  }

/**
 * @brief Define a field from a function for hydro.
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param conversion_func The conversion function.
 * @param field_size The size of the field to write.
 */
#define csds_define_field_from_function_hydro(csds_field, field_name,     \
                                              conversion_func, size)      \
  {                                                                       \
    csds_define_field_from_function_general(csds_field, field_name,       \
                                            conversion_func, size, hydro) \
  }

/**
 * @brief Define a field from a function for stars.
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param conversion_func The conversion function.
 * @param field_size The size of the field to write.
 */
#define csds_define_field_from_function_stars(csds_field, field_name,     \
                                              conversion_func, size)      \
  {                                                                       \
    csds_define_field_from_function_general(csds_field, field_name,       \
                                            conversion_func, size, stars) \
  }

/**
 * @brief Define a field from a function for gravity.
 *
 * @param csds_field A pointer to the #csds_field to initialize.
 * @param field_name A string containing the name of the field.
 * @param conversion_func The conversion function.
 * @param field_size The size of the field to write.
 */
#define csds_define_field_from_function_gravity(csds_field, field_name,  \
                                                conversion_func, size)   \
  {                                                                      \
    csds_define_field_from_function_general(csds_field, field_name,      \
                                            conversion_func, size, grav) \
  }

void csds_write_description(struct csds_writer *log, struct engine *e);

#endif

#endif /* SWIFT_CSDS_IO_H */
