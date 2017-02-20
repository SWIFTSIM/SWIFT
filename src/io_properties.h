/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016  Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#ifndef SWIFT_IO_PROPERTIES_H
#define SWIFT_IO_PROPERTIES_H

/* Config parameters. */
#include "../config.h"

/**
 * @brief The two sorts of data present in the GADGET IC files: compulsory to
 * start a run or optional.
 */
enum DATA_IMPORTANCE { COMPULSORY = 1, OPTIONAL = 0, UNUSED = -1 };

/**
 * @brief The properties of a given dataset for i/o
 */
struct io_props {

  /* Name */
  char name[FIELD_BUFFER_SIZE];

  /* Type of the field */
  enum IO_DATA_TYPE type;

  /* Dimension (1D, 3D, ...) */
  int dimension;

  /* Is it compulsory ? (input only) */
  enum DATA_IMPORTANCE importance;

  /* Units of the quantity */
  enum UnitConversionFactor units;

  /* Pointer to the field of the first particle in the array */
  char* field;

  /* The size of the particles */
  size_t partSize;

  /* The particle arrays */
  struct part* parts;
  struct gpart* gparts;

  /* Conversion function for part */
  float (*convert_part)(struct engine*, struct part*);

  /* Conversion function for gpart */
  float (*convert_gpart)(struct engine*, struct gpart*);
};

/**
 * @brief Constructs an #io_props from its parameters
 */
#define io_make_input_field(name, type, dim, importance, units, part, field) \
  io_make_input_field_(name, type, dim, importance, units,                   \
                       (char*)(&(part[0]).field), sizeof(part[0]))

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param importance Is this dataset compulsory ?
 * @param units The units of the dataset
 * @param field Pointer to the field of the first particle
 * @param partSize The size in byte of the particle
 *
 * Do not call this function directly. Use the macro defined above.
 */
struct io_props io_make_input_field_(char name[FIELD_BUFFER_SIZE],
                                     enum IO_DATA_TYPE type, int dimension,
                                     enum DATA_IMPORTANCE importance,
                                     enum UnitConversionFactor units,
                                     char* field, size_t partSize) {
  struct io_props r;
  strcpy(r.name, name);
  r.type = type;
  r.dimension = dimension;
  r.importance = importance;
  r.units = units;
  r.field = field;
  r.partSize = partSize;
  r.parts = NULL;
  r.gparts = NULL;
  r.convert_part = NULL;
  r.convert_gpart = NULL;

  return r;
}

/**
 * @brief Constructs an #io_props from its parameters
 */
#define io_make_output_field(name, type, dim, units, part, field)          \
  io_make_output_field_(name, type, dim, units, (char*)(&(part[0]).field), \
                        sizeof(part[0]))

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param field Pointer to the field of the first particle
 * @param partSize The size in byte of the particle
 *
 * Do not call this function directly. Use the macro defined above.
 */
struct io_props io_make_output_field_(char name[FIELD_BUFFER_SIZE],
                                      enum IO_DATA_TYPE type, int dimension,
                                      enum UnitConversionFactor units,
                                      char* field, size_t partSize) {
  struct io_props r;
  strcpy(r.name, name);
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.field = field;
  r.partSize = partSize;
  r.parts = NULL;
  r.gparts = NULL;
  r.convert_part = NULL;
  r.convert_gpart = NULL;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_part(name, type, dim, units, part, field, \
                                          convert)                             \
  io_make_output_field_convert_part_(name, type, dim, units,                   \
                                     (char*)(&(part[0]).field),                \
                                     sizeof(part[0]), part, convert)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param field Pointer to the field of the first particle
 * @param partSize The size in byte of the particle
 * @param parts The particle array
 * @param functionPtr The function used to convert a particle to a float
 *
 * Do not call this function directly. Use the macro defined above.
 */
struct io_props io_make_output_field_convert_part_(
    char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum UnitConversionFactor units, char* field, size_t partSize,
    struct part* parts, float (*functionPtr)(struct engine*, struct part*)) {

  struct io_props r;
  strcpy(r.name, name);
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.field = field;
  r.partSize = partSize;
  r.parts = parts;
  r.gparts = NULL;
  r.convert_part = functionPtr;
  r.convert_gpart = NULL;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_gpart(name, type, dim, units, part, \
                                           field, convert)               \
  io_make_output_field_convert_gpart_(name, type, dim, units,            \
                                      (char*)(&(part[0]).field),         \
                                      sizeof(part[0]), gpart, convert)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param field Pointer to the field of the first particle
 * @param partSize The size in byte of the particle
 * @param gparts The particle array
 * @param functionPtr The function used to convert a g-particle to a float
 *
 * Do not call this function directly. Use the macro defined above.
 */
struct io_props io_make_output_field_convert_gpart_(
    char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum UnitConversionFactor units, char* field, size_t partSize,
    struct gpart* gparts, float (*functionPtr)(struct engine*, struct gpart*)) {

  struct io_props r;
  strcpy(r.name, name);
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.field = field;
  r.partSize = partSize;
  r.parts = NULL;
  r.gparts = gparts;
  r.convert_part = NULL;
  r.convert_gpart = functionPtr;

  return r;
}

#endif /* SWIFT_IO_PROPERTIES_H */
