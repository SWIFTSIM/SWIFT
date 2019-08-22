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

/* Local includes. */
#include "common_io.h"
#include "error.h"
#include "inline.h"
#include "part.h"

/* Standard includes. */
#include <string.h>

/**
 * @brief The two sorts of data present in the GADGET IC files: compulsory to
 * start a run or optional.
 */
enum DATA_IMPORTANCE { COMPULSORY = 1, OPTIONAL = 0, UNUSED = -1 };

/* Helper typedefs */
typedef void (*conversion_func_part_float)(const struct engine*,
                                           const struct part*,
                                           const struct xpart*, float*);
typedef void (*conversion_func_part_int)(const struct engine*,
                                         const struct part*,
                                         const struct xpart*, int*);
typedef void (*conversion_func_part_double)(const struct engine*,
                                            const struct part*,
                                            const struct xpart*, double*);
typedef void (*conversion_func_part_long_long)(const struct engine*,
                                               const struct part*,
                                               const struct xpart*, long long*);
typedef void (*conversion_func_gpart_float)(const struct engine*,
                                            const struct gpart*, float*);
typedef void (*conversion_func_gpart_int)(const struct engine*,
                                          const struct gpart*, int*);
typedef void (*conversion_func_gpart_double)(const struct engine*,
                                             const struct gpart*, double*);
typedef void (*conversion_func_gpart_long_long)(const struct engine*,
                                                const struct gpart*,
                                                long long*);
typedef void (*conversion_func_spart_float)(const struct engine*,
                                            const struct spart*, float*);
typedef void (*conversion_func_spart_int)(const struct engine*,
                                          const struct spart*, int*);
typedef void (*conversion_func_spart_double)(const struct engine*,
                                             const struct spart*, double*);
typedef void (*conversion_func_spart_long_long)(const struct engine*,
                                                const struct spart*,
                                                long long*);
typedef void (*conversion_func_bpart_float)(const struct engine*,
                                            const struct bpart*, float*);
typedef void (*conversion_func_bpart_int)(const struct engine*,
                                          const struct bpart*, int*);
typedef void (*conversion_func_bpart_double)(const struct engine*,
                                             const struct bpart*, double*);
typedef void (*conversion_func_bpart_long_long)(const struct engine*,
                                                const struct bpart*,
                                                long long*);

/**
 * @brief The properties of a given dataset for i/o
 */
struct io_props {

  /* Name */
  char name[FIELD_BUFFER_SIZE];

  /* Description of the variable to write to the field's meta-data */
  char description[DESCRIPTION_BUFFER_SIZE];

  /* Type of the field */
  enum IO_DATA_TYPE type;

  /* Dimension (1D, 3D, ...) */
  int dimension;

  /* Is it compulsory ? (input only) */
  enum DATA_IMPORTANCE importance;

  /* Units of the quantity */
  enum unit_conversion_factor units;

  /* Scale-factor exponent to apply for unit conversion to physical */
  float scale_factor_exponent;

  /* Pointer to the field of the first particle in the array */
  char* field;

  /* Pointer to the start of the temporary buffer used in i/o */
  char* start_temp_c;
  int* start_temp_i;
  float* start_temp_f;
  double* start_temp_d;
  long long* start_temp_l;

  /* Pointer to the engine */
  const struct engine* e;

  /* The size of the particles */
  size_t partSize;

  /* The particle arrays */
  const struct part* parts;
  const struct xpart* xparts;
  const struct gpart* gparts;
  const struct spart* sparts;
  const struct bpart* bparts;

  /* Are we converting? */
  int conversion;

  /* Conversion function for part */
  conversion_func_part_float convert_part_f;
  conversion_func_part_int convert_part_i;
  conversion_func_part_double convert_part_d;
  conversion_func_part_long_long convert_part_l;

  /* Conversion function for gpart */
  conversion_func_gpart_float convert_gpart_f;
  conversion_func_gpart_int convert_gpart_i;
  conversion_func_gpart_double convert_gpart_d;
  conversion_func_gpart_long_long convert_gpart_l;

  /* Conversion function for spart */
  conversion_func_spart_float convert_spart_f;
  conversion_func_spart_int convert_spart_i;
  conversion_func_spart_double convert_spart_d;
  conversion_func_spart_long_long convert_spart_l;

  /* Conversion function for bpart */
  conversion_func_bpart_float convert_bpart_f;
  conversion_func_bpart_int convert_bpart_i;
  conversion_func_bpart_double convert_bpart_d;
  conversion_func_bpart_long_long convert_bpart_l;
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
INLINE static struct io_props io_make_input_field_(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum DATA_IMPORTANCE importance, enum unit_conversion_factor units,
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
  r.xparts = NULL;
  r.gparts = NULL;
  r.sparts = NULL;
  r.bparts = NULL;
  r.conversion = 0;
  r.convert_part_f = NULL;
  r.convert_part_d = NULL;
  r.convert_part_l = NULL;
  r.convert_gpart_f = NULL;
  r.convert_gpart_d = NULL;
  r.convert_gpart_l = NULL;
  r.convert_spart_f = NULL;
  r.convert_spart_d = NULL;
  r.convert_spart_l = NULL;
  r.convert_bpart_f = NULL;
  r.convert_bpart_d = NULL;
  r.convert_bpart_l = NULL;

  return r;
}

/**
 * @brief Constructs an #io_props from its parameters
 */
#define io_make_output_field(name, type, dim, units, a_exponent, part, field, \
                             desc)                                            \
  io_make_output_field_(name, type, dim, units, a_exponent,                   \
                        (char*)(&(part[0]).field), sizeof(part[0]), desc)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param field Pointer to the field of the first particle
 * @param partSize The size in byte of the particle
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, char* field,
    size_t partSize, const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.field = field;
  r.partSize = partSize;
  r.conversion = 0;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_part(name, type, dim, units, a_exponent,  \
                                          part, xpart, convert, desc)          \
  io_make_output_field_convert_part_##type(name, type, dim, units, a_exponent, \
                                           sizeof(part[0]), part, xpart,       \
                                           convert, desc)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param partSize The size in byte of the particle
 * @param parts The particle array
 * @param xparts The xparticle array
 * @param functionPtr The function used to convert a particle to an int
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_part_INT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t partSize,
    const struct part* parts, const struct xpart* xparts,
    conversion_func_part_int functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = partSize;
  r.parts = parts;
  r.xparts = xparts;
  r.conversion = 1;
  r.convert_part_i = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param partSize The size in byte of the particle
 * @param parts The particle array
 * @param xparts The xparticle array
 * @param functionPtr The function used to convert a particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_part_FLOAT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t partSize,
    const struct part* parts, const struct xpart* xparts,
    conversion_func_part_float functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = partSize;
  r.parts = parts;
  r.xparts = xparts;
  r.conversion = 1;
  r.convert_part_f = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param partSize The size in byte of the particle
 * @param parts The particle array
 * @param xparts The xparticle array
 * @param functionPtr The function used to convert a particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_part_DOUBLE(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t partSize,
    const struct part* parts, const struct xpart* xparts,
    conversion_func_part_double functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = partSize;
  r.parts = parts;
  r.xparts = xparts;
  r.conversion = 1;
  r.convert_part_d = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param partSize The size in byte of the particle
 * @param parts The particle array
 * @param xparts The xparticle array
 * @param functionPtr The function used to convert a particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_part_LONGLONG(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t partSize,
    const struct part* parts, const struct xpart* xparts,
    conversion_func_part_long_long functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = partSize;
  r.parts = parts;
  r.xparts = xparts;
  r.conversion = 1;
  r.convert_part_l = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_gpart(name, type, dim, units, a_exponent, \
                                           gpart, convert, desc)               \
  io_make_output_field_convert_gpart_##type(name, type, dim, units,            \
                                            a_exponent, sizeof(gpart[0]),      \
                                            gpart, convert, desc)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param gpartSize The size in byte of the particle
 * @param gparts The particle array
 * @param functionPtr The function used to convert a g-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_gpart_INT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t gpartSize,
    const struct gpart* gparts, conversion_func_gpart_int functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = gpartSize;
  r.gparts = gparts;
  r.conversion = 1;
  r.convert_gpart_i = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param gpartSize The size in byte of the particle
 * @param gparts The particle array
 * @param functionPtr The function used to convert a g-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_gpart_FLOAT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t gpartSize,
    const struct gpart* gparts, conversion_func_gpart_float functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = gpartSize;
  r.gparts = gparts;
  r.conversion = 1;
  r.convert_gpart_f = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param gpartSize The size in byte of the particle
 * @param gparts The particle array
 * @param functionPtr The function used to convert a g-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_gpart_DOUBLE(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t gpartSize,
    const struct gpart* gparts, conversion_func_gpart_double functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = gpartSize;
  r.gparts = gparts;
  r.conversion = 1;
  r.convert_gpart_d = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param gpartSize The size in byte of the particle
 * @param gparts The particle array
 * @param functionPtr The function used to convert a g-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_gpart_LONGLONG(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t gpartSize,
    const struct gpart* gparts, conversion_func_gpart_long_long functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = gpartSize;
  r.gparts = gparts;
  r.conversion = 1;
  r.convert_gpart_l = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_spart(name, type, dim, units, a_exponent, \
                                           spart, convert, desc)               \
  io_make_output_field_convert_spart_##type(name, type, dim, units,            \
                                            a_exponent, sizeof(spart[0]),      \
                                            spart, convert, desc)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param spartSize The size in byte of the particle
 * @param sparts The particle array
 * @param functionPtr The function used to convert a s-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_spart_INT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t spartSize,
    const struct spart* sparts, conversion_func_spart_int functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = spartSize;
  r.sparts = sparts;
  r.conversion = 1;
  r.convert_spart_i = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param spartSize The size in byte of the particle
 * @param sparts The particle array
 * @param functionPtr The function used to convert a g-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_spart_FLOAT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t spartSize,
    const struct spart* sparts, conversion_func_spart_float functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = spartSize;
  r.sparts = sparts;
  r.conversion = 1;
  r.convert_spart_f = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param spartSize The size in byte of the particle
 * @param sparts The particle array
 * @param functionPtr The function used to convert a s-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_spart_DOUBLE(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t spartSize,
    const struct spart* sparts, conversion_func_spart_double functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = spartSize;
  r.sparts = sparts;
  r.conversion = 1;
  r.convert_spart_d = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param spartSize The size in byte of the particle
 * @param sparts The particle array
 * @param functionPtr The function used to convert a s-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_spart_LONGLONG(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t spartSize,
    const struct spart* sparts, conversion_func_spart_long_long functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = spartSize;
  r.sparts = sparts;
  r.conversion = 1;
  r.convert_spart_l = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_bpart(name, type, dim, units, a_exponent, \
                                           bpart, convert, desc)               \
  io_make_output_field_convert_bpart_##type(name, type, dim, units,            \
                                            a_exponent, sizeof(bpart[0]),      \
                                            bpart, convert, desc)

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param bpartSize The size in byte of the particle
 * @param bparts The particle array
 * @param functionPtr The function used to convert a b-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_bpart_INT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t bpartSize,
    const struct bpart* bparts, conversion_func_bpart_int functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = bpartSize;
  r.bparts = bparts;
  r.conversion = 1;
  r.convert_bpart_i = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param bpartSize The size in byte of the particle
 * @param bparts The particle array
 * @param functionPtr The function used to convert a g-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_bpart_FLOAT(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t bpartSize,
    const struct bpart* bparts, conversion_func_bpart_float functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = bpartSize;
  r.bparts = bparts;
  r.conversion = 1;
  r.convert_bpart_f = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param bpartSize The size in byte of the particle
 * @param bparts The particle array
 * @param functionPtr The function used to convert a s-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_bpart_DOUBLE(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t bpartSize,
    const struct bpart* bparts, conversion_func_bpart_double functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = bpartSize;
  r.bparts = bparts;
  r.conversion = 1;
  r.convert_bpart_d = functionPtr;

  return r;
}

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param bpartSize The size in byte of the particle
 * @param bparts The particle array
 * @param functionPtr The function used to convert a s-particle to a double
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_bpart_LONGLONG(
    const char name[FIELD_BUFFER_SIZE], enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t bpartSize,
    const struct bpart* bparts, conversion_func_bpart_long_long functionPtr,
    const char description[DESCRIPTION_BUFFER_SIZE]) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  strcpy(r.name, name);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    strcpy(r.description, description);
  }
  r.type = type;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = bpartSize;
  r.bparts = bparts;
  r.conversion = 1;
  r.convert_bpart_l = functionPtr;

  return r;
}

#endif /* SWIFT_IO_PROPERTIES_H */
