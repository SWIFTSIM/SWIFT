/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016  Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include <config.h>

/* Local includes. */
#include "common_io.h"
#include "error.h"
#include "inline.h"
#include "io_compression.h"
#include "part.h"
#include "units.h"

/* Standard includes. */
#include <string.h>

/**
 * @brief The two sorts of data present in the GADGET IC files: compulsory to
 * start a run or optional.
 */
enum DATA_IMPORTANCE { COMPULSORY = 1, OPTIONAL = 0, UNUSED = -1 };

/* Helper typedefs */
typedef void (*conversion_func_part_float)(const struct engine *,
                                           const struct part *,
                                           const struct xpart *, float *);
typedef void (*conversion_func_part_int)(const struct engine *,
                                         const struct part *,
                                         const struct xpart *, int *);
typedef void (*conversion_func_part_double)(const struct engine *,
                                            const struct part *,
                                            const struct xpart *, double *);
typedef void (*conversion_func_part_long_long)(const struct engine *,
                                               const struct part *,
                                               const struct xpart *,
                                               long long *);
typedef void (*conversion_func_gpart_float)(const struct engine *,
                                            const struct gpart *, float *);
typedef void (*conversion_func_gpart_int)(const struct engine *,
                                          const struct gpart *, int *);
typedef void (*conversion_func_gpart_double)(const struct engine *,
                                             const struct gpart *, double *);
typedef void (*conversion_func_gpart_long_long)(const struct engine *,
                                                const struct gpart *,
                                                long long *);
typedef void (*conversion_func_spart_float)(const struct engine *,
                                            const struct spart *, float *);
typedef void (*conversion_func_spart_int)(const struct engine *,
                                          const struct spart *, int *);
typedef void (*conversion_func_spart_double)(const struct engine *,
                                             const struct spart *, double *);
typedef void (*conversion_func_spart_long_long)(const struct engine *,
                                                const struct spart *,
                                                long long *);
typedef void (*conversion_func_bpart_float)(const struct engine *,
                                            const struct bpart *, float *);
typedef void (*conversion_func_bpart_int)(const struct engine *,
                                          const struct bpart *, int *);
typedef void (*conversion_func_bpart_double)(const struct engine *,
                                             const struct bpart *, double *);
typedef void (*conversion_func_bpart_long_long)(const struct engine *,
                                                const struct bpart *,
                                                long long *);
typedef void (*conversion_func_sink_float)(const struct engine *,
                                           const struct sink *, float *);
typedef void (*conversion_func_sink_int)(const struct engine *,
                                         const struct sink *, int *);
typedef void (*conversion_func_sink_double)(const struct engine *,
                                            const struct sink *, double *);
typedef void (*conversion_func_sink_long_long)(const struct engine *,
                                               const struct sink *,
                                               long long *);

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

  /* Has this entry been filled */
  int is_used;

  /* Is the entry actually written in the physical frame? */
  int is_physical;

  /* Is the entry not convertible to comoving (only meaningful if the entry is
   * physical) */
  int is_convertible_to_comoving;

  /* Is it compulsory ? (input only) */
  enum DATA_IMPORTANCE importance;

  /* Units of the quantity */
  enum unit_conversion_factor units;

  /* Default value to apply for optional fields when not found in the ICs */
  float default_value;

  /* Scale-factor exponent to apply for unit conversion to physical */
  float scale_factor_exponent;

  /* Pointer to the field of the first particle in the array */
  char *field;

  /* Lossy compression scheme to use for this field */
  enum lossy_compression_schemes lossy_compression;

  /* Pointer to the start of the temporary buffer used in i/o */
  char *start_temp_c;
  int *start_temp_i;
  float *start_temp_f;
  double *start_temp_d;
  long long *start_temp_l;

  /* Pointer to the engine */
  const struct engine *e;

  /* The size of the particles */
  size_t partSize;

  /* The particle arrays */
  const struct part *parts;
  const struct xpart *xparts;
  const struct gpart *gparts;
  const struct spart *sparts;
  const struct bpart *bparts;
  const struct sink *sinks;

  /* Are we converting? */
  int conversion;

  union {

    void *ptr_func;

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

    /* Conversion function for sink */
    conversion_func_sink_float convert_sink_f;
    conversion_func_sink_int convert_sink_i;
    conversion_func_sink_double convert_sink_d;
    conversion_func_sink_long_long convert_sink_l;
  };
};

/**
 * @brief Copies a string safely (avoids buffer overrun).
 *
 * @param dst Pointer to the destination array where the content is to be
 * copied.
 * @param src String to copy.
 * @param dst_len Length of the destination array.
 */
INLINE static void safe_strcpy(char *restrict dst, const char *restrict src,
                               size_t dst_len) {
  strncpy(dst, src, dst_len - 1);
  dst[dst_len - 1] = '\0';
}

/**
 * @brief Constructs an #io_props from its parameters
 *
 * @param name The name of the field in the ICs.
 * @param type The data type.
 * @param dim The dimensionality of the field.
 * @param importance Is this field compulsory or optional?
 * @param units The units used for this field.
 * @param part Pointer to the particle array where to write.
 * @param field Name of the field in the particle structure to write to.
 */
#define io_make_input_field(name, type, dim, importance, units, part, field) \
  io_make_input_field_(name, type, dim, importance, units,                   \
                       (char *)(&(part[0]).field), sizeof(part[0]), 0.)

/**
 * @brief Constructs an #io_props from its parameters with a user-defined
 * default value to use for optional fields.
 *
 * @param name The name of the field in the ICs.
 * @param type The data type.
 * @param dim The dimensionality of the field.
 * @param importance Is this field compulsory or optional?
 * @param units The units used for this field.
 * @param part Pointer to the particle array where to write.
 * @param field Name of the field in the particle structure to write to.
 * @param def The value to use as a default if the field is optional and not
 * found in the ICs.
 */
#define io_make_input_field_default(name, type, dim, importance, units, part, \
                                    field, def)                               \
  io_make_input_field_(name, type, dim, importance, units,                    \
                       (char *)(&(part[0]).field), sizeof(part[0]), def)

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
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum DATA_IMPORTANCE importance, enum unit_conversion_factor units,
    char *field, size_t partSize, const float default_value) {
  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  r.type = type;
  r.is_used = 1;
  r.dimension = dimension;
  r.importance = importance;
  r.units = units;
  r.field = field;
  r.partSize = partSize;
  r.default_value = default_value;

  if (default_value != 0.f && importance != OPTIONAL)
    error("Cannot set a non-zero default value for a compulsory field!");
  if (default_value != 0.f && type != FLOAT)
    error(
        "Can only set non-zero default value for a field using a FLOAT type!");

  return r;
}

/**
 * @brief Constructs an #io_props from its parameters
 */
#define io_make_output_field(name, type, dim, units, a_exponent, part, field, \
                             desc)                                            \
  io_make_output_field_(name, type, dim, units, a_exponent,                   \
                        (char *)(&(part[0]).field), sizeof(part[0]), desc,    \
                        /*physical=*/0, /*convertible_to_physical=*/1);

/**
 * @brief Constructs an #io_props from its parameters for a comoving quantity
 *
 * An alias of io_make_output_field().
 */
#define io_make_comoving_output_field(name, type, dim, units, a_exponent, \
                                      part, field, desc)                  \
  io_make_output_field(name, type, dim, units, a_exponent, part, field, desc)

/**
 * @brief Constructs an #io_props from its parameters for a physical quantity
 */
#define io_make_physical_output_field(name, type, dim, units, a_exponent,  \
                                      part, field, convertible, desc)      \
  io_make_output_field_(name, type, dim, units, a_exponent,                \
                        (char *)(&(part[0]).field), sizeof(part[0]), desc, \
                        /*physical=*/1,                                    \
                        /*convertible_to_physical=*/convertible);

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
 * @param is_physical Is the quantity written in the physical frame?
 * @param is_convertible_to_comoving Is the quantity convertible to comoving?
 *
 * The last argument only makes sense if the quantity is written to
 * in physical.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, char *field,
    size_t partSize, const char *description, const int is_physical,
    const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
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
#define io_make_output_field_convert_part(name, type, dim, units, a_exponent, \
                                          part, xpart, convert, desc)         \
  io_make_output_field_convert_part_(                                         \
      name, type, dim, units, a_exponent, sizeof(part[0]), part, xpart,       \
      convert, desc, /*physical=*/0, /*convertible_to_physical=*/1);

#define io_make_comoving_output_field_convert_part(                           \
    name, type, dim, units, a_exponent, part, xpart, convert, desc)           \
  io_make_output_field_convert_part(name, type, dim, units, a_exponent, part, \
                                    xpart, convert, desc)

#define io_make_physical_output_field_convert_part(name, type, dim, units,     \
                                                   a_exponent, part, xpart,    \
                                                   convertible, convert, desc) \
  io_make_output_field_convert_part_(name, type, dim, units, a_exponent,       \
                                     sizeof(part[0]), part, xpart, convert,    \
                                     desc, /*physical=*/1, convertible);

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
INLINE static struct io_props io_make_output_field_convert_part_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t partSize,
    const struct part *parts, const struct xpart *xparts, void *functionPtr,
    const char *description, const int is_physical,
    const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = partSize;
  r.parts = parts;
  r.xparts = xparts;
  r.conversion = 1;
  r.ptr_func = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_gpart(name, type, dim, units, a_exponent, \
                                           gpart, convert, desc)               \
  io_make_output_field_convert_gpart_(                                         \
      name, type, dim, units, a_exponent, sizeof(gpart[0]), gpart, convert,    \
      desc, /*physical=*/0, /*convertible_to_physical=*/1)

#define io_make_comoving_output_field_convert_gpart(                     \
    name, type, dim, units, a_exponent, gpart, convert, desc)            \
  io_make_output_field_convert_gpart(name, type, dim, units, a_exponent, \
                                     gpart, convert, desc)

#define io_make_physical_output_field_convert_gpart(                          \
    name, type, dim, units, a_exponent, gpart, convertible, convert, desc)    \
  io_make_output_field_convert_gpart_(name, type, dim, units, a_exponent,     \
                                      sizeof(gpart[0]), gpart, convert, desc, \
                                      /*physical=*/1, convertible);

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
INLINE static struct io_props io_make_output_field_convert_gpart_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t gpartSize,
    const struct gpart *gparts, void *functionPtr, const char *description,
    const int is_physical, const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = gpartSize;
  r.gparts = gparts;
  r.conversion = 1;
  r.ptr_func = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_spart(name, type, dim, units, a_exponent, \
                                           spart, convert, desc)               \
  io_make_output_field_convert_spart_(                                         \
      name, type, dim, units, a_exponent, sizeof(spart[0]), spart, convert,    \
      desc, /*physical=*/0, /*convertible_to_physical=*/1)

#define io_make_comoving_output_field_convert_spart(                     \
    name, type, dim, units, a_exponent, spart, convert, desc)            \
  io_make_output_field_convert_spart(name, type, dim, units, a_exponent, \
                                     spart, convert, desc)

#define io_make_physical_output_field_convert_spart(                          \
    name, type, dim, units, a_exponent, spart, convertible, convert, desc)    \
  io_make_output_field_convert_spart_(name, type, dim, units, a_exponent,     \
                                      sizeof(spart[0]), spart, convert, desc, \
                                      /*physical=*/1, convertible);

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
INLINE static struct io_props io_make_output_field_convert_spart_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t spartSize,
    const struct spart *sparts, void *functionPtr, const char *description,
    const int is_physical, const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = spartSize;
  r.sparts = sparts;
  r.conversion = 1;
  r.ptr_func = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_bpart(name, type, dim, units, a_exponent, \
                                           bpart, convert, desc)               \
  io_make_output_field_convert_bpart_(                                         \
      name, type, dim, units, a_exponent, sizeof(bpart[0]), bpart, convert,    \
      desc, /*physical=*/0, /*convertible_to_physical=*/1)

#define io_make_comoving_output_field_convert_bpart(                     \
    name, type, dim, units, a_exponent, bpart, convert, desc)            \
  io_make_output_field_convert_bpart(name, type, dim, units, a_exponent, \
                                     bpart, convert, desc)

#define io_make_physical_output_field_convert_bpart(                          \
    name, type, dim, units, a_exponent, bpart, convertible, convert, desc)    \
  io_make_output_field_convert_bpart_(name, type, dim, units, a_exponent,     \
                                      sizeof(bpart[0]), bpart, convert, desc, \
                                      /*physical=*/1, convertible);

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
INLINE static struct io_props io_make_output_field_convert_bpart_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t bpartSize,
    const struct bpart *bparts, void *functionPtr, const char *description,
    const int is_physical, const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = bpartSize;
  r.bparts = bparts;
  r.conversion = 1;
  r.ptr_func = functionPtr;

  return r;
}

/**
 * @brief Constructs an #io_props (with conversion) from its parameters
 */
#define io_make_output_field_convert_sink(name, type, dim, units, a_exponent, \
                                          sink, convert, desc)                \
  io_make_output_field_convert_sink_(                                         \
      name, type, dim, units, a_exponent, sizeof(sink[0]), sink, convert,     \
      desc, /*physical=*/0, /*convertible_to_physical=*/1)

#define io_make_comoving_output_field_convert_sink(                           \
    name, type, dim, units, a_exponent, ink, convert, desc)                   \
  io_make_output_field_convert_sink(name, type, dim, units, a_exponent, sink, \
                                    convert, desc)

#define io_make_physical_output_field_convert_sink(                        \
    name, type, dim, units, a_exponent, sink, convertible, convert, desc)  \
  io_make_output_field_convert_sink_(name, type, dim, units, a_exponent,   \
                                     sizeof(sink[0]), sink, convert, desc, \
                                     /*physical=*/1, convertible);

/**
 * @brief Construct an #io_props from its parameters
 *
 * @param name Name of the field to read
 * @param type The type of the data
 * @param dimension Dataset dimension (1D, 3D, ...)
 * @param units The units of the dataset
 * @param a_exponent Exponent of the scale-factor to convert to physical units.
 * @param sinkSize The size in byte of the particle
 * @param sinks The particle array
 * @param functionPtr The function used to convert a sink-particle to a float
 * @param description Description of the field added to the meta-data.
 *
 * Do not call this function directly. Use the macro defined above.
 */
INLINE static struct io_props io_make_output_field_convert_sink_(
    const char *name, enum IO_DATA_TYPE type, int dimension,
    enum unit_conversion_factor units, float a_exponent, size_t sinkSize,
    const struct sink *sinks, void *functionPtr, const char *description,
    const int is_physical, const int is_convertible_to_comoving) {

  struct io_props r;
  bzero(&r, sizeof(struct io_props));

  safe_strcpy(r.name, name, FIELD_BUFFER_SIZE);
  if (strlen(description) == 0) {
    sprintf(r.description, "No description given");
  } else {
    safe_strcpy(r.description, description, DESCRIPTION_BUFFER_SIZE);
  }
  r.type = type;
  r.is_used = 1;
  r.is_physical = is_physical;
  r.is_convertible_to_comoving = is_convertible_to_comoving;
  r.dimension = dimension;
  r.importance = UNUSED;
  r.units = units;
  r.scale_factor_exponent = a_exponent;
  r.partSize = sinkSize;
  r.sinks = sinks;
  r.conversion = 1;
  r.ptr_func = functionPtr;

  return r;
}

#endif /* SWIFT_IO_PROPERTIES_H */
