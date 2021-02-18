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
#ifndef SWIFT_LOGGER_PYTHON_TOOLS_H
#define SWIFT_LOGGER_PYTHON_TOOLS_H

#include "../config.h"

/* Local includes. */
#include "logger_tools.h"

#if HAVE_PYTHON
/* Include numpy */
#include <numpy/arrayobject.h>

#define CUSTOM_NPY_TYPE -1

/* Structure that allows the user to easily define
   subfield in a numpy array */
struct logger_python_subfield {
  /* Name of the subfield */
  char name[STRING_SIZE];

  /* Offset in the bytes */
  size_t offset;

  /* Numpy typenum of the subfield (e.g. NPY_FLOAT32, NPY_INT32). */
  int typenum;

  /* Number of dimensions for the subfield */
  int dimension;
};

/* Structure that allows the user to easily define the
   fields in a numpy array. */
struct logger_python_field {
  /* Dimension of the field (e.g. 1 for density and 3 for coordinates). */
  int dimension;

  /* Numpy typenum of the array (e.g. NPY_FLOAT32, NPY_INT32 or
   * CUSTOM_NPY_TYPE). */
  int typenum;

  /* Number of subfields registred */
  int subfields_registred;

  /* List of subtypes for custom fields */
  struct logger_python_subfield *subfields;
};

/**
 * @brief Generate a #logger_python_field structure.
 *
 * @param dimension The number of dimension for the field.
 * @param type The numpy data type (e.g. NPY_FLOAT32, NPY_DOUBLE, NPY_LONGLONG,
 * ...)
 *
 * @return The initialized structure.
 */
__attribute__((always_inline)) INLINE static struct logger_python_field
logger_loader_python_field(int dimension, int numpy_type) {
  struct logger_python_field ret;

  ret.typenum = numpy_type;
  ret.dimension = dimension;
  ret.subfields_registred = 0;
  ret.subfields = NULL;

  /* We are done with basic types */
  if (numpy_type != CUSTOM_NPY_TYPE) return ret;

  /* Allocate the memory of the subfields if needed */
  ret.subfields = (struct logger_python_subfield *)malloc(
      dimension * sizeof(struct logger_python_subfield));

  if (ret.subfields == NULL) {
    error_python("Failed to allocate memory for a subfield");
  }

  return ret;
}

/**
 * @brief Define a subtype for a numpy array.
 *
 * @param field The #logger_python_field where to add the subfield.
 * @param dimension The number of dimension for the field.
 * @param type The numpy data type (e.g. NPY_FLOAT32, NPY_DOUBLE, NPY_LONGLONG,
 * ...)
 * @param name The name of the subfield.
 */
__attribute__((always_inline)) INLINE static void
logger_loader_python_field_add_subfield(struct logger_python_field *field,
                                        const char *name, size_t offset,
                                        int dimension, int numpy_type) {

  /* Ensure that we still have enough place */
  if (field->subfields_registred == field->dimension) {
    error_python(
        "Trying to register too many subfields"
        " (currently registering %s)",
        name);
  }

  /* Get the current subfield and register it */
  struct logger_python_subfield *subfield =
      &field->subfields[field->subfields_registred];
  field->subfields_registred += 1;

  /* Set the variables for the subfield */
  subfield->dimension = dimension;
  subfield->typenum = numpy_type;
  subfield->offset = offset;
  strcpy(subfield->name, name);
}

/**
 * @brief Free the memory allocated for the field.
 *
 * @param field The #logger_python_field to clean.
 */
__attribute__((always_inline)) INLINE static void
logger_loader_python_field_free(struct logger_python_field *field) {
  if (field->typenum == CUSTOM_NPY_TYPE) {
    free(field->subfields);
    field->subfields = NULL;
  }
}

/**
 * @brief Generate a python string for the format provided.
 *
 * @param format (output) The format string.
 * @param typenum The C typenum (e.g. NPY_FLOAT32, NPY_DOUBLE, ...).
 * @param dimension The number of dimensions.
 */
__attribute__((always_inline)) INLINE static PyObject *
logger_python_tools_get_format_string(int typenum, int dimension) {

  const char *type_format = NULL;
  switch (typenum) {
    case NPY_FLOAT32:
      type_format = "f4";
      break;
    case NPY_FLOAT64:
      type_format = "f8";
      break;

    case NPY_LONGLONG:
      switch (sizeof(long long)) {
        case 8:
          type_format = "i8";
          break;
        case 4:
          type_format = "i4";
          break;

        default:
          error_python("This size of long long (%li) is not yet supported",
                       sizeof(long long));
      }
      break;

    default:
      error_python("Format not implemented");
  }

  /* Now construct the complete string */
  char format[STRING_SIZE];
  if (dimension == 1) {
    strcpy(format, type_format);
  } else {
    sprintf(format, "%i%s", dimension, type_format);
  }

  /* Convert it into python */
  return PyUnicode_FromString(format);
}

#endif  // HAVE_PYTHON

#endif  // SWIFT_LOGGER_PYTHON_TOOLS_H
