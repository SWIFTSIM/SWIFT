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

/* Structure that allows the user to easily define the
   fields in a numpy array. */
struct logger_python_field {
  /* Dimension of the field (e.g. 1 for density and 3 for coordinates). */
  int dimension;

  /* Numpy typenum of the array (e.g. NPY_FLOAT32, NPY_INT32). */
  int typenum;
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

  return ret;
}

#endif  // HAVE_PYTHON

#endif  // SWIFT_LOGGER_PYTHON_TOOLS_H
