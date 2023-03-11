/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_HDF5_FUNCTIONS_GEAR_H
#define SWIFT_HDF5_FUNCTIONS_GEAR_H

/* Local includes. */
#include "chemistry.h"
#include "inline.h"

/**
 * @brief Reads a string attribute (array) from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the attribute to read.
 * @param data (output) The attribute read from the HDF5 group (need to be
 * allocated).
 * @param number_element Number of elements in the attribute.
 * @param size_per_element Maximal size per element in data.
 *
 * Calls #error() if an error occurs.
 */
__attribute__((always_inline)) INLINE static void
io_read_string_array_attribute(hid_t grp, const char *name, char *data,
                               hsize_t number_element,
                               hsize_t size_per_element) {

  /* Open attribute */
  const hid_t h_attr = H5Aopen(grp, name, H5P_DEFAULT);
  if (h_attr < 0) error("Error while opening attribute '%s'", name);

  /* Get the number of elements */
  hsize_t count = io_get_number_element_in_attribute(h_attr);

  /* Check if correct number of element */
  if (count != number_element) {
    error(
        "Error found a different number of elements than expected (%lli != "
        "%lli) in attribute %s",
        count, number_element, name);
  }

  /* Get the string length */
  const hid_t type = H5Aget_type(h_attr);
  if (type < 0) error("Failed to get attribute type");

  size_t sdim = H5Tget_size(type);

  /* Check if the size is correct */
  if (sdim > size_per_element) {
    error("Cannot read string longer than %lli in %s", size_per_element, name);
  }

  /* Allocate the temporary array */
  char *tmp = malloc(sizeof(char) * sdim * number_element);
  if (tmp == NULL) {
    error("Failed to allocate the temporary array.");
  }

  /* Read attribute */
  const hid_t h_err = H5Aread(h_attr, type, tmp);
  if (h_err < 0) error("Error while reading attribute '%s'", name);

  /* Copy the attribute correctly */
  for (hsize_t i = 0; i < number_element; i++) {
    char *src = tmp + i * sdim;
    char *dest = data + i * size_per_element;
    memcpy(dest, src, sdim);
  }

  /* Cleanup */
  free(tmp);
  H5Aclose(h_attr);
}

/**
 * @brief Open a group in the yields table (#h5_close_group needs to be called).
 *
 * @param filename The filename.
 * @param group_name The name of the group to open.
 * @param file_id (output) The id of the file opened.
 * @param group_id (output) The id of the group opened.
 *
 */
__attribute__((always_inline)) INLINE static void h5_open_group(
    const char *filename, char *group_name, hid_t *file_id, hid_t *group_id) {

  /* Open file. */
  *file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (*file_id < 0) error("unable to open file %s.\n", filename);

  /* Open group. */
  *group_id = H5Gopen(*file_id, group_name, H5P_DEFAULT);
  if (*group_id < 0) error("unable to open group %s.\n", group_name);
}

/**
 * @brief Close a group in the yields table.
 *
 * @param file_id The id of the file opened.
 * @param group_id The id of the group opened.
 *
 */
__attribute__((always_inline)) INLINE static void h5_close_group(
    hid_t file_id, hid_t group_id) {

  /* Close group */
  hid_t status = H5Gclose(group_id);
  if (status < 0) error("error closing group.");

  /* Close file */
  status = H5Fclose(file_id);
  if (status < 0) error("error closing file.");
}
#endif  // SWIFT_HDF5_FUNCTIONS_GEAR_H
