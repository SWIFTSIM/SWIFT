/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Local headers */
#include "error.h"

/* This object's header */
#include "hdf5_object_to_blob.h"

#ifdef HAVE_HDF5
#include <hdf5.h>

#define MAX_FILENAME_LENGTH 500

/**
 * @brief Generates a file name using the process ID and host name
 *
 * @param prefix Prefix to add to the file name
 * @param name Buffer to write the file name to
 * @param len Length of the buffer
 *
 */
static void unique_filename(char *prefix, char *name, size_t len) {

  /* Get our process ID */
  long pid = (long)getpid();

  /* Get our host name. May be truncated or empty. */
  char hostname[MAX_FILENAME_LENGTH];
  if (gethostname(hostname, MAX_FILENAME_LENGTH) != 0) {
    /* gethostname() call failed. Make hostname an empty string. */
    hostname[0] = (char)0;
  } else {
    /* gethostname() successful. Add null terminator in case of truncation . */
    hostname[MAX_FILENAME_LENGTH - 1] = (char)0;
  }

  /* Try to create a unique name. If truncated, try to use it anyway. */
  int ret = snprintf(name, len, "%s.%s.%ld.tmp", prefix, hostname, pid);
  if (ret < 0) {
    error("snprintf call failed generating filename");
  } else if (((size_t)ret) >= len) {
    /* Output has been truncated */
    name[len - 1] = (char)0;
  }
}

/**
 * @brief Copy a HDF5 object and its members to a byte array
 *
 * @param group_id Group containing the object to copy
 * @param name Name of the object
 * @param len Returns the size of the output array
 * @param data Returns a pointer to the output array
 *
 * This creates a HDF5 file image in memory, copies the required
 * object to the file image, then returns a copy of the image as
 * a newly allocated array.
 *
 */
void hdf5_object_to_blob(hid_t group_id, char *name, size_t *len, void **data) {

  /* Set up property list */
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  const size_t increment = 1024 * 1024;
  const hbool_t backing_store = 0;
  if (H5Pset_fapl_core(fapl_id, increment, backing_store) < 0)
    error("Unable to set core driver");

  /* Generate a file name which (hopefully!) doesn't exist */
  char filename[MAX_FILENAME_LENGTH];
  unique_filename("/tmp/__SWIFT_HDF5_CORE_WRITE", filename,
                  MAX_FILENAME_LENGTH);

  /* Create the file image in memory (file with this name must not exist) */
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
  H5Pclose(fapl_id);
  if (file_id < 0) error("Failed to create HDF5 file in memory");

  /* Copy the input object to the file image */
  if (H5Ocopy(group_id, name, file_id, name, H5P_DEFAULT, H5P_DEFAULT) < 0)
    error("Failed to copy HDF5 object %s", name);
  H5Fflush(file_id, H5F_SCOPE_GLOBAL);

  /* Store the size of the file image */
  ssize_t size = H5Fget_file_image(file_id, NULL, 0);
  if (size < 0) error("Failed to get HDF5 file image size");
  *len = (size_t)size;

  /* Copy the file image to a new data array */
  *data = malloc(*len);
  if (H5Fget_file_image(file_id, *data, *len) < 0)
    error("Failed to copy HDF5 file image to array");

  /* Close the file image */
  H5Fclose(file_id);
}

/**
 * @brief Copy a HDF5 object from a byte array to a HDF5 file
 *
 * @param len Size of the input byte array
 * @param data Pointer to the input byte array
 * @param dest_id HDF5 group in which to write the object
 * @param name Name of the object
 *
 * This opens the HDF5 file image supplied in data and copies
 * the named object to the HDF5 group specified by dest_id.
 *
 */
void blob_to_hdf5_object(size_t len, void *data, hid_t dest_id, char *name) {

  /* Set up property list */
  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  const size_t increment = 1024 * 1024;
  const hbool_t backing_store = 0;
  if (H5Pset_fapl_core(fapl_id, increment, backing_store) < 0)
    error("Unable to set core driver");
  if (H5Pset_file_image(fapl_id, data, len) < 0)
    error("Unable to set file image");

  /* Generate a file name which (hopefully!) doesn't exist */
  char filename[MAX_FILENAME_LENGTH];
  unique_filename("/tmp/__SWIFT_HDF5_CORE_READ", filename, MAX_FILENAME_LENGTH);

  /* Open the file image (file with this name must not exist) */
  hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
  if (file_id < 0)
    error(
        "Opening HDF5 file image failed, possibly because file %s already "
        "exists",
        filename);
  H5Pclose(fapl_id);

  /* Copy the object from the file image to the destination */
  if (H5Ocopy(file_id, name, dest_id, name, H5P_DEFAULT, H5P_DEFAULT) < 0)
    error("Failed to copy HDF5 object %s", name);

  /* Close the file image */
  H5Fclose(file_id);
}

#endif
