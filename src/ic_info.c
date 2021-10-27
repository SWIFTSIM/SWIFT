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
#include <stdlib.h>
#include <stdio.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers */
#include "error.h"

/* This object's header */
#include "ic_info.h"

/* Includes */
#include "hdf5_object_to_blob.h"
#include "restart.h"

extern int engine_rank;

/**
 * @brief Initialize the ic_info struct
 *
 * @param ics Pointer to the struct
 * @param params The Swift parameter set
 */
void ic_info_init(struct ic_info *ics, struct swift_params *params) {
  
  /* Initially we have no file image */
  ics->file_image_length = 0;
  ics->file_image_data = NULL;

  /* Store the name of the HDF5 group to copy */
  parser_get_opt_param_string(params,
                              "InitialConditions:metadata_group_name",
                              ics->group_name, "ICs_parameters");
}

/**
 * @brief Clean up the ic_info struct
 *
 * @param ics Pointer to the struct
 */
void ic_info_clean(struct ic_info *ics) {

  /* Deallocate file image, if it exists */
  if(ics->file_image_data) {
    free(ics->file_image_data);
    ics->file_image_data = NULL;
  }
  ics->file_image_length = 0;

}

/**
 * @brief Dump the ic_info struct to the supplied file
 *
 * @param ics Pointer to the struct
 * @param stream The file to write to
 */
void ic_info_struct_dump(struct ic_info *ics, FILE *stream) {
  
  /* Write the struct */
  restart_write_blocks(ics, sizeof(struct ic_info), 1, stream,
                       "ic_info", "ic_info struct");

  /* Write the HDF5 file image if there is one */
  if(ics->file_image_data) {
    restart_write_blocks(ics->file_image_data, ics->file_image_length,
                         1, stream, "ic_info", "ic_info file image");
  }
}

/**
 * @brief Restore the ic_info struct from the supplied file
 *
 * @param ics Pointer to the struct
 * @param stream The file to read from
 */
void ic_info_struct_restore(struct ic_info *ics, FILE *stream) {

  /* Read in the struct */
  restart_read_blocks(ics, sizeof(struct ic_info), 1, stream,
                      NULL, "ic_info struct");

  /* Read in the HDF5 file image, if there is one */
  if(ics->file_image_data) {
    ics->file_image_data = malloc(ics->file_image_length);
    restart_read_blocks(ics->file_image_data, ics->file_image_length,
                        1, stream, NULL, "ic_info file image");
  }
}

/**
 * @brief Broadcast the ic_info struct to all MPI ranks
 *
 * @param ics Pointer to the struct
 * @param root The root rank for the broadcast operation
 */
void ic_info_struct_broadcast(struct ic_info *ics, int root) {

#ifdef WITH_MPI
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Broadcast the struct */
  MPI_Bcast(ics, sizeof(struct ic_info), MPI_BYTE, root, MPI_COMM_WORLD);
  if(ics->file_image_data) {

    /* Allocate file image data on other ranks */
    if(rank != root)
      ics->file_image_data = malloc(ics->file_image_length);
    
    /* Broadcast the file image data */
    MPI_Bcast(ics->file_image_data, ics->file_image_length,
              MPI_BYTE, root, MPI_COMM_WORLD);
  }
#endif
}


#ifdef HAVE_HDF5
/**
 * @brief Read the ICs metadata from a HDF5 file
 *
 * @param ics Pointer to the struct
 * @param file_id The HDF5 file (or group) to write to
 */
void ic_info_read_hdf5(struct ic_info *ics, hid_t file_id) {
  
  /* Check that the named group exists as a link in the ICs file */
  if(H5Lexists(file_id, ics->group_name, H5P_DEFAULT) > 0) {
    /* Check that the link resolves to an object */
    if(H5Oexists_by_name(file_id, ics->group_name, H5P_DEFAULT) > 0) {      
      /* Read the HDF5 object into memory */
      hdf5_object_to_blob(file_id, ics->group_name, &(ics->file_image_length),
                          &(ics->file_image_data));
      if(engine_rank==0)
        message("Read metadata group %s from ICs file", ics->group_name);
    } else {
      if(engine_rank==0)
        message("Metadata group %s in ICs file appears to be a broken link",
                ics->group_name);
    }
  } else {
    if(engine_rank==0)
      message("Metadata group %s not found in ICs file", ics->group_name);
  }
}


/**
 * @brief Write the ICs metadata to a HDF5 file
 *
 * @param ics Pointer to the struct
 * @param file_id The HDF5 file (or group) to read from
 */
void ic_info_write_hdf5(struct ic_info *ics, hid_t file_id) {
  
  if(ics->file_image_data) {
      blob_to_hdf5_object(ics->file_image_length, ics->file_image_data, file_id,
                          ics->group_name);
    }
}
#endif

