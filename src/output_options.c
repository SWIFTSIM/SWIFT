
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Josh Borrow (josh.borrow@durham.ac.uk)
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

#include <stdlib.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include "output_options.h"
#include "parser.h"
#include "swift.h"

/**
 * @brief Initialise the output options struct with the information read
 *        from file. Only rank 0 reads from file; this data is then broadcast
 *        to all nodes.
 *
 * @param parameter_file the pre-parsed parameter file.
 * @param mpi_rank the MPI rank of this node.
 * @param output_options the empty output_options struct (return).
 **/
void output_options_init(struct swift_params* parameter_file, int mpi_rank,
                         struct output_options* output_options) {

  /* Load select_output */
  struct swift_params* select_output =
      (struct swift_params*)malloc(sizeof(struct swift_params));
  if (select_output == NULL)
    error("Error allocating memory for select output options.");

  if (mpi_rank == 0) {
    const int select_output_on = parser_get_opt_param_int(
        parameter_file, "Snapshots:select_output_on", 0);

    if (select_output_on) {
      char select_param_filename[PARSER_MAX_LINE_SIZE];
      parser_get_param_string(parameter_file, "Snapshots:select_output",
                              select_param_filename);
      message("Reading select output parameters from file '%s'",
              select_param_filename);
      parser_read_file(select_param_filename, select_output);
      parser_print_params(select_output);
    }
  }

#ifdef WITH_MPI
  MPI_Bcast(select_output, sizeof(struct swift_params), MPI_BYTE, 0,
            MPI_COMM_WORLD);
#endif

  output_options->select_output = select_output;

  /* Repeat the same process for the compression_options file */
  struct swift_params* compression_options =
      (struct swift_params*)malloc(sizeof(struct swift_params));
  if (compression_options == NULL)
    error("Error allocating memory for compression options parameters.");

  if (mpi_rank == 0) {
    const int compression_options_on = parser_get_opt_param_int(
        parameter_file, "Snapshots:compression_options_on", 0);

    if (compression_options_on) {
      char select_param_filename[PARSER_MAX_LINE_SIZE];
      parser_get_param_string(parameter_file, "Snapshots:compression_options",
                              select_param_filename);
      message("Reading select output parameters from file '%s'",
              select_param_filename);
      parser_read_file(select_param_filename, compression_options);
      parser_print_params(compression_options);
    }
  }

#ifdef WITH_MPI
  MPI_Bcast(compression_options, sizeof(struct swift_params), MPI_BYTE, 0,
            MPI_COMM_WORLD);
#endif

  output_options->compression_options = compression_options;
}

/**
 * @breif Destroys an output_options instance. To be used in tandem
 *        with output_options_init_no_output_list.
 *
 * @param output_options the output_options struct to free the contents of.
 **/
void output_options_clean(struct output_options** output_options) {
  if (*output_options) {
    free((void*)(*output_options)->select_output);
    free((void*)(*output_options)->compression_options);
    *output_options = NULL;
  }
}

/**
 * @brief Dumps the output_options struct to restart file
 *
 * @param output_options pointer to output options struct
 * @param stream bytestream to write to
 **/
void output_options_struct_dump(struct output_options* output_options,
                                FILE* stream) {
  parser_struct_dump(output_options->select_output, stream);
  parser_struct_dump(output_options->compression_options, stream);
}

/**
 * @brief Loads the output_options struct from the restart file
 *
 * @param output_options pointer to the output options struct
 * @param stream bytestream to read from
 **/
void output_options_struct_restore(struct output_options* output_options,
                                   FILE* stream) {
  struct swift_params* select_output =
      (struct swift_params*)malloc(sizeof(struct swift_params));
  parser_struct_restore(select_output, stream);

  struct swift_params* compression_options =
      (struct swift_params*)malloc(sizeof(struct swift_params));
  parser_struct_restore(compression_options, stream);

  output_options->select_output = select_output;
  output_options->compression_options = compression_options;
}