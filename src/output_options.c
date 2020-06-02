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
#include "../config.h"

/* Some standard headers. */
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "output_options.h"

/* Local headers. */
#include "parser.h"
#include "part_type.h"
#include "swift.h"

/* Compression level names. */
const char* compression_level_names[compression_level_count] = {
    "off", "on", "low", "med", "high"};

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

  output_options->select_output = select_output;
}

/**
 * @brief Decides whether or not a given field should be written. Returns
 *        a truthy value if yes, falsey if not.
 *
 * @param output_options pointer to the output options struct
 * @param snapshot_type pointer to a char array containing the type of
 *        output (i.e. top level section in the yaml).
 * @param field_name pointer to a char array containing the name of the
 *        relevant field.
 * @param part_type integer particle type
 *
 * @return should_write integer determining whether this field should be
 *         written
 **/
int output_options_should_write_field(struct output_options* output_options,
                                      char* snapshot_type, char* field_name,
                                      enum part_type part_type) {
  /* Full name for the field path */
  char field[PARSER_MAX_LINE_SIZE];
  sprintf(field, "%s:%.*s_%s", snapshot_type, FIELD_BUFFER_SIZE, field_name,
          part_type_names[part_type]);

  char compression_level[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(
      output_options->select_output, field, compression_level,
      compression_level_names[compression_level_default]);

  int should_write = strcmp(compression_level_names[compression_do_not_write],
                            compression_level);

#ifdef SWIFT_DEBUG_CHECKS
  message("Determining if %s is to be written. Returning %d from %s.", field,
          should_write, compression_level);
#endif

  return should_write;
}
