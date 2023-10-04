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
#include <config.h>

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "output_options.h"

/* Local headers. */
#include "common_io.h"
#include "error.h"
#include "parser.h"

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

  /* Start by zero-ing everything */
  bzero(output_options, sizeof(struct output_options));

  /* Load select_output */
  struct swift_params* select_output =
      (struct swift_params*)malloc(sizeof(struct swift_params));
  if (select_output == NULL) {
    error("Error allocating memory for select output options.");
  }
  bzero(select_output, sizeof(struct swift_params));

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
    }
  }

#ifdef WITH_MPI
  MPI_Bcast(select_output, sizeof(struct swift_params), MPI_BYTE, 0,
            MPI_COMM_WORLD);
#endif

  output_options->select_output = select_output;
}

/**
 * @brief Destroys an output_options instance.
 *
 * @param output_options the output_options struct to free the contents of.
 **/
void output_options_clean(struct output_options* output_options) {
  if (output_options) {
    free(output_options->select_output);
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

  const size_t count_fields =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * swift_type_count;
  restart_write_blocks(output_options->num_fields_to_write,
                       count_fields * sizeof(int), 1, stream,
                       "output_options_fields", "output options");

  const size_t count_chars =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * FILENAME_BUFFER_SIZE;
  restart_write_blocks(output_options->basenames, count_chars * sizeof(char), 1,
                       stream, "output_options_basenames", "output options");
  restart_write_blocks(output_options->subdir_names, count_chars * sizeof(char),
                       1, stream, "output_options_subdir_names",
                       "output options");

  const size_t count_sub =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * swift_type_count;
  restart_write_blocks(output_options->subsample, count_sub * sizeof(int), 1,
                       stream, "output_options_subsample", "output options");
  restart_write_blocks(output_options->subsample_fractions,
                       count_sub * sizeof(float), 1, stream,
                       "output_options_subsample_fractions", "output options");
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
  bzero(select_output, sizeof(struct swift_params));
  parser_struct_restore(select_output, stream);

  output_options->select_output = select_output;

  const size_t count_fields =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * swift_type_count;
  restart_read_blocks(output_options->num_fields_to_write,
                      count_fields * sizeof(int), 1, stream, NULL,
                      "output options_fields");

  const size_t count_chars =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * FILENAME_BUFFER_SIZE;
  restart_read_blocks(output_options->basenames, count_chars * sizeof(char), 1,
                      stream, NULL, "output_options_basenames");
  restart_read_blocks(output_options->subdir_names, count_chars * sizeof(char),
                      1, stream, NULL, "output_options_subdir_names");

  const size_t count_sub =
      (OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1) * swift_type_count;
  restart_read_blocks(output_options->subsample, count_sub * sizeof(int), 1,
                      stream, NULL, "output_options_subsample");
  restart_read_blocks(output_options->subsample_fractions,
                      count_sub * sizeof(float), 1, stream, NULL,
                      "output_options_subsample_fractions");
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
 * @param ptype integer particle type
 * @param compression_level_current_default The default output strategy
 *        based on the snapshot_type and part_type.
 * @param verbose The verbose level.
 *
 * @return should_write integer determining whether this field should be
 *         written
 **/
enum lossy_compression_schemes output_options_get_field_compression(
    const struct output_options* output_options, const char* snapshot_type,
    const char* field_name, const enum part_type ptype,
    const enum lossy_compression_schemes compression_level_current_default,
    const int verbose) {

  /* Full name for the field path */
  char field[PARSER_MAX_LINE_SIZE];
  sprintf(field, "%.*s:%.*s_%s", FIELD_BUFFER_SIZE, snapshot_type,
          FIELD_BUFFER_SIZE, field_name, part_type_names[ptype]);

  char compression_level[FIELD_BUFFER_SIZE];
  parser_get_opt_param_string(
      output_options->select_output, field, compression_level,
      lossy_compression_schemes_names[compression_level_current_default]);

#ifdef SWIFT_DEBUG_CHECKS

  const int should_write =
      strcmp(lossy_compression_schemes_names[compression_do_not_write],
             compression_level);

  if (verbose) {
    message("Field: %s Write: %s Comp level: \"%s\"", field,
            should_write ? "True" : "False", compression_level);
  }
#endif

  return compression_scheme_from_name(compression_level);
}

/**
 * @brief Return the default output strategy of a given particle type.
 *
 * This can only be "on" or "off". No lossy compression strategy can be
 * applied at the level of an entire particle type.
 *
 * @param output_params The parsed select output file.
 * @param snapshot_type The type of snapshot we are writing
 * @param ptype The #part_type we are considering.
 * @param verbose The verbose level.
 */
enum lossy_compression_schemes output_options_get_ptype_default_compression(
    struct swift_params* output_params, const char* snapshot_type,
    const enum part_type ptype, const int verbose) {

  /* Full name for the default path */
  char field[PARSER_MAX_LINE_SIZE];
  sprintf(field, "%.*s:Standard_%s", FIELD_BUFFER_SIZE, snapshot_type,
          part_type_names[ptype]);

  char compression_level[FIELD_BUFFER_SIZE];
  parser_get_opt_param_string(
      output_params, field, compression_level,
      lossy_compression_schemes_names[compression_level_default]);

  /* Need to find out which of the entries this corresponds to... */
  int level_index;
  for (level_index = 0; level_index < compression_level_count; level_index++) {
    if (!strcmp(lossy_compression_schemes_names[level_index],
                compression_level))
      break;
  }

  /* Make sure that the supplied default option is either on or off, not a
   * compression strategy (these should only be set on a per-field basis) */
  if (!(level_index == compression_do_not_write ||
        level_index == compression_write_lossless)) {
    error(
        "A lossy default compression strategy was specified for snapshot "
        "type %s and particle type %d. This is not allowed, lossy "
        "compression must be set on a field-by-field basis.",
        snapshot_type, ptype);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check whether we could translate the level string to a known entry. */
  if (level_index >= compression_level_count) {
    error(
        "Could not resolve compression level \"%s\" as default compression "
        "level of particle type %s in snapshot type %s.",
        compression_level, part_type_names[ptype], snapshot_type);
  }

  if (verbose) {
    message(
        "Default compression %s, snapshot type %s, compression level %s and "
        "level code %d",
        part_type_names[ptype], snapshot_type, compression_level, level_index);
  }
#endif

  return (enum lossy_compression_schemes)level_index;
}

/**
 * @brief Return the number of fields to be written for a ptype.
 *
 * @param output_options The output_options struct.
 * @param selection_name The current output selection name.
 * @param ptype The particle type index.
 */
int output_options_get_num_fields_to_write(
    const struct output_options* output_options, const char* selection_name,
    const int ptype) {

  /* Get the ID of the output selection in the structure */
  int selection_id =
      parser_get_section_id(output_options->select_output, selection_name);

#ifdef SWIFT_DEBUG_CHECKS
  /* The only situation where we might legitimately not find the selection
   * name is if it is the default. Everything else means trouble. */
  if (strcmp(selection_name, select_output_header_default_name) &&
      selection_id < 0)
    error(
        "Output selection '%s' could not be located in output_options "
        "structure. Please investigate.",
        selection_name);

  /* While we're at it, make sure the selection ID is not impossibly high */
  if (selection_id >= output_options->select_output->sectionCount) {
    error(
        "Output selection '%s' was apparently located in index %d of the "
        "output_options structure, but this only has %d sections.",
        selection_name, selection_id,
        output_options->select_output->sectionCount);
  }
#endif

  /* Special treatment for absent `Default` section */
  if (selection_id < 0) {
    selection_id = output_options->select_output->sectionCount;
  }

  return output_options->num_fields_to_write[selection_id][ptype];
}

/**
 * @brief Return the sub-directory and snapshot basename for the current output
 * selection.
 *
 * @param output_options The #output_options structure
 * @param selection_name The current output selection name.
 * @param default_subdirname The default general sub-directory name.
 * @param default_basename The default general snapshot base name.
 * @param subdir_name (return) The sub-directory name to use for this dump.
 * @param basename (return) The snapshot base name to use for this dump.
 */
void output_options_get_basename(const struct output_options* output_options,
                                 const char* selection_name,
                                 const char* default_subdirname,
                                 const char* default_basename,
                                 char subdir_name[FILENAME_BUFFER_SIZE],
                                 char basename[FILENAME_BUFFER_SIZE]) {

  /* Get the ID of the output selection in the structure */
  int selection_id =
      parser_get_section_id(output_options->select_output, selection_name);

  /* Special treatment for absent `Default` section */
  if (selection_id < 0) {
    selection_id = output_options->select_output->sectionCount;
  }

  /* If the default keyword is found, we use the name provided
   * in the param file (not the select output!), aka. the argument
   * of the function. */
  if (strcmp(output_options->basenames[selection_id],
             select_output_default_basename) == 0) {
    sprintf(basename, "%s", default_basename);
  } else {
    sprintf(basename, "%s", output_options->basenames[selection_id]);
  }

  /* If the default keyword is found, we use the subdir name provided
   * in the param file (not the select output!), aka. the argument
   * of the function. */
  if (strcmp(output_options->subdir_names[selection_id],
             select_output_default_subdir_name) == 0) {
    sprintf(subdir_name, "%s", default_subdirname);
  } else {
    sprintf(subdir_name, "%s", output_options->subdir_names[selection_id]);
  }
}

/**
 * @brief Return the subsampling fraction for the current output selection.
 *
 * @param output_options The #output_options structure.
 * @param selection_name The current output selection name.
 * @param default_subsample The default general subsampling status.
 * @param default_subsample_fraction The default general subsampling fraction.
 * @param subsample (return) The subsampling status to use for this dump.
 * @param subsample_fraction (return) The subsampling fraction to use for this
 * dump.
 */
void output_options_get_subsampling(
    const struct output_options* output_options, const char* selection_name,
    const int default_subsample[swift_type_count],
    const float default_subsample_fraction[swift_type_count],
    int subsample[swift_type_count],
    float subsample_fraction[swift_type_count]) {

  /* Get the ID of the output selection in the structure */
  int selection_id =
      parser_get_section_id(output_options->select_output, selection_name);

  /* Special treatment for absent `Default` section */
  if (selection_id < 0) {
    selection_id = output_options->select_output->sectionCount;
  }

  /* If the default keyword is found, we use the subsample flag provided
   * in the param file (not the select output!), aka. the argument
   * of the function. */
  if (output_options->subsample[selection_id][0] ==
      select_output_default_subsample) {
    memcpy(subsample, default_subsample, sizeof(int) * swift_type_count);
  } else {
    memcpy(subsample, output_options->subsample[selection_id],
           sizeof(int) * swift_type_count);
  }

  /* If the default keyword is found, we use the subsample fraction provided
   * in the param file (not the select output!), aka. the argument
   * of the function. */
  if (output_options->subsample_fractions[selection_id][0] ==
      select_output_default_subsample_fraction) {
    memcpy(subsample_fraction, default_subsample_fraction,
           sizeof(float) * swift_type_count);
  } else {
    memcpy(subsample_fraction,
           output_options->subsample_fractions[selection_id],
           sizeof(float) * swift_type_count);
  }
}
