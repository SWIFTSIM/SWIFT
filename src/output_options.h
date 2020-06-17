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
#ifndef SWIFT_OUTPUT_OPTIONS_H
#define SWIFT_OUTPUT_OPTIONS_H

/* Local headers. */
#include "output_list.h"
#include "part_type.h"
#include "restart.h"

/**
 * @brief Compression levels for snapshot fields
 */
enum compression_levels {
  compression_do_not_write = 0,
  compression_write_lossless,
  compression_write_low_lossy,
  compression_write_med_lossy,
  compression_write_high_lossy,
  /* Counter, always leave last */
  compression_level_count,
};

/*! Default value for SelectOutput */
#define compression_level_default compression_write_lossless

/*! Default name for the SelectOutput header */
#define select_output_header_default_name "Default"

/**
 * @brief Names of the compression levels, used in the select_output.yml
 *        parameter file.
 **/
extern const char* compression_level_names[];

/**
 * @brief Output selection properties, including the parsed files.
 **/
struct output_options {

  /*! Select output file, parsed */
  struct swift_params* select_output;

  /* Number of fields to write for each output selection and ptype.
   * We need one more than max num of output styles, in case the Default
   * output style is used but not specified. */
  int num_fields_to_write[OUTPUT_LIST_MAX_NUM_OF_SELECT_OUTPUT_STYLES + 1]
                         [swift_type_count];
};

/* Create and destroy */
void output_options_init(struct swift_params* parameter_file, int mpi_rank,
                         struct output_options* output_options);
void output_options_clean(struct output_options* output_options);

/* Restart files */
void output_options_struct_dump(struct output_options* output_options,
                                FILE* stream);
void output_options_struct_restore(struct output_options* output_options,
                                   FILE* stream);

/* Logic functions */
int output_options_should_write_field(
    const struct output_options* output_options, const char* snapshot_type,
    const char* field_name, const enum part_type part_type,
    const enum compression_levels comp_level_current_default);

enum compression_levels output_options_get_ptype_default(
    struct swift_params* output_params, const char* snapshot_type,
    const enum part_type part_type);

int output_options_get_num_fields_to_write(
    const struct output_options* output_options, const char* selection_name,
    const int ptype);

#endif
