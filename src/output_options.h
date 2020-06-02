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

#include "parser.h"
#include "part_type.h"

/* Compression level names */
enum compression_levels {
  compression_do_not_write = 0,
  compression_write_lossless,
  compression_write_low_lossy,
  compression_write_med_lossy,
  compression_write_high_lossy,
  /* Counter, always leave last */
  compression_level_count,
};

/* Default value for SelectOutput */
#define compression_level_default compression_write_lossless

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

  /* Pass-through struct for now but may need more later. */
};

/* Create and destroy */
void output_options_init(struct swift_params* parameter_file, int mpi_rank,
                         struct output_options* output_options);
void output_options_clean(struct output_options** output_options);

/* Restart files */
void output_options_struct_dump(struct output_options* output_options,
                                FILE* stream);
void output_options_struct_restore(struct output_options* output_options,
                                   FILE* stream);

/* Logic functions */
int output_options_should_write_field(struct output_options* output_options,
                                      char* snapshot_type, char* field_name,
                                      enum part_type part_type);

#endif