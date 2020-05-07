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

/**
 * @brief Output selection properties, including the parsed files.
 **/
struct output_options {

  /*! Select output file, parsed */
  struct swift_params* select_output;

  /*! Compression output file, parsed */
  struct swift_params* compression_options;
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

#endif