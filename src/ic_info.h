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

#ifndef SWIFT_IC_INFO_H
#define SWIFT_IC_INFO_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdlib.h>
#include <stdio.h>

/* Includes. */
#include "parser.h"

struct ic_info {
  char group_name[PARSER_MAX_LINE_SIZE];
  size_t file_image_length;
  void *file_image_data;
};

void ic_info_init(struct ic_info *ics, struct swift_params *params);
void ic_info_clean(struct ic_info *ics);
void ic_info_struct_dump(struct ic_info *ics, FILE *stream);
void ic_info_struct_restore(struct ic_info *ics, FILE *stream);
void ic_info_struct_broadcast(struct ic_info *ics, int root);

#ifdef HAVE_HDF5
void ic_info_read_hdf5(struct ic_info *ics, hid_t file_id);
void ic_info_write_hdf5(struct ic_info *ics, hid_t file_id);
#endif

#endif /* SWIFT_IC_INFO_H */
