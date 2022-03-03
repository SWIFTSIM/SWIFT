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

#ifndef SWIFT_PIXEL_INDEX_H
#define SWIFT_PIXEL_INDEX_H

#include <limits.h>
#include <stdint.h>

/* Config parameters. */
#include "../config.h"

/* Type to use for HEALPix pixel indexes */
typedef int64_t pixel_index_t;

/* Maximum pixel index (determines maximum map size) */
#define MAX_PIXEL_INDEX INT64_MAX

/* Corresponding MPI type */
#define MPI_PIXEL_INDEX_T MPI_INT64_T

#endif
