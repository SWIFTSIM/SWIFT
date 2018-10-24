/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_EAGLE_COOL_TABLES_H
#define SWIFT_EAGLE_COOL_TABLES_H

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function
 */

/* Config parameters. */
#include "../config.h"

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cooling_struct.h"
#include "error.h"
#include "interpolate.h"

void get_cooling_redshifts(struct cooling_function_data *);

void read_cooling_header(char *, struct cooling_function_data *);

void allocate_cooling_tables(struct cooling_function_data *restrict);

void eagle_check_cooling_tables(struct cooling_function_data *restrict, const int);

#endif
