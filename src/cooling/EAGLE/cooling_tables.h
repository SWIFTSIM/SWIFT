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

#include "cooling_struct.h"

/*! Number of different bins along the redhsift axis of the tables */
#define eagle_cooling_N_redshifts 49

/*! Number of redshift bins loaded at any given point int time */
#define eagle_cooling_N_loaded_redshifts 2

/*! Number of different bins along the temperature axis of the tables */
#define eagle_cooling_N_temperature 176

/*! Number of different bins along the density axis of the tables */
#define eagle_cooling_N_density 41

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_metal 9

/*! Number of different bins along the metal axis of the tables */
#define eagle_cooling_N_He_frac 7

/*! Number of different bins along the abundances axis of the tables */
#define eagle_cooling_N_abundances 11

void get_cooling_redshifts(struct cooling_function_data *cooling);

void read_cooling_header(const char *fname,
                         struct cooling_function_data *cooling);

void allocate_cooling_tables(struct cooling_function_data *restrict cooling);

void get_redshift_invariant_table(
    struct cooling_function_data *restrict cooling, const int photodis);
void get_cooling_table(struct cooling_function_data *restrict cooling,
                       const int low_z_index, const int high_z_index);

#endif
