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
#ifndef SWIFT_INTERPOL_EAGLE_H
#define SWIFT_INTERPOL_EAGLE_H

/**
 * @file src/cooling/EAGLE/interpolate.h
 * @brief EAGLE interpolation function declarations
 */

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <hdf5.h>
#include <math.h>
#include <time.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling_struct.h"
#include "error.h"
#include "hydro.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

int row_major_index_2d(int, int, int, int);

int row_major_index_3d(int, int, int, int, int, int); 

int row_major_index_4d(int, int, int, int, int, int, int, int);

void get_index_1d(float *, int, double, int *, float *);

void get_redshift_index(
    float, int *, float *,
    const struct cooling_function_data *restrict);

float interpol_1d(float *, int, float);

double interpol_1d_dbl(double *, int, float);

float interpol_2d(float *, int, int, float, float, int, int, double *, double *);

double interpol_2d_dbl(double *, int, int, double, double, int, int, double *, double *);

float interpol_3d(float *, int, int, int, float, float, float, int, int, int, double *, double *);

float interpol_4d(float *, int, int, int, int, float, float, float, float, int, int, int, int, double *, double *);

void construct_1d_table_from_2d(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, double *, float *, float *);

void construct_1d_table_from_3d(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, float, int, int, double *, float *, float *);

void
construct_1d_print_table_from_3d_elements(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, double *, float *, float *, float *);

void
construct_1d_table_from_3d_elements(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, double *, float *, float *, float *);

void construct_1d_table_from_4d(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float, int,
    int, float, int, int, float, int, int, double *, float *, float *);

void
construct_1d_print_table_from_4d_elements(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, float, int, int, double *, float *, float *, float *);

void
construct_1d_table_from_4d_elements(
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, float *, int, float,
    int, int, float, int, int, double *, float *, float *, float *);

double eagle_convert_temp_to_u_1d_table(double, float *, const struct cooling_function_data *restrict);

double
eagle_convert_u_to_temp(
    double, float *,
    int, int, int, 
    float, float, float,
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *);

double
eagle_convert_u_to_temp_1d_table(
    double, float *, double *,
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *);

void construct_1d_tables(
		int, float, int, float,
		int, float,
                const struct phys_const *restrict,
                const struct unit_system *restrict,
                const struct cosmology *restrict,
                const struct cooling_function_data *restrict,
                const struct part *restrict,
		float *,
		double *,
		double *,
		double *,
		double *,
		double *,
		double *,
		float *, float *);


#endif /* SWIFT_INTERPOL_EAGLE_H */
