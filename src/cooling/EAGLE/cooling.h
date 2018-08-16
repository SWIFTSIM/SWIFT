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
#ifndef SWIFT_COOLING_EAGLE_H
#define SWIFT_COOLING_EAGLE_H

/**
 * @file src/cooling/EAGLE/cooling.h
 * @brief EAGLE cooling function declarations
 */

/* Local includes. */
#include "cooling_struct.h"
#include "cosmology.h"
#include "eagle_cool_tables.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

enum hdf5_allowed_types {
  hdf5_short,
  hdf5_int,
  hdf5_long,
  hdf5_float,
  hdf5_double,
  hdf5_char
};

double eagle_helium_reionization_extraheat(double, double, const struct cooling_function_data *restrict);

double eagle_metal_cooling_rate(
    double, double *, int, float, int, float,
    int, float, const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *,
    double *,
    float *);

double
eagle_metal_cooling_rate_1d_table(
    double, double *, double *,
    double *, double *,
    double *, double *,
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *, 
    double *);

double eagle_cooling_rate(
    double, double *, int,
    float, int, float, int, float,
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *,
    float *);

double eagle_cooling_rate_1d_table(
    double, double *, double *,
    double *, double *,
    double *, double *, int,
    float, int, float, int, float,
    const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *,
    float *);

double
eagle_print_metal_cooling_rate(
    int, float, int, float, int,
    float, const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *,
    float *);

double
eagle_print_metal_cooling_rate_1d_table(
    double *, double *,
    double *, double *,
    double *, const struct part *restrict,
    const struct cooling_function_data *restrict,
    const struct cosmology *restrict,
    const struct phys_const *);

float bisection_iter(
    float, double, int,
    float, int, float, int, float, float,
    struct part *restrict, const struct cosmology *restrict,
    const struct cooling_function_data *restrict,
    const struct phys_const *restrict,
    float *, float);

float newton_iter(
    float, double, int,
    float, int, float, int, float,
    float, struct part *restrict, 
    const struct cosmology *restrict,
    const struct cooling_function_data *restrict,
    const struct phys_const *restrict, 
    float *, float, int *);

void abundance_ratio_to_solar(
		const struct part *restrict, 
		const struct cooling_function_data *restrict,
		float *);

void cooling_cool_part(
		const struct phys_const *restrict,
		const struct unit_system *restrict,
		const struct cosmology *restrict,
		const struct cooling_function_data *restrict,
		struct part *restrict, struct xpart *restrict, float);

void cooling_write_flavour(
    hid_t);

float cooling_timestep(
    const struct cooling_function_data *restrict,
    const struct phys_const *restrict,
    const struct cosmology *restrict,
    const struct unit_system *restrict, const struct part *restrict, const struct xpart *restrict); 

void cooling_first_init_part(
    const struct phys_const* restrict,
    const struct unit_system* restrict,
    const struct cosmology* restrict,
    const struct cooling_function_data* restrict,
    const struct part* restrict, struct xpart* restrict);

float cooling_get_radiated_energy(
    const struct xpart *restrict);

void eagle_check_cooling_tables(struct cooling_function_data *, int);

void cooling_update(const struct phys_const*,
                                  const struct unit_system*,
                                  const struct cosmology*,
                                  struct cooling_function_data*);

void cooling_init_backend(
    struct swift_params *, const struct unit_system *,
    const struct phys_const *,
    struct cooling_function_data *);

void cooling_print_backend(const struct cooling_function_data *);

#endif /* SWIFT_COOLING_EAGLE_H */
