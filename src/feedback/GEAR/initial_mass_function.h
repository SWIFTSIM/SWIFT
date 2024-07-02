/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
#define SWIFT_INITIAL_MASS_FUNCTION_GEAR_H

#include "hdf5_functions.h"
#include "stellar_evolution_struct.h"

float initial_mass_function_get_exponent(
    const struct initial_mass_function *imf, float mass_min, float mass_max);
void initial_mass_function_print(const struct initial_mass_function *imf);
float initial_mass_function_sample(const struct initial_mass_function *imf,
                                   float f);

void initial_mass_function_integrate(const struct initial_mass_function *imf,
                                     float *data, size_t count,
                                     float log_mass_min, float step_size);
float initial_mass_function_get_coefficient(
    const struct initial_mass_function *imf, float mass_min, float mass_max);
float initial_mass_function_get_integral_xi(
    const struct initial_mass_function *imf, float m1, float m2);
float initial_mass_function_get_imf(const struct initial_mass_function *imf,
                                    float m);
float initial_mass_function_get_imf_mass_fraction(
    const struct initial_mass_function *imf, const float m1, const float m2);
float initial_mass_function_get_imf_number_fraction(
    const struct initial_mass_function *imf, const float m1, const float m2);

void initial_mass_function_compute_coefficients(
    struct initial_mass_function *imf);

void initial_mass_function_read_from_table(struct initial_mass_function *imf,
                                           struct swift_params *params,
                                           const char *filename);

void initial_mass_function_init(struct initial_mass_function *imf,
                                const struct phys_const *phys_const,
                                const struct unit_system *us,
                                struct swift_params *params,
                                const char *filename);

void initial_mass_function_dump(const struct initial_mass_function *imf,
                                FILE *stream, const struct stellar_model *sm);

void initial_mass_function_restore(struct initial_mass_function *imf,
                                   FILE *stream,
                                   const struct stellar_model *sm);

void initial_mass_function_clean(struct initial_mass_function *imf);

double initial_mass_function_sample_power_law(double min_mass, double max_mass,
                                              double exp, double x);

void initial_mass_function_compute_Mc_Md_Mtot(
    const struct initial_mass_function *imf, const double minimal_discrete_mass,
    const double stellar_particle_mass, double *M_continuous,
    double *M_discrete, double *M_tot);

#endif  // SWIFT_INITIAL_MASS_FUNCTION_GEAR_H
