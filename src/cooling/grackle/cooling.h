/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_GRACKLE_H
#define SWIFT_COOLING_GRACKLE_H

/**
 * @file src/cooling/grackle/cooling.h
 * @brief Cooling using the GRACKLE 3.0 library.
 */

/* Some standard headers. */
#include <fenv.h>
#include <float.h>
#include <math.h>

/* The grackle library itself */
#include <grackle.h>

/* Local includes. */
#include "chemistry.h"
#include "cooling_io.h"
#include "entropy_floor.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3

void cooling_update(const struct cosmology* cosmo,
                    struct cooling_function_data* cooling, struct space* s);
void cooling_print_fractions(const struct xpart* restrict xp);
int cooling_converged(const struct xpart* restrict xp,
                      const struct xpart* restrict old, const float limit);
void cooling_compute_equilibrium(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct hydro_props* hydro_properties,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp);
void cooling_first_init_part(const struct phys_const* restrict phys_const,
                             const struct unit_system* restrict us,
                             const struct hydro_props* hydro_properties,
                             const struct cosmology* restrict cosmo,
                             const struct cooling_function_data* cooling,
                             const struct part* restrict p,
                             struct xpart* restrict xp);

float cooling_get_radiated_energy(const struct xpart* restrict xp);
void cooling_print_backend(const struct cooling_function_data* cooling);

void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho);
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho);
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              struct xpart* xp, gr_float rho);
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho);
void cooling_copy_from_grackle2(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho);
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho);
void cooling_copy_to_grackle(grackle_field_data* data, const struct part* p,
                             struct xpart* xp, gr_float rho);
void cooling_copy_from_grackle(grackle_field_data* data, const struct part* p,
                               struct xpart* xp, gr_float rho);
void cooling_apply_self_shielding(
    const struct cooling_function_data* restrict cooling,
    chemistry_data* restrict chemistry, const struct part* restrict p,
    const struct cosmology* cosmo);
gr_float cooling_new_energy(
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct hydro_props* hydro_properties,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, struct xpart* restrict xp, double dt,
    double dt_therm);

gr_float cooling_time(const struct phys_const* restrict phys_const,
                      const struct unit_system* restrict us,
                      const struct hydro_props* hydro_properties,
                      const struct cosmology* restrict cosmo,
                      const struct cooling_function_data* restrict cooling,
                      const struct part* restrict p, struct xpart* restrict xp);

void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_properties,
                       const struct entropy_floor_properties* floor_props,
                       const struct cooling_function_data* restrict cooling,
                       struct part* restrict p, struct xpart* restrict xp,
                       const double dt, const double dt_therm,
                       const double time);

float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* hydro_properties,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp);
float cooling_timestep(const struct cooling_function_data* restrict cooling,
                       const struct phys_const* restrict phys_const,
                       const struct cosmology* restrict cosmo,
                       const struct unit_system* restrict us,
                       const struct hydro_props* hydro_properties,
                       const struct part* restrict p,
                       const struct xpart* restrict xp);

void cooling_split_part(struct part* p, struct xpart* xp, double n);

void cooling_init_units(const struct unit_system* us,
                        const struct phys_const* phys_const,
                        struct cooling_function_data* cooling);
void cooling_init_grackle(struct cooling_function_data* cooling);

void cooling_init_backend(struct swift_params* parameter_file,
                          const struct unit_system* us,
                          const struct phys_const* phys_const,
                          const struct hydro_props* hydro_props,
                          struct cooling_function_data* cooling);

void cooling_clean(struct cooling_function_data* cooling);
void cooling_struct_dump(const struct cooling_function_data* cooling,
                         FILE* stream);
void cooling_struct_restore(struct cooling_function_data* cooling, FILE* stream,
                            const struct cosmology* cosmo);
#endif /* SWIFT_COOLING_GRACKLE_H */
