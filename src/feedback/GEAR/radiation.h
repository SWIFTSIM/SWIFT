/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_GEAR_H
#define SWIFT_RADIATION_GEAR_H

/**
 * @file src/feedback/GEAR/radiation.h
 * @brief Subgrid radiation feedback for GEAR. This files contains functions to
 * compute quantities for the radiation feedback.
 */

#include "cooling.h"
#include "hydro.h"
#include "part.h"
#include "physical_constants.h"
#include "stellar_evolution_struct.h"
#include "units.h"

#define BB_NU_INTEGRATION_STEPS 500
#define RADIATION_N_IONIZATION_STEPS 1000

double radiation_get_star_ionization_rate(const struct spart* sp);

double radiation_get_part_number_hydrogen_atoms(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp);

double radiation_get_part_rate_to_fully_ionize(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp);

void radiation_tag_part_as_ionized(struct part* p, struct xpart* xpj);
void radiation_reset_part_ionized_tag(struct part* p, struct xpart* xpj);
int radiation_is_part_tagged_as_ionized(struct part* p, struct xpart* xpj);
void radiation_consume_ionizing_photons(struct spart* sp,
                                        double Delta_dot_N_ion);
float radiation_get_comoving_gas_column_density_at_star(const struct spart* sp);

float radiation_get_physical_IR_opacity(const struct spart* sp,
                               const struct unit_system* us,
					const struct phys_const* phys_const,
					const struct cosmology* cosmo);

float radiation_get_physical_IR_optical_depth(const struct spart* sp,
					      const struct unit_system* us,
					      const struct phys_const* phys_const,
					      const struct cosmology* cosmo);

float radiation_get_star_physical_radiation_pressure(const struct spart* sp,
						     const float Delta_t,
						     const struct phys_const* phys_const,
						     const struct unit_system* us,
						     const struct cosmology* cosmo);

int radiation_is_part_ionized(const struct phys_const* phys_const,
                              const struct hydro_props* hydro_props,
                              const struct unit_system* us,
                              const struct cosmology* cosmo,
                              const struct cooling_function_data* cooling,
                              const struct part* p, const struct xpart* xp);

float radiation_get_individual_star_radius(const struct spart* sp,
                                           const struct unit_system* us,
                                           const struct phys_const* phys_const);
float radiation_get_individual_star_temperature(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

float radiation_get_individual_star_luminosity(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

double radiation_get_individual_star_ionizing_photon_emission_rate_fit(
    const struct spart* sp, const struct unit_system* us,
    const struct phys_const* phys_const);

#endif /* SWIFT_RADIATION_GEAR_H */
