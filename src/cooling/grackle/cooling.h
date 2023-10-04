/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
 * @brief Cooling using the GRACKLE 3.1.1 library.
 */

/* Local includes */
#include "cooling_properties.h"
#include "error.h"
#include "inline.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct pressure_floor_props;
struct space;
struct phys_const;
struct unit_system;
struct swift_params;

/* need to rework (and check) code if changed */
#define GRACKLE_NPART 1
#define GRACKLE_RANK 3

void cooling_update(const struct phys_const* phys_const,
                    const struct cosmology* cosmo,
                    const struct pressure_floor_props* pressure_floor,
                    struct cooling_function_data* cooling, struct space* s,
                    const double time);

void cooling_first_init_part(const struct phys_const* phys_const,
                             const struct unit_system* us,
                             const struct hydro_props* hydro_properties,
                             const struct cosmology* cosmo,
                             const struct cooling_function_data* cooling,
                             const struct part* p, struct xpart* xp);

/**
 * @brief Sets the cooling properties of the (x-)particles to a valid start
 * state. The function requires the density to be defined and thus must
 * be called after its computation.
 *
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param hydro_props The #hydro_props.
 * @param cosmo The #cosmology.
 * @param cooling The properties of the cooling function.
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 */
void cooling_post_init_part(const struct phys_const* phys_const,
                            const struct unit_system* us,
                            const struct hydro_props* hydro_properties,
                            const struct cosmology* cosmo,
                            const struct cooling_function_data* cooling,
                            const struct part* p, struct xpart* xp);

/**
 * @brief Returns the subgrid temperature of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_temperature(const struct part* p,
                                                    const struct xpart* xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

/**
 * @brief Returns the subgrid density of a particle.
 *
 * This model has no subgrid quantity. We return an error.
 *
 * @param p The particle.
 * @param xp The extended particle data.
 */
INLINE static float cooling_get_subgrid_density(const struct part* p,
                                                const struct xpart* xp) {
  error("This cooling model does not use subgrid quantities!");
  return -1.f;
}

float cooling_get_radiated_energy(const struct xpart* restrict xp);
void cooling_print_backend(const struct cooling_function_data* cooling);

void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_properties,
                       const struct entropy_floor_properties* floor_props,
                       const struct pressure_floor_props* pressure_floor,
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

double cooling_get_ycompton(const struct phys_const* phys_const,
                            const struct hydro_props* hydro_props,
                            const struct unit_system* us,
                            const struct cosmology* cosmo,
                            const struct cooling_function_data* cooling,
                            const struct part* p, const struct xpart* xp);

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

/**
 * @brief Compute the electron pressure of a #part based on the cooling
 * function.
 *
 * Does not exist in this model. We return 0.
 *
 * @param phys_const #phys_const data structure.
 * @param hydro_props The properties of the hydro scheme.
 * @param us The internal system of units.
 * @param cosmo #cosmology data structure.
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 * @param xp Pointer to the #xpart data.
 */
INLINE static double cooling_get_electron_pressure(
    const struct phys_const* phys_const, const struct hydro_props* hydro_props,
    const struct unit_system* us, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling, const struct part* p,
    const struct xpart* xp) {
  return 0;
}

#endif /* SWIFT_COOLING_GRACKLE_H */
