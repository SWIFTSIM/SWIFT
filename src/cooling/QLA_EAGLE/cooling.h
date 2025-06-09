/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_QLA_EAGLE_H
#define SWIFT_COOLING_QLA_EAGLE_H

/**
 * @file src/cooling/QLA_EAGLE/cooling.h
 * @brief QLA_EAGLE cooling function declarations
 */

/* Local includes. */
#include "cooling_properties.h"
#include "error.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct pressure_floor_props;
struct space;
struct phys_const;

void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time);

void cooling_cool_part(const struct phys_const *phys_const,
                       const struct unit_system *us,
                       const struct cosmology *cosmo,
                       const struct hydro_props *hydro_properties,
                       const struct entropy_floor_properties *floor_props,
                       const struct pressure_floor_props *pressure_floor,
                       const struct cooling_function_data *cooling,
                       struct part *restrict p, struct xpart *restrict xp,
                       const float dt, const float dt_therm, const double time);

float cooling_timestep(const struct cooling_function_data *restrict cooling,
                       const struct phys_const *restrict phys_const,
                       const struct cosmology *restrict cosmo,
                       const struct unit_system *restrict us,
                       const struct hydro_props *hydro_props,
                       const struct part *restrict p,
                       const struct xpart *restrict xp);

float cooling_get_electron_density(const struct phys_const *phys_const,
                                   const struct hydro_props *hydro_props,
                                   const struct unit_system *us,
                                   const struct cosmology *cosmo,
                                   const struct cooling_function_data *cooling,
                                   const struct part *p,
                                   const struct xpart *xp);

double cooling_get_electron_pressure(
    const struct phys_const *phys_const, const struct hydro_props *hydro_props,
    const struct unit_system *us, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double cooling_get_ycompton(const struct phys_const *phys_const,
                            const struct hydro_props *hydro_props,
                            const struct unit_system *us,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            const struct part *p, const struct xpart *xp);

void cooling_first_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp);

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
void cooling_post_init_part(
    const struct phys_const *restrict phys_const,
    const struct unit_system *restrict us,
    const struct hydro_props *hydro_props,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, struct xpart *restrict xp);

float cooling_get_temperature_from_gas(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const float rho_phys,
    const float XH, const float logZZsol, const float u_phys,
    const int HII_region);

float cooling_get_temperature(
    const struct phys_const *restrict phys_const,
    const struct hydro_props *restrict hydro_props,
    const struct unit_system *restrict us,
    const struct cosmology *restrict cosmo,
    const struct cooling_function_data *restrict cooling,
    const struct part *restrict p, const struct xpart *restrict xp);

float cooling_get_particle_subgrid_HI_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

float cooling_get_particle_subgrid_HII_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

float cooling_get_particle_subgrid_H2_fraction(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

double compute_subgrid_property(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const float rho_phys, const float logZZsol,
    const float XH, const float P_phys, const float log10_T,
    const float log10_T_EOS_max, const int HII_region,
    const float *abundance_ratio, const double log_u_cgs,
    const enum cooling_subgrid_properties isub);

float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp);

float cooling_get_subgrid_density(const struct part *p, const struct xpart *xp);

float cooling_get_radiated_energy(const struct xpart *restrict xp);

void cooling_split_part(struct part *p, struct xpart *xp, double n);

void cooling_Hydrogen_reionization(
    const struct cooling_function_data *cooling, const struct cosmology *cosmo,
    const struct pressure_floor_props *pressure_floor, struct space *s);

void cooling_init_backend(struct swift_params *parameter_file,
                          const struct unit_system *us,
                          const struct phys_const *phys_const,
                          const struct hydro_props *hydro_props,
                          struct cooling_function_data *cooling);

void cooling_print_backend(const struct cooling_function_data *cooling);

void cooling_clean(struct cooling_function_data *data);

#endif /* SWIFT_COOLING_QLA_EAGLE_H */
