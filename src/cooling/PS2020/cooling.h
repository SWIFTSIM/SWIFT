/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_COOLING_PS2020_H
#define SWIFT_COOLING_PS2020_H

/**
 * @file src/cooling/PS2020/cooling.h
 * @brief PS2020 cooling function declarations
 */

/* Local includes. */
#include "cooling_properties.h"
#include "cooling_tables.h"

struct part;
struct xpart;
struct cosmology;
struct hydro_props;
struct entropy_floor_properties;
struct pressure_floor_props;
struct feedback_props;
struct space;

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
                       struct part *p, struct xpart *xp, const float dt,
                       const float dt_therm, const double time);

float cooling_timestep(const struct cooling_function_data *cooling,
                       const struct phys_const *phys_const,
                       const struct cosmology *cosmo,
                       const struct unit_system *us,
                       const struct hydro_props *hydro_props,
                       const struct part *p, const struct xpart *xp);

void cooling_first_init_part(const struct phys_const *phys_const,
                             const struct unit_system *us,
                             const struct hydro_props *hydro_props,
                             const struct cosmology *cosmo,
                             const struct cooling_function_data *cooling,
                             struct part *p, struct xpart *xp);

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
void cooling_post_init_part(const struct phys_const *phys_const,
                            const struct unit_system *us,
                            const struct hydro_props *hydro_props,
                            const struct cosmology *cosmo,
                            const struct cooling_function_data *cooling,
                            struct part *p, struct xpart *xp);

float cooling_get_temperature_from_gas(
    const struct phys_const *phys_const, const struct cosmology *cosmo,
    const struct cooling_function_data *cooling, const float rho_phys,
    const float XH, const float logZZsol, const float u_phys,
    const int HII_region);

float cooling_get_temperature(const struct phys_const *phys_const,
                              const struct hydro_props *hydro_props,
                              const struct unit_system *us,
                              const struct cosmology *cosmo,
                              const struct cooling_function_data *cooling,
                              const struct part *p, const struct xpart *xp);

double cooling_get_electron_density(const struct phys_const *phys_const,
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

float cooling_get_particle_subgrid_density(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

float cooling_get_particle_subgrid_temperature(
    const struct unit_system *us, const struct phys_const *phys_const,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, const struct part *p,
    const struct xpart *xp);

void cooling_set_particle_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, struct part *p,
    struct xpart *xp);

double get_thermal_equilibrium_pressure(
    const struct cooling_function_data *cooling,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const double log_u_cgs, const float log10nH_local, const double rho_cgs,
    const double redshift, const int ired, const int imet, const float dred,
    const float dmet);

double compute_subgrid_property(
    const struct cooling_function_data *cooling,
    const struct phys_const *phys_const,
    const struct entropy_floor_properties *floor_props,
    const struct cosmology *cosmo, const float rho_phys, const float logZZsol,
    const float XH, const float P_phys, const float log10_T,
    const float log10_u_EOS_max_cgs, const int HII_region,
    const float abundance_ratio[colibre_cooling_N_elementtypes],
    const double log_u_cgs, const enum cooling_subgrid_properties isub);

float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp);

float cooling_get_subgrid_density(const struct part *p, const struct xpart *xp);

float cooling_get_radiated_energy(const struct xpart *xp);

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

/**
 * @brief Converts cooling quantities of a particle at the start of a run
 *
 * This function is called once at the end of the engine_init_particle()
 * routine (at the start of a calculation) after the densities of
 * particles have been computed.
 *
 * For this cooling module, this routine does not do anything.
 *
 * @param p The particle to act upon
 * @param xp The extended particle to act upon
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param phys_const #phys_const data structure.
 * @param us Internal system of units data structure.
 * @param floor_props Properties of the entropy floor.
 * @param cooling #cooling_function_data data structure.
 */
__attribute__((always_inline)) INLINE static void cooling_convert_quantities(
    struct part *restrict p, struct xpart *restrict xp,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling) {}

/**
 * @brief Updates cooling properties of particle hit by feedback.
 *
 * For this cooling module, this routine does not do anything.
 *
 * @param xp Pointer to the extended particle data.
 */
__attribute__((always_inline)) INLINE static void
cooling_update_feedback_particle(struct xpart *restrict xp) {}

#endif /* SWIFT_COOLING_PS2020_H */
