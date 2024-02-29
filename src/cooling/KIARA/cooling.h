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
#ifndef SWIFT_COOLING_KIARA_H
#define SWIFT_COOLING_KIARA_H

/**
 * @file src/cooling/KIARA/cooling.h
 * @brief Cooling using the GRACKLE 3.1.1 library.
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
#include "cooling_properties.h"
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
#if COOLING_GRACKLE_MODE >= 2
  #define N_SPECIES 46  /* This further includes properties for dust model */
#else
  #define N_SPECIES 21  /* This includes extra values at end to hold rho,u,dudt,vx,vy,vz,u_floor,mZ,dummyvar */
#endif


void cooling_update(const struct phys_const *phys_const,
                    const struct cosmology *cosmo,
                    const struct pressure_floor_props *pressure_floor,
                    struct cooling_function_data *cooling, struct space *s,
                    const double time);

void cooling_print_fractions(const struct xpart* restrict xp);
void cooling_first_init_part(const struct phys_const* restrict phys_const,
                             const struct unit_system* restrict us,
                             const struct hydro_props* hydro_properties,
                             const struct cosmology* restrict cosmo,
                             const struct cooling_function_data* cooling,
                             struct part* restrict p,
                             struct xpart* restrict xp);
void cooling_post_init_part(const struct phys_const *restrict phys_const,
    			    const struct unit_system *restrict us,
    			    const struct hydro_props *hydro_props,
    			    const struct cosmology *restrict cosmo,
    			    const struct cooling_function_data *restrict cooling,
    			    const struct part *restrict p, struct xpart *restrict xp);

void cooling_print_backend(const struct cooling_function_data* cooling);

void cooling_copy_to_grackle1(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]);
void cooling_copy_to_grackle2(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, 
                              const struct cooling_function_data* restrict cooling,
			      gr_float rho,
                              gr_float species_densities[N_SPECIES]);
void cooling_copy_to_grackle3(grackle_field_data* data, const struct part* p,
                              const struct xpart* xp, gr_float rho,
                              gr_float species_densities[N_SPECIES]);
void cooling_copy_from_grackle1(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho);
void cooling_copy_from_grackle2(grackle_field_data* data, struct part* p,
                                struct xpart* xp,
                                const struct cooling_function_data* restrict cooling,
                                gr_float rho);
void cooling_copy_from_grackle3(grackle_field_data* data, const struct part* p,
                                struct xpart* xp, gr_float rho);
void cooling_copy_to_grackle(grackle_field_data* data,
                             const struct unit_system* restrict us,
                             const struct cosmology* restrict cosmo,
			     const struct cooling_function_data* restrict cooling,
                             const struct part* p, const struct xpart* xp,
			     const double dt, const double u_floor,
			     gr_float species_densities[N_SPECIES],
			     chemistry_data* my_chemistry);
void cooling_copy_from_grackle(grackle_field_data* data, struct part* p,
                               struct xpart* xp, 
			       const struct cooling_function_data* restrict cooling, 
			       gr_float rho);
void cooling_grackle_free_data(grackle_field_data* data);
gr_float cooling_grackle_driver(const struct phys_const* restrict phys_const,
                                const struct unit_system* restrict us,
                                const struct cosmology* restrict cosmo,
                                const struct hydro_props* hydro_properties,
                                const struct cooling_function_data* restrict
                                    cooling,
                                struct part* restrict p,
                                struct xpart* restrict xp, double dt, double u_floor,
				int mode);

gr_float cooling_time(const struct phys_const* restrict phys_const,
                      const struct unit_system* restrict us,
                      const struct hydro_props* hydro_properties,
                      const struct cosmology* restrict cosmo,
                      const struct cooling_function_data* restrict cooling,
                      const struct part* restrict p, struct xpart* restrict xp);

float cooling_get_temperature(
    const struct phys_const* restrict phys_const,
    const struct hydro_props* hydro_properties,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo,
    const struct cooling_function_data* restrict cooling,
    const struct part* restrict p, const struct xpart* restrict xp);

void cooling_cool_part(const struct phys_const* restrict phys_const,
                       const struct unit_system* restrict us,
                       const struct cosmology* restrict cosmo,
                       const struct hydro_props* hydro_properties,
                       const struct entropy_floor_properties* floor_props,
		       const struct pressure_floor_props *pressure_floor_props,
                       const struct cooling_function_data* restrict cooling,
                       struct part* restrict p, struct xpart* restrict xp,
                       const double dt, const double dt_therm,
                       const double time);

void cooling_set_particle_subgrid_properties(
    const struct phys_const *phys_const, const struct unit_system *us,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const struct cooling_function_data *cooling, 
    struct part *p, struct xpart *xp);

float cooling_get_subgrid_temperature(const struct part *p,
                                      const struct xpart *xp);

float cooling_get_subgrid_density(const struct part *p, const struct xpart *xp);

float cooling_get_radiated_energy(const struct xpart* restrict xp);

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

/**
 * @brief Compute the specific thermal energy (physical) for a given temperature.
 *
 * Converts T to u (internal physical units) for a given particle.
 *
 * @param temperature Particle temperature in K
 * @param ne Electron number density relative to H atom density
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 */
INLINE static double cooling_convert_temp_to_u(
    const double temperature, const double ne, const struct cooling_function_data* cooling, 
    const struct part* p) {

  const float X_H = chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
  const float yhelium = (1. - X_H) / (4. * X_H);
  const float mu = (1. + yhelium) / (1. + ne + 4. * yhelium);

  return (temperature * mu * cooling->temp_to_u_factor);
}

/**
 * @brief Compute the temperature for a given physical specific energy
 *
 * Converts T to u (internal physical units) for a given particle.
 *
 * @param u Physical specific energy
 * @param ne Electron number density relative to H atom density
 * @param cooling #cooling_function_data struct.
 * @param p #part data.
 */
INLINE static double cooling_convert_u_to_temp(
    const double u, const double ne, const struct cooling_function_data* cooling, 
    const struct part* p) {

  const float X_H = chemistry_get_metal_mass_fraction_for_cooling(p)[chemistry_element_H];
  const float yhelium = (1. - X_H) / (4. * X_H);
  const float mu = (1. + yhelium) / (1. + ne + 4. * yhelium);

  return u / (mu * cooling->temp_to_u_factor);
}

#endif /* SWIFT_COOLING_KIARA_H */
