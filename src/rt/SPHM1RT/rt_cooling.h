/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c)    2022 Tsang Keung Chan (chantsangkeung@gmail.com)
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
#ifndef SWIFT_RT_SPHM1RT_COOLING_H
#define SWIFT_RT_SPHM1RT_COOLING_H
/**
 * @file src/rt/SPHM1RT/rt_cooling.h
 * @brief Main header file for the SPHM1RT radiative transfer scheme
 * thermochemistry related functions.
 */

#include "rt_properties.h"
#include "rt_unphysical.h"

/**
 * @brief initialize particle quantities relevant for the thermochemistry.
 *
 * @param p part to work with
 * @param rt_props rt_properties struct
 * @param phys_const physical constants struct
 * @param us unit system struct
 * @param cosmo cosmology struct
 */
__attribute__((always_inline)) INLINE static void rt_tchem_first_init_part(
    struct part* restrict p, const struct rt_props* rt_props,
    const struct phys_const* restrict phys_const,
    const struct unit_system* restrict us,
    const struct cosmology* restrict cosmo) {

  struct rt_part_data* rpd = &p->rt_data;

  /* Initialize mass fractions for total metals and each metal individually */
  if (rt_props->initial_metal_mass_fraction_total != -1.f) {
    for (int elem = 0; elem < rt_chemistry_element_count; ++elem) {
      rpd->tchem.metal_mass_fraction[elem] =
          rt_props->initial_metal_mass_fraction[elem];
    }
  }

  /* Initialize species from parameter files */
  if (rt_props->useabundances != 0) {
    for (int spec = 0; spec < rt_species_count; ++spec) {
      rpd->tchem.abundances[spec] = rt_props->initial_species_abundance[spec];
    }
  }

  if (rt_props->skip_thermochemistry != 1) {
    /* Check that we didn't do something stupid */
    rt_check_unphysical_elem_spec(p, rt_props);
  }
}

/**
 * @brief Main function for the thermochemistry step.
 *
 * @param p Particle to work on.
 * @param xp Pointer to the particle' extended data.
 * @param rt_props RT properties struct
 * @param cosmo The current cosmological model.
 * @param hydro_props The #hydro_props.
 * @param phys_const The physical constants in internal units.
 * @param us The internal system of units.
 * @param dt The time-step of this particle.
 */
void rt_do_thermochemistry(struct part* restrict p, struct xpart* restrict xp,
                           struct rt_props* rt_props,
                           const struct cosmology* restrict cosmo,
                           const struct hydro_props* hydro_props,
                           const struct phys_const* restrict phys_const,
                           const struct unit_system* restrict us,
                           const double dt);

#endif /* SWIFT_RT_SPHM1RT_COOLING_H */
