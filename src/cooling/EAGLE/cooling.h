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
#include "hydro_properties.h"
#include "interpolate.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Calculate ratio of particle element abundances
 * to solar abundance. This replaces set_Cooling_SolarAbundances
 * function in EAGLE.
 * Multiple if statements are necessary because order of elements
 * in tables is different from chemistry_element enum.
 * Tables: H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe
 * Enum: H, He, C, N, O, Ne, Mg, Si, Fe
 * The order in ratio_solar is:
 * H, He, C, N, O, Ne, Mg, Si, Fe, S, Ca
 * Hence Fe, S, Ca need to be treated separately to be put in the
 * correct place in the output array.
 *
 * @param p Pointer to #part struct
 * @param cooling #cooling_function_data struct
 * @param ratio_solar Pointer to array or ratios
 */
__attribute__((always_inline)) INLINE void abundance_ratio_to_solar(
    const struct part *restrict p,
    const struct cooling_function_data *restrict cooling, float *ratio_solar) {

  /* compute ratios for all elements */
  for (enum chemistry_element elem = chemistry_element_H;
       elem < chemistry_element_count; elem++) {
    if (elem == chemistry_element_Fe) {
      /* NOTE: solar abundances have iron last with calcium and sulphur directly
       * before, hence +2 */
      ratio_solar[elem] = p->chemistry_data.metal_mass_fraction[elem] /
                          cooling->SolarAbundances[elem + 2];
    } else {
      ratio_solar[elem] = p->chemistry_data.metal_mass_fraction[elem] /
                          cooling->SolarAbundances[elem];
    }
  }

  /* assign ratios for Ca and S, note positions of these elements occur before
   * Fe */
  ratio_solar[chemistry_element_count] =
      p->chemistry_data.metal_mass_fraction[chemistry_element_Si] *
      cooling->sulphur_over_silicon_ratio /
      cooling->SolarAbundances[chemistry_element_count - 1];
  ratio_solar[chemistry_element_count + 1] =
      p->chemistry_data.metal_mass_fraction[chemistry_element_Si] *
      cooling->calcium_over_silicon_ratio /
      cooling->SolarAbundances[chemistry_element_count];
}

void cooling_update(const struct cosmology *, struct cooling_function_data *,
                    const int);

void cooling_cool_part(const struct phys_const *restrict,
                       const struct unit_system *restrict,
                       const struct cosmology *restrict,
                       const struct hydro_props *restrict,
                       const struct cooling_function_data *restrict,
                       struct part *restrict, struct xpart *restrict, float,
                       float);

float cooling_timestep(const struct cooling_function_data *restrict,
                       const struct phys_const *restrict,
                       const struct cosmology *restrict,
                       const struct unit_system *restrict,
                       const struct hydro_props *, const struct part *restrict,
                       const struct xpart *restrict);

void cooling_first_init_part(const struct phys_const *restrict,
                             const struct unit_system *restrict,
                             const struct cosmology *restrict,
                             const struct cooling_function_data *restrict,
                             const struct part *restrict,
                             struct xpart *restrict);

float cooling_get_radiated_energy(const struct xpart *restrict);

void cooling_init_backend(struct swift_params *, const struct unit_system *,
                          const struct phys_const *,
                          struct cooling_function_data *);

void cooling_print_backend(const struct cooling_function_data *);

void cooling_restore_tables(struct cooling_function_data *,
                            const struct cosmology *);

void cooling_clean(struct cooling_function_data *data);

#endif /* SWIFT_COOLING_EAGLE_H */
