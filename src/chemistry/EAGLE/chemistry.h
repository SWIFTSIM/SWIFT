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
#ifndef SWIFT_CHEMISTRY_EAGLE_H
#define SWIFT_CHEMISTRY_EAGLE_H

/**
 * @file src/chemistry/gear/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "error.h"
#include "hydro.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Hydrogen", "Helium",    "Carbon",  "Nitrogen", "Oxygen",
      "Neon",     "Magnesium", "Silicon", "Iron"};

  return chemistry_element_names[elem];
}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    struct part* restrict p, struct xpart* restrict xp,
    const struct chemistry_data* data) {

  p->chemistry_data.metal_mass_fraction_total =
      data->initial_metal_mass_fraction_total;
  for (int elem = 0; elem < chemistry_element_count; ++elem)
    p->chemistry_data.metal_mass_fraction[elem] =
        data->initial_metal_mass_fraction[elem];
}

/**
 * @brief Initialises the chemistry properties.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The properties to initialise.
 */
static INLINE void chemistry_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const, struct chemistry_data* data) {

  /* Read the total metallicity */
  data->initial_metal_mass_fraction_total =
      parser_get_param_float(parameter_file, "EAGLEChemistry:InitMetallicity");

  /* Read the individual mass fractions */
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[50];
    sprintf(buffer, "EAGLEChemistry:InitAbundance_%s",
            chemistry_get_element_name(elem));

    data->initial_metal_mass_fraction[elem] =
        parser_get_param_float(parameter_file, buffer);
  }

  /* Read the constant ratios */
  data->calcium_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:CalciumOverSilicon");
  data->sulphur_over_silicon_ratio = parser_get_param_float(
      parameter_file, "EAGLEChemistry:SulphurOverSilicon");
}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_data containing information about the current model.
 */
static INLINE void chemistry_print_backend(const struct chemistry_data* data) {

  message("Chemistry model is 'EAGLE' tracking %d elements.",
          chemistry_element_count);
}

#endif /* SWIFT_CHEMISTRY_EAGLE_H */
