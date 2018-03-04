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
#ifndef SWIFT_CHEMISTRY_NONE_H
#define SWIFT_CHEMISTRY_NONE_H

/**
 * @file src/chemistry/none/chemistry.h
 * @brief Empty infrastructure for the cases without chemistry function
 */

/* Some standard headers. */
#include <float.h>
#include <math.h>

/* Local includes. */
#include "chemistry_struct.h"
#include "cosmology.h"
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

  static const char* chemistry_element_names[chemistry_element_count] = {};

  return chemistry_element_names[elem];
}

/**
 * @brief Initialises the chemistry properties.
 *
 * Nothing to do here.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param data The global chemistry information (to be filled).
 */
static INLINE void chemistry_init_backend(
    const struct swift_params* parameter_file, const struct unit_system* us,
    const struct phys_const* phys_const, struct chemistry_data* data) {}

/**
 * @brief Prints the properties of the chemistry model to stdout.
 *
 * @brief The #chemistry_data containing information about the current model.
 */
static INLINE void chemistry_print_backend(const struct chemistry_data* data) {

  message("Chemistry function is 'No chemistry'.");
}

/**
 * @brief Finishes the density calculation.
 *
 * @param p The particle to act upon
 * @param cd The global chemistry information.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void chemistry_end_density(
    struct part* restrict p, const struct chemistry_data* cd,
    const struct cosmology* cosmo) {}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param xp Pointer to the extended particle data.
 * @param data The global chemistry information used for this run.
 */
__attribute__((always_inline)) INLINE static void chemistry_first_init_part(
    const struct part* restrict p, struct xpart* restrict xp,
    const struct chemistry_data* data) {}

/**
 * @brief Sets the chemistry properties of the (x-)particles to a valid start
 * state.
 *
 * Nothing to do here.
 *
 * @param p Pointer to the particle data.
 * @param data The global chemistry information.
 */
__attribute__((always_inline)) INLINE static void chemistry_init_part(
    struct part* restrict p, const struct chemistry_data* data) {}

#endif /* SWIFT_CHEMISTRY_NONE_H */
