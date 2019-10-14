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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "cooling.h"

/* Local includes */
#include "error.h"
#include "hydro_properties.h"
#include "restart.h"

/**
 * @brief Initialises the cooling properties.
 *
 * Calls cooling_init_backend for the chosen cooling function.
 *
 * @param parameter_file The parsed parameter file.
 * @param us The current internal system of units.
 * @param phys_const The physical constants in internal units.
 * @param hydro_props The properties of the hydro scheme.
 * @param cooling The cooling properties to initialize
 */
void cooling_init(struct swift_params* parameter_file,
                  const struct unit_system* us,
                  const struct phys_const* phys_const,
                  const struct hydro_props* hydro_props,
                  struct cooling_function_data* cooling) {

  /* Verify that we are not doing something stupid here */
  if (hydro_props->minimal_temperature <= 0.)
    error(
        "ERROR: Cannot run with cooling switched on and no minimal "
        "temperature.");

  cooling_init_backend(parameter_file, us, phys_const, hydro_props, cooling);
}

/**
 * @brief Prints the properties of the cooling model to stdout.
 *
 * Calls cooling_print_backend for the chosen cooling function.
 *
 * @param cooling The properties of the cooling function.
 */
void cooling_print(const struct cooling_function_data* cooling) {

  cooling_print_backend(cooling);
}
