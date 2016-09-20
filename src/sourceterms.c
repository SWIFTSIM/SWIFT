/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
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

/* Local includes. */
#include "const.h"
#include "parser.h"
#include "units.h"

/* This object's header. */
#include "sourceterms.h"

/**
 * @brief Initialises the source terms
 *
 * @param parameter_file The parsed parameter file
 * @param us The current internal system of units
 * @param source the structure that has all the source term properties
 */
void source_terms_init(const struct swift_params* parameter_file,
                       struct UnitSystem* us, struct sourceterms* source) {

#ifdef SN_FEEDBACK
  source->supernova.time = parser_get_param_double(parameter_file, "SN:time");
  source->supernova.energy =
      parser_get_param_double(parameter_file, "SN:energy");
  source->supernova.x = parser_get_param_double(parameter_file, "SN:x");
  source->supernova.y = parser_get_param_double(parameter_file, "SN:y");
  source->supernova.z = parser_get_param_double(parameter_file, "SN:z");
  source->supernova.status = supernova_is_not_done;
#endif /* SN_FEEDBCK */
};

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param source the structure that has all the source term properties
 */
void source_terms_print(const struct sourceterms* source) {

#ifdef SN_FEEDBACK
  message(
      " Single SNe of energy= %e will explode at time= %e at location "
      "(%e,%e,%e)",
      source->supernova.energy, source->supernova.time, source->supernova.x,
      source->supernova.y, source->supernova.z);
#endif /* SN_FEEDBACK */
};
