/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "potentials.h"

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param us The current internal system of units
 * @param potential The external potential properties to initialize
 */
void initPotentialProperties(const struct swift_params * parameter_file,
									  struct UnitSystem* us,
                             struct external_potential* potential) {
  message(" %e\t %e", PARSEC_IN_CGS, units_conversion_factor(us, UNIT_CONV_LENGTH));


  potential->point_mass.x    = parser_get_param_double(parameter_file, "PointMass:position_x");
  potential->point_mass.y    = parser_get_param_double(parameter_file, "PointMass:position_y");
  potential->point_mass.z    = parser_get_param_double(parameter_file, "PointMass:position_z");
  potential->point_mass.mass = parser_get_param_double(parameter_file, "PointMass:mass");
			 
  /* potential->point_mass.x =  */
  /*     50000 * PARSEC_IN_CGS / units_conversion_factor(us, UNIT_CONV_LENGTH); */
  /* potential->point_mass.y = */
  /*     50000 * PARSEC_IN_CGS / units_conversion_factor(us, UNIT_CONV_LENGTH); */
  /* potential->point_mass.z = */
  /*     50000 * PARSEC_IN_CGS / units_conversion_factor(us, UNIT_CONV_LENGTH); */
  /* potential->point_mass.mass = */
  /*     1e10 * SOLAR_MASS_IN_CGS / units_conversion_factor(us, UNIT_CONV_MASS); */
}
