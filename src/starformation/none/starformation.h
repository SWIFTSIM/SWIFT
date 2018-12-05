/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
 *******************************************************************************/

#ifndef SWIFT_NO_STARFORMATION_H
#define SWIFT_NO_STARFORMATION_H

/* Some standard headers */
#include <stdlib.h>

/* Local includes */
#include "cosmology.h"
#include "physical_constants.h"
#include "units.h"
#include "parser.h"
#include "equation_of_state.h"

/* Starformation struct */
struct star_formation {};

/*
 * @brief Calculate if the gas has the potential of becoming
 * a star.
 *
 * @param starform the star formation law properties to use.
 * @param p the gas particles
 * @param xp the additional properties of the gas particles
 * @param phys_const the physical constants in internal units
 * @param cosmo the cosmological parameters and properties
 *
 * */
static int starformation_potential_to_become_star(
    const struct star_formation* starform, const struct parts* p,
    const struct xparts* xp, const struct phys_const* const phys_const,
    const struct cosmology* cosmo){

  return 0;
}

/*
 * @brief Calculate if the gas particle is converted 
 *
 * @param starform the star formation struct
 * @param p the gas particles with their properties
 * @param xp the additional gas particle properties
 * @param cosmo the cosmological properties
 *
 * */
static void starformation_convert_to_gas( 
    const struct star_formation* starform, const struct parts* p,
    const struct xparts* xp, const struct cosmology* cosmo) {}

/* 
 * @brief initialization of the star formation law 
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param starform the star formation law properties to initialize
 *
 * */
static void starformation_init_backend(
  struct swift_params* parameter_file, const struct phys_const* phys_const,
  const struct unit_system* us, const struct star_formation* starform) {}

/* @brief Prints the used parameters of the star formation law 
 *
 * @param starform the star formation law properties.
 * */
static void starformation_print_backend(
    const struct star_formation* starform){ 

  message("Star formation law is 'No Star Formation'");

}


#endif /* SWIFT_NO_STARFORMATION_H */
