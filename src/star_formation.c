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
 ******************************************************************************/

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "part.h"
#include "restart.h"
#include "star_formation.h"
#include "units.h"

/**
 * @brief  Initialises the star formation law properties in the internal
 * unit system.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us the current internal system of units
 * @param hydro_props The propoerties of the hydro scheme.
 * @param starform the properties of the star formation law
 */
void starformation_init(struct swift_params* parameter_file,
                        const struct phys_const* phys_const,
                        const struct unit_system* us,
                        const struct hydro_props* hydro_props,
                        struct star_formation* starform) {

  starformation_init_backend(parameter_file, phys_const, us, hydro_props,
                             starform);
}

/**
 * @brief Print the properties of the star fromation law
 *
 * @param starform the star formation properties.
 */
void starformation_print(const struct star_formation* starform) {

  starformation_print_backend(starform);
}

/**
 * @brief Write an star_formation struct to the given FILE as a stream of
 * bytes.
 *
 * @param starform the star formation struct
 * @param stream the file stream
 */
void starformation_struct_dump(const struct star_formation* starform,
                               FILE* stream) {
  restart_write_blocks((void*)starform, sizeof(struct star_formation), 1,
                       stream, "starformation", "star formation");
}

/**
 * @brief Restore a star_formation struct from the given FILE as a stream of
 * bytes.
 *
 * @param starform the star formation struct
 * @param stream the file stream
 */
void starformation_struct_restore(const struct star_formation* starform,
                                  FILE* stream) {
  restart_read_blocks((void*)starform, sizeof(struct star_formation), 1, stream,
                      NULL, "star formation");
}
