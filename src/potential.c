/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Tom Theuns (tom.theuns@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "potential.h"
#include "restart.h"

/**
 * @brief Initialises the external potential properties in the internal system
 * of units.
 *
 * @param parameter_file The parsed parameter file
 * @param phys_const Physical constants in internal units
 * @param us The current internal system of units
 * @param s The #space we run in.
 * @param potential The external potential properties to initialize
 */
void potential_init(struct swift_params* parameter_file,
                    const struct phys_const* phys_const,
                    const struct unit_system* us, const struct space* s,
                    struct external_potential* potential) {

  potential_init_backend(parameter_file, phys_const, us, s, potential);
}

/**
 * @brief Prints the properties of the external potential to stdout.
 *
 * @param  potential The external potential properties.
 */
void potential_print(const struct external_potential* potential) {

  potential_print_backend(potential);
}

/**
 * @brief Write an external_potential struct to the given FILE as a stream of
 * bytes.
 *
 * @param potential the struct
 * @param stream the file stream
 */
void potential_struct_dump(const struct external_potential* potential,
                           FILE* stream) {
  restart_write_blocks((void*)potential, sizeof(struct external_potential), 1,
                       stream, "externalpotential", "external potential");
}

/**
 * @brief Restore a external_potential struct from the given FILE as a stream of
 * bytes.
 *
 * @param potential the struct
 * @param stream the file stream
 */
void potential_struct_restore(const struct external_potential* potential,
                              FILE* stream) {
  restart_read_blocks((void*)potential, sizeof(struct external_potential), 1,
                      stream, NULL, "external potential");
}
