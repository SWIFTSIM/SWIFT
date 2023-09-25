/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* This object's header. */
#include "forcing.h"
#include "restart.h"

/**
 * @brief Write an forcing_terms struct to the given FILE as a stream of
 * bytes.
 *
 * @param terms the struct
 * @param stream the file stream
 */
void forcing_terms_struct_dump(const struct forcing_terms* terms,
                               FILE* stream) {
  restart_write_blocks((void*)terms, sizeof(struct forcing_terms), 1, stream,
                       "forcingterms", "forcing terms");
}

/**
 * @brief Restore a forcing_terms struct from the given FILE as a stream of
 * bytes.
 *
 * @param terms the struct
 * @param stream the file stream
 */
void forcing_terms_struct_restore(const struct forcing_terms* terms,
                                  FILE* stream) {
  restart_read_blocks((void*)terms, sizeof(struct forcing_terms), 1, stream,
                      NULL, "forcing terms");
}
