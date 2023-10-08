/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_FORCING_H
#define SWIFT_FORCING_H

/**
 * @file src/potential.h
 * @brief Branches between the different external gravitational potentials.
 */

/* Config parameters. */
#include <config.h>

/* Import the right external potential definition */
#if defined(FORCING_NONE)
#include "./forcing/none/forcing.h"
#elif defined(FORCING_ROBERTS_FLOW)
#include "./forcing/roberts_flow/forcing.h"
#elif defined(FORCING_ROBERTS_FLOW_ACCELERATION)
#include "./forcing/roberts_flow_acceleration/forcing.h"
#else
#error "Invalid choice of forcing terms"
#endif

void forcing_terms_struct_dump(const struct forcing_terms* terms, FILE* stream);
void forcing_terms_struct_restore(const struct forcing_terms* terms,
                                  FILE* stream);

#endif /* SWIFT_FORCING_H */
