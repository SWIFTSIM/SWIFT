/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 James Willis (james.s.willis@durham.ac.uk)
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
#ifndef SWIFT_VELOCIRAPTOR_INTERFACE_H
#define SWIFT_VELOCIRAPTOR_INTERFACE_H

/* Config parameters. */
#include "../config.h"

/* Forward declaration */
struct engine;

/* VELOCIraptor wrapper functions. */
void velociraptor_init(struct engine *e);
void velociraptor_invoke(struct engine *e, const int linked_with_snap);

#endif /* SWIFT_VELOCIRAPTOR_INTERFACE_H */
