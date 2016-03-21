/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_PART_H
#define SWIFT_PART_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "const.h"

/* Some constants. */
#define part_align 64
#define gpart_align 32
#define xpart_align 32

/* Import the right particle definition */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_part.h"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_part.h"
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_part.h"
#elif defined(NO_SPH)
#include "./hydro/Gadget2/hydro_part.h"
#else
#error "Invalid choice of SPH variant"
#endif

#if defined(GRAVITY)
#if defined(DEFAULT_GRAVITY)
#include "./gravity/Default/gravity_part.h"
#elif defined(EXTERNAL_POTENTIAL)
#include "./gravity/ExternalPotential/gravity_part.h"
#else
#error "Invalid choice of gravity variant"
#endif
#endif

#ifdef WITH_MPI
void part_create_mpi_type(MPI_Datatype* part_type);
void xpart_create_mpi_type(MPI_Datatype* xpart_type);
#endif

#endif /* SWIFT_PART_H */
