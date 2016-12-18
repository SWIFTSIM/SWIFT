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

/* Standard headers. */
#include <stddef.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "align.h"

/* Some constants. */
#define part_align 128
#define xpart_align 128
#define gpart_align 128

/* Import the right hydro particle definition */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(HOPKINS_PE_SPH)
#include "./hydro/PressureEntropy/hydro_part.h"
#define hydro_need_extra_init_loop 1
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(GIZMO_SPH)
#include "./hydro/Gizmo/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#else
#error "Invalid choice of SPH variant"
#endif

/* Import the right gravity particle definition */
#include "./gravity/Default/gravity_part.h"

void part_relink_gparts(struct part *parts, size_t N, ptrdiff_t offset);
void part_relink_parts(struct gpart *gparts, size_t N, struct part *parts);
#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype part_mpi_type;
extern MPI_Datatype xpart_mpi_type;
extern MPI_Datatype gpart_mpi_type;

void part_create_mpi_types();
#endif

#endif /* SWIFT_PART_H */
