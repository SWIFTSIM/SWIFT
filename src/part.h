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
#include <config.h>

/* Standard headers. */
#include <stddef.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "align.h"
#include "part_type.h"

/* Pre-declarations */
struct threadpool;

/* Some constants. */
#define part_align 128
#define xpart_align 128
#define spart_align 128
#define gpart_align 128
#define bpart_align 128
#define sink_align 128

/* Import the right hydro particle definition */
#if defined(NONE_SPH)
#include "./hydro/None/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(HOPKINS_PE_SPH)
#include "./hydro/PressureEntropy/hydro_part.h"
#define hydro_need_extra_init_loop 1
#elif defined(HOPKINS_PU_SPH)
#include "./hydro/PressureEnergy/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(HOPKINS_PU_SPH_MONAGHAN)
#include "./hydro/PressureEnergyMorrisMonaghanAV/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(PHANTOM_SPH)
#include "./hydro/Phantom/hydro_part.h"
#define EXTRA_HYDRO_LOOP
#define hydro_need_extra_init_loop 0
#elif defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
#include "./hydro/Gizmo/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#define MPI_SYMMETRIC_FORCE_INTERACTION
#elif defined(SHADOWFAX_SPH)
#include "./hydro/Shadowswift/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#elif defined(PLANETARY_SPH)
#include "./hydro/Planetary/hydro_part.h"
#define hydro_need_extra_init_loop 0
#elif defined(SPHENIX_SPH)
#include "./hydro/SPHENIX/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#elif defined(GASOLINE_SPH)
#include "./hydro/Gasoline/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#elif defined(ANARCHY_PU_SPH)
#include "./hydro/AnarchyPU/hydro_part.h"
#define hydro_need_extra_init_loop 0
#define EXTRA_HYDRO_LOOP
#else
#error "Invalid choice of SPH variant"
#endif

/* Import the right gravity particle definition */
#if defined(DEFAULT_GRAVITY)
#include "./gravity/Default/gravity_part.h"
#elif defined(MULTI_SOFTENING_GRAVITY)
#include "./gravity/MultiSoftening/gravity_part.h"
#else
#error "Invalid choice of gravity variant"
#endif

/* Import the right star particle definition */
#if defined(STARS_NONE)
#include "./stars/None/stars_part.h"
#elif defined(STARS_BASIC)
#include "./stars/Basic/stars_part.h"
#elif defined(STARS_EAGLE)
#include "./stars/EAGLE/stars_part.h"
#elif defined(STARS_GEAR)
#include "./stars/GEAR/stars_part.h"
#else
#error "Invalid choice of star particle"
#endif

/* Import the right black hole particle definition */
#if defined(BLACK_HOLES_NONE)
#include "./black_holes/Default/black_holes_part.h"
#elif defined(BLACK_HOLES_EAGLE)
#include "./black_holes/EAGLE/black_holes_part.h"
#elif defined(BLACK_HOLES_SPIN_JET)
#include "./black_holes/SPIN_JET/black_holes_part.h"
#else
#error "Invalid choice of black hole particle"
#endif

/* Import the right sink particle definition */
#if defined(SINK_NONE)
#include "./sink/Default/sink_part.h"
#elif defined(SINK_GEAR)
#include "./sink/GEAR/sink_part.h"
#else
#error "Invalid choice of sink particle"
#endif

void part_relink_gparts_to_parts(struct part *parts, const size_t N,
                                 const ptrdiff_t offset);
void part_relink_gparts_to_sparts(struct spart *sparts, const size_t N,
                                  const ptrdiff_t offset);
void part_relink_gparts_to_bparts(struct bpart *bparts, const size_t N,
                                  const ptrdiff_t offset);
void part_relink_gparts_to_sinks(struct sink *sinks, const size_t N,
                                 const ptrdiff_t offset);
void part_relink_parts_to_gparts(struct gpart *gparts, const size_t N,
                                 struct part *parts);
void part_relink_sparts_to_gparts(struct gpart *gparts, const size_t N,
                                  struct spart *sparts);
void part_relink_bparts_to_gparts(struct gpart *gparts, const size_t N,
                                  struct bpart *bparts);
void part_relink_sinks_to_gparts(struct gpart *gparts, const size_t N,
                                 struct sink *sinks);
void part_relink_all_parts_to_gparts(struct gpart *gparts, const size_t N,
                                     struct part *parts, struct sink *sinks,
                                     struct spart *sparts, struct bpart *bparts,
                                     struct threadpool *tp);
void part_verify_links(struct part *parts, struct gpart *gparts,
                       struct sink *sinks, struct spart *sparts,
                       struct bpart *bparts, size_t nr_parts, size_t nr_gparts,
                       size_t nr_sinks, size_t nr_sparts, size_t nr_bparts,
                       int verbose);

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
extern MPI_Datatype part_mpi_type;
extern MPI_Datatype xpart_mpi_type;
extern MPI_Datatype gpart_mpi_type;
extern MPI_Datatype spart_mpi_type;
extern MPI_Datatype bpart_mpi_type;

void part_create_mpi_types(void);
void part_free_mpi_types(void);
#endif

#endif /* SWIFT_PART_H */
