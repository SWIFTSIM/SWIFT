/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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
#ifndef SWIFT_PARALLEL_IO_H
#define SWIFT_PARALLEL_IO_H

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5) && defined(WITH_MPI) && defined(HAVE_PARALLEL_HDF5)

/* MPI headers. */
#include <mpi.h>

/* Includes. */
#include "part.h"

struct engine;
struct unit_system;

void read_ic_parallel(char* fileName, const struct unit_system* internal_units,
                      double dim[3], struct part** parts, struct gpart** gparts,
                      struct spart** sparts, struct bpart** bparts,
                      size_t* Ngas, size_t* Ngparts, size_t* Ngparts_background,
                      size_t* Nsparts, size_t* Nbparts, int* flag_entropy,
                      int with_hydro, int with_gravity, int with_stars,
                      int with_black_holes, int with_cosmology, int cleanup_h,
                      int cleanup_sqrt_a, double h, double a, int mpi_rank,
                      int mpi_size, MPI_Comm comm, MPI_Info info,
                      int nr_threads, int dry_run);

void write_output_parallel(struct engine* e,
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units,
                           int mpi_rank, int mpi_size, MPI_Comm comm,
                           MPI_Info info);
#endif

#endif /* SWIFT_PARALLEL_IO_H */
