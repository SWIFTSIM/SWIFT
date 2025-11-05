/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#ifndef SWIFT_SERIAL_IO_H
#define SWIFT_SERIAL_IO_H

/* Config parameters. */
#include <config.h>

#if defined(HAVE_HDF5) && defined(WITH_MPI) && !defined(HAVE_PARALLEL_HDF5)

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Includes. */
#include "ic_info.h"
#include "part.h"

struct engine;
struct unit_system;

void read_ic_serial(char* fileName, const struct unit_system* internal_units,
                    double dim[3], struct part** parts, struct gpart** gparts,
                    struct sink** sinks, struct spart** sparts,
                    struct bpart** bparts, size_t* Ngas, size_t* Ngparts,
                    size_t* Ngparts_background, size_t* Nnuparts,
                    size_t* Nsinks, size_t* Nstars, size_t* Nblackholes,
                    int* flag_entropy, const int with_hydro,
                    const int with_gravity, const int with_sinks,
                    const int with_stars, const int with_black_holes,
                    const int with_cosmology, const int cleanup_h,
                    const int cleanup_sqrt_a, const double h, const double a,
                    const int mpi_rank, int mpi_size, MPI_Comm comm,
                    MPI_Info info, const int n_threads, const int dry_run,
                    const int remap_ids, struct ic_info* ics_metadata);

void write_output_serial(struct engine* e,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units,
                         const int fof, const int mpi_rank, const int mpi_size,
                         MPI_Comm comm, MPI_Info info);

#endif /* HAVE_HDF5 && WITH_MPI && !HAVE_PARALLEL_HDF5 */

#endif /* SWIFT_SERIAL_IO_H */
