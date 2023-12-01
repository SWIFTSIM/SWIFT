/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DISTRIBUTED_IO_H
#define SWIFT_DISTRIBUTED_IO_H

/* Config parameters. */
#include <config.h>

#if defined(HAVE_HDF5) && defined(WITH_MPI)

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

struct engine;
struct unit_system;

void write_output_distributed(struct engine* e,
                              const struct unit_system* internal_units,
                              const struct unit_system* snapshot_units,
                              const int fof, int mpi_rank, int mpi_size,
                              MPI_Comm comm, MPI_Info info);

#endif /* HAVE_HDF5 && WITH_MPI */

#endif /* SWIFT_DISTRIBUTED_IO_H */
