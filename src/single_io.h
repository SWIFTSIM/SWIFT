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
#ifndef SWIFT_SINGLE_IO_H
#define SWIFT_SINGLE_IO_H

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Local includes */
#include "part.h"

struct engine;
struct unit_system;

void read_ic_single(const char* fileName,
                    const struct unit_system* internal_units, double dim[3],
                    struct part** parts, struct gpart** gparts,
                    struct sink** sinks, struct spart** sparts,
                    struct bpart** bparts, size_t* Ngas, size_t* Ndm,
                    size_t* Ndm_background, size_t* Nnuparts, size_t* Nsinks,
                    size_t* Nstars, size_t* Nblackholes, int* flag_entropy,
                    int with_hydro, int with_gravity, int with_sinks,
                    int with_stars, int with_black_holes, int with_cosmology,
                    int cleanup_h, int cleanup_sqrt_a, double h, double a,
                    int nr_threads, int dry_run, int remap_ids);

void write_output_single(struct engine* e,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units);

#endif /* HAVE_HDF5 && !WITH_MPI */

#endif /* SWIFT_SINGLE_IO_H */
