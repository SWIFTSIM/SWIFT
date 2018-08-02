/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk).
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
#ifndef SWIFT_VERSION_H
#define SWIFT_VERSION_H

const char* package_description(void);
const char* package_version(void);
const char* hostname(void);
const char* git_revision(void);
const char* git_branch(void);
const char* git_date(void);
const char* configuration_options(void);
const char* compilation_cflags(void);
const char* compiler_name(void);
const char* compiler_version(void);
const char* mpi_version(void);
const char* metis_version(void);
const char* parmetis_version(void);
const char* hdf5_version(void);
const char* fftw3_version(void);
const char* libgsl_version(void);
const char* thread_barrier_version(void);
const char* allocator_version(void);
void greetings(void);

#endif /* SWIFT_VERSION_H */
