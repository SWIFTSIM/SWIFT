/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_MPIUSE_H
#define SWIFT_MPIUSE_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "cycle.h"

/* Includes. */
#include <stdlib.h>

/* API. */
#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)
void mpiuse_log_dump(const char *filename, ticks stepticks);
void mpiuse_log_allocation(int type, int subtype, void *ptr, int activation,
                           size_t size, int otherrank, int tag);
#else

/* No-op when not reporting. */
#define mpiuse_log_allocation(type, subtype, ptr, activation, size, otherrank, \
                              tag)                                             \
  ;
#endif /* defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI) */

#endif /* SWIFT_MPIUSE_H */
