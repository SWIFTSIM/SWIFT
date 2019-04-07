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
#ifndef SWIFT_DEBUG_H
#define SWIFT_DEBUG_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include "cell.h"
#include "part.h"
#include "space.h"

void printParticle(const struct part *parts, const struct xpart *xparts,
                   long long int id, size_t N);
void printgParticle(const struct gpart *gparts, const struct part *parts,
                    long long int id, size_t N);
void printParticle_single(const struct part *p, const struct xpart *xp);
void printgParticle_single(struct gpart *gp);

int checkSpacehmax(struct space *s);
int checkCellhdxmax(const struct cell *c, int *depth);
void dumpCells(const char *prefix, int super, int active, int mpiactive,
               int pactive, struct space *s, int rank, int step);

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
#include "metis.h"
void dumpMETISGraph(const char *prefix, idx_t nvtxs, idx_t ncon, idx_t *xadj,
                    idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt);
#endif

#ifdef HAVE_MPI
void dumpCellRanks(const char *prefix, struct cell *cells_top, int nr_cells);
#endif

#endif /* SWIFT_DEBUG_H */
