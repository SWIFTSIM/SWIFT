/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_MPICACHE_H
#define SWIFT_MPICACHE_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "cycle.h"
#include "memuse_rnodes.h"

/* Includes. */
#include <stdlib.h>

/* A cache. */
struct mpicache {
  struct memuse_rnode root;
};

/* API. */
struct mpicache *mpicache_init(int nr_ranks);
void mpicache_add(struct mpicache *cache, int rank, int subtype, int tag,
                  size_t size, void *data);
void mpicache_fetch(struct mpicache *cache, int rank, int subtype, int tag,
                    size_t *size, void **data);
void mpicache_destroy(struct mpicache *cache);

#endif /* SWIFT_MPICACHE_H */
