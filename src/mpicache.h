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
#include "task.h"

/* Includes. */
#include <stdlib.h>

/* A single cache entry. */
struct mpicache_entry {
  int send_node;
  int recv_node;
  struct task *task;
};

/* Shift needed for subtype bits. */
extern int mpicache_subtype_shift;

/* A cache. */
struct mpicache {
  volatile size_t entries_size;
  volatile size_t nr_entries;
  volatile size_t entries_done;
  struct mpicache_entry *volatile entries;

  int *window_sizes;
  int *window_nodes;
  int *window_subtypes;
  int window_size;
  int nr_windows;
};

/* API. */
struct mpicache *mpicache_init(void);
void mpicache_add(struct mpicache *cache, int send_node, int recv_node, struct task *t);
void mpicache_destroy(struct mpicache *cache);
void mpicache_apply(struct mpicache *cache, int sort_by_send);

#endif /* SWIFT_MPICACHE_H */
