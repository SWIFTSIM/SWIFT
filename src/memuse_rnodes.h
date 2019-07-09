/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#ifndef SWIFT_MEMUSE_RNODE_H
#define SWIFT_MEMUSE_RNODE_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include <stdlib.h>

/* A radix node, this has a single byte key and a pointer to some related
 * resource. It also holds a sorted list of children, if any. */
struct memuse_rnode {

  /* Byte key of this node. */
  uint8_t keypart;

  /* Value of this node, if set. */
  void *ptr;

  /* Sorted pointers to children of this node. */
  struct memuse_rnode **children;
  unsigned int count;
};

void memuse_rnode_dump(int depth, struct memuse_rnode *node, int full);
void memuse_rnode_insert_child(struct memuse_rnode *node, uint8_t depth,
                               uint8_t *key, uint8_t keylen, void *value);
struct memuse_rnode *memuse_rnode_find_child(struct memuse_rnode *node,
                                             uint8_t depth, uint8_t *key,
                                             uint8_t keylen);
void memuse_rnode_cleanup(struct memuse_rnode *node);

#endif /* SWIFT_MEMUSE_RNODE_H */
